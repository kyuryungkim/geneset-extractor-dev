#!/usr/bin/env python3
"""Build the merged cohort segments TSV + sample metadata for TCGA CNV models.

Consumes open-access GDC "Masked Copy Number Segment" files (DNAcopy;
``*.nocnv_grch38.seg.v2.txt`` with columns GDC_Aliquot, Chromosome, Start, End,
Num_Probes, Segment_Mean) plus a GDC sample sheet, and emits:

  - segments_tsv         merged GDC-style seg: sample_id, Chromosome, Start, End, Segment_Mean
  - sample_metadata_tsv  sample_id, project_id, sample_type

By default keeps tumor sample types (Primary Tumor, Metastatic). cnv_gene_extractor
maps segments -> genes via the GTF; downloaded inputs live OUTSIDE the committed tree
(e.g. inputs/NCI_GDC_TCGA_CNV/). This is a thin data-prep wrapper; DIG owns the workflow.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--segments_dir", required=True, help="Directory of GDC Masked Copy Number Segment *.txt files (recursive).")
    p.add_argument("--sample_sheet_tsv", required=True, help="GDC sample sheet TSV.")
    p.add_argument("--out_segments_tsv", required=True)
    p.add_argument("--out_sample_metadata_tsv", required=True)
    p.add_argument("--keep_sample_types", default="Primary Tumor,Metastatic")
    p.add_argument("--file_name_column", default="File Name")
    p.add_argument("--sample_id_column", default="Sample ID")
    p.add_argument("--project_id_column", default="Project ID")
    p.add_argument("--sample_type_column", default="Sample Type")
    # GDC seg column names (case-insensitive resolved below)
    p.add_argument("--seg_chrom_column", default="Chromosome")
    p.add_argument("--seg_start_column", default="Start")
    p.add_argument("--seg_end_column", default="End")
    p.add_argument("--seg_mean_column", default="Segment_Mean")
    return p.parse_args()


def read_sample_sheet(path, *, file_col, sample_col, project_col, type_col, keep_types):
    keep = {t.strip() for t in keep_types.split(",") if t.strip()}
    by_file = {}
    with open(path, "r", encoding="utf-8", newline="") as h:
        for row in csv.DictReader(h, delimiter="\t"):
            stype = str(row.get(type_col, "")).strip()
            if keep and stype not in keep:
                continue
            fn = str(row.get(file_col, "")).strip()
            sid = str(row.get(sample_col, "")).strip()
            proj = str(row.get(project_col, "")).strip()
            if fn and sid and proj:
                by_file[fn] = {"sample_id": sid, "project_id": proj, "sample_type": stype}
    if not by_file:
        raise SystemExit("No samples selected from sample sheet (check --keep_sample_types / column names).")
    return by_file


def resolve_col(fieldnames, want):
    low = {f.lower(): f for f in fieldnames}
    return low.get(want.lower())


def read_seg_file(path, *, chrom_c, start_c, end_c, mean_c):
    with open(path, "r", encoding="utf-8", newline="") as h:
        reader = csv.DictReader(h, delimiter="\t")
        fn = reader.fieldnames or []
        c = {k: resolve_col(fn, v) for k, v in {"chrom": chrom_c, "start": start_c, "end": end_c, "mean": mean_c}.items()}
        if not all(c.values()):
            raise SystemExit(f"{path}: missing expected seg columns (have {fn})")
        out = []
        for row in reader:
            out.append((row[c["chrom"]], row[c["start"]], row[c["end"]], row[c["mean"]]))
        return out


def main() -> int:
    args = parse_args()
    segments_dir = Path(args.segments_dir).expanduser().resolve()
    sheet = read_sample_sheet(
        Path(args.sample_sheet_tsv).expanduser().resolve(),
        file_col=args.file_name_column, sample_col=args.sample_id_column,
        project_col=args.project_id_column, type_col=args.sample_type_column,
        keep_types=args.keep_sample_types,
    )
    file_by_name = {p.name: p for p in segments_dir.rglob("*.txt")}

    out_seg = Path(args.out_segments_tsv).expanduser().resolve()
    out_seg.parent.mkdir(parents=True, exist_ok=True)
    meta_rows = []
    n_seg = 0
    with out_seg.open("w", encoding="utf-8", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(["sample_id", "Chromosome", "Start", "End", "Segment_Mean"])
        for fn, meta in sheet.items():
            path = file_by_name.get(fn)
            if path is None:
                print(f"  WARN: seg file not found for {fn}; skipping", flush=True)
                continue
            rows = read_seg_file(path, chrom_c=args.seg_chrom_column, start_c=args.seg_start_column,
                                 end_c=args.seg_end_column, mean_c=args.seg_mean_column)
            for chrom, start, end, mean in rows:
                w.writerow([meta["sample_id"], chrom, start, end, mean]); n_seg += 1
            meta_rows.append(meta)
    if not meta_rows:
        raise SystemExit("No seg files matched the sample sheet.")

    out_meta = Path(args.out_sample_metadata_tsv).expanduser().resolve()
    out_meta.parent.mkdir(parents=True, exist_ok=True)
    with out_meta.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, delimiter="\t", fieldnames=["sample_id", "project_id", "sample_type"], lineterminator="\n")
        w.writeheader(); w.writerows(meta_rows)

    print(f"Wrote {len(meta_rows)} samples / {n_seg} segments -> {out_seg}", flush=True)
    print(f"Wrote sample metadata -> {out_meta}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
