#!/usr/bin/env python3
"""Build the merged 450K beta matrix + sample metadata for TCGA methylation models.

Consumes open-access GDC "Methylation Beta Value" files (SeSAMe; 2-column:
Composite Element + Beta Value) plus a GDC sample sheet, and emits:
  - beta_matrix_tsv   probe-by-sample: probe_id + one column per sample_id
  - sample_metadata_tsv  sample_id, project_id, sample_type

Keeps Primary Tumor + Solid Tissue Normal by default (tumor-vs-normal design). Downloaded
inputs live OUTSIDE the committed tree (e.g. inputs/NCI_GDC_TCGA_Methylation/). Thin data-prep;
DIG owns the differential workflow (methylation_diff_prepare).
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--beta_dir", required=True, help="Directory of GDC Methylation Beta Value *.txt files (recursive).")
    p.add_argument("--sample_sheet_tsv", required=True)
    p.add_argument("--out_beta_matrix_tsv", required=True)
    p.add_argument("--out_sample_metadata_tsv", required=True)
    p.add_argument("--keep_sample_types", default="Primary Tumor,Solid Tissue Normal")
    p.add_argument("--file_name_column", default="File Name")
    p.add_argument("--sample_id_column", default="Sample ID")
    p.add_argument("--project_id_column", default="Project ID")
    p.add_argument("--sample_type_column", default="Sample Type")
    return p.parse_args()


def read_sample_sheet(path, *, file_col, sample_col, project_col, type_col, keep_types):
    keep = {t.strip() for t in keep_types.split(",") if t.strip()}
    by_file = {}
    with open(path, "r", encoding="utf-8", newline="") as h:
        for row in csv.DictReader(h, delimiter="\t"):
            st = str(row.get(type_col, "")).strip()
            if keep and st not in keep:
                continue
            fn = str(row.get(file_col, "")).strip(); sid = str(row.get(sample_col, "")).strip(); pr = str(row.get(project_col, "")).strip()
            if fn and sid and pr:
                by_file[fn] = {"sample_id": sid, "project_id": pr, "sample_type": st}
    if not by_file:
        raise SystemExit("No samples selected from sample sheet.")
    return by_file


def read_beta_file(path):
    """Return {probe_id: beta_str} from a 2-column GDC SeSAMe beta file (skips header if present)."""
    out = {}
    with open(path, "r", encoding="utf-8", newline="") as h:
        for i, line in enumerate(h):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            probe, beta = parts[0].strip(), parts[1].strip()
            if i == 0 and (probe.lower() in {"composite element ref", "composite element", "probe_id"} or beta.lower() in {"beta value", "beta_value"}):
                continue
            if probe:
                out[probe] = beta
    return out


def main() -> int:
    args = parse_args()
    beta_dir = Path(args.beta_dir).expanduser().resolve()
    sheet = read_sample_sheet(Path(args.sample_sheet_tsv).expanduser().resolve(),
                              file_col=args.file_name_column, sample_col=args.sample_id_column,
                              project_col=args.project_id_column, type_col=args.sample_type_column,
                              keep_types=args.keep_sample_types)
    file_by_name = {p.name: p for p in beta_dir.rglob("*.txt")}

    sample_betas = {}
    probe_order = None
    meta_rows = []
    for fn, meta in sheet.items():
        path = file_by_name.get(fn)
        if path is None:
            print(f"  WARN: beta file not found for {fn}; skipping", flush=True)
            continue
        betas = read_beta_file(path)
        if probe_order is None:
            probe_order = list(betas.keys())
        sample_betas[meta["sample_id"]] = betas
        meta_rows.append(meta)
    if probe_order is None or not sample_betas:
        raise SystemExit("No beta files matched the sample sheet.")

    sample_ids = [m["sample_id"] for m in meta_rows]
    out_beta = Path(args.out_beta_matrix_tsv).expanduser().resolve(); out_beta.parent.mkdir(parents=True, exist_ok=True)
    with out_beta.open("w", encoding="utf-8", newline="") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(["probe_id", *sample_ids])
        for probe in probe_order:
            w.writerow([probe, *[sample_betas[s].get(probe, "") for s in sample_ids]])
    out_meta = Path(args.out_sample_metadata_tsv).expanduser().resolve(); out_meta.parent.mkdir(parents=True, exist_ok=True)
    with out_meta.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, delimiter="\t", fieldnames=["sample_id", "project_id", "sample_type"], lineterminator="\n")
        w.writeheader(); w.writerows(meta_rows)
    print(f"Wrote {len(sample_ids)} samples x {len(probe_order)} probes -> {out_beta}", flush=True)
    print(f"Wrote sample metadata -> {out_meta}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
