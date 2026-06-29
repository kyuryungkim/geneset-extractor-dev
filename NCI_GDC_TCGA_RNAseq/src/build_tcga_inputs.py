#!/usr/bin/env python3
"""Build the merged counts matrix + sample metadata for TCGA tumor-vs-rest models.

Consumes open-access GDC STAR-Counts files (``*.rna_seq.augmented_star_gene_counts.tsv``)
plus a GDC sample sheet, and emits:

  - counts_tsv          gene_by_sample matrix: gene_id, gene_symbol, <sample_id...>
  - sample_metadata_tsv sample_id, project_id, sample_type

By default only ``Primary Tumor`` samples are kept (the tumor-vs-rest design pools all
other primary tumors as "rest"). Downloaded inputs live OUTSIDE the committed source tree
(e.g. inputs/NCI_GDC_TCGA_RNAseq/); never hard-code local paths into committed code.

This is a thin data-prep wrapper: the DE/contrast logic lives entirely in DIG
(rna_de_prepare). We only assemble a standard matrix here.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

# GDC STAR count files begin with a header then 4 summary rows (N_unmapped, ...).
SUMMARY_PREFIX = "N_"
COUNT_COLUMN = "unstranded"  # raw counts column used for DE


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts_dir", required=True, help="Directory of GDC STAR-Counts *.tsv files (recursively).")
    p.add_argument("--sample_sheet_tsv", required=True, help="GDC sample sheet TSV.")
    p.add_argument("--out_counts_tsv", required=True)
    p.add_argument("--out_sample_metadata_tsv", required=True)
    p.add_argument("--keep_sample_types", default="Primary Tumor",
                   help="Comma-separated GDC sample_type values to keep (default 'Primary Tumor').")
    p.add_argument("--file_name_column", default="File Name")
    p.add_argument("--sample_id_column", default="Sample ID")
    p.add_argument("--project_id_column", default="Project ID")
    p.add_argument("--sample_type_column", default="Sample Type")
    return p.parse_args()


def read_sample_sheet(path: Path, *, file_col, sample_col, project_col, type_col, keep_types) -> dict[str, dict[str, str]]:
    keep = {t.strip() for t in keep_types.split(",") if t.strip()}
    by_file: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            stype = str(row.get(type_col, "")).strip()
            if keep and stype not in keep:
                continue
            file_name = str(row.get(file_col, "")).strip()
            sample_id = str(row.get(sample_col, "")).strip()
            project_id = str(row.get(project_col, "")).strip()
            if not file_name or not sample_id or not project_id:
                continue
            by_file[file_name] = {"sample_id": sample_id, "project_id": project_id, "sample_type": stype}
    if not by_file:
        raise SystemExit("No samples selected from sample sheet (check --keep_sample_types and column names).")
    return by_file


def read_star_counts(path: Path) -> tuple[list[tuple[str, str]], dict[str, int]]:
    """Return (ordered [(gene_id, gene_symbol)], {gene_id: count})."""
    genes: list[tuple[str, str]] = []
    counts: dict[str, int] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene_id = str(row.get("gene_id", "")).strip()
            if not gene_id or gene_id.startswith(SUMMARY_PREFIX):
                continue
            symbol = str(row.get("gene_name", "")).strip() or gene_id
            try:
                value = int(float(row.get(COUNT_COLUMN, "0") or "0"))
            except ValueError:
                value = 0
            genes.append((gene_id, symbol))
            counts[gene_id] = value
    if not genes:
        raise SystemExit(f"No gene rows parsed from {path} (unexpected format).")
    return genes, counts


def main() -> int:
    args = parse_args()
    counts_dir = Path(args.counts_dir).expanduser().resolve()
    sheet = read_sample_sheet(
        Path(args.sample_sheet_tsv).expanduser().resolve(),
        file_col=args.file_name_column, sample_col=args.sample_id_column,
        project_col=args.project_id_column, type_col=args.sample_type_column,
        keep_types=args.keep_sample_types,
    )

    # Index count files by basename for sample-sheet lookup.
    file_by_name = {p.name: p for p in counts_dir.rglob("*.tsv")}

    gene_order: list[tuple[str, str]] | None = None
    sample_counts: dict[str, dict[str, int]] = {}
    metadata_rows: list[dict[str, str]] = []
    for file_name, meta in sheet.items():
        path = file_by_name.get(file_name)
        if path is None:
            print(f"  WARN: count file not found for {file_name}; skipping", flush=True)
            continue
        genes, counts = read_star_counts(path)
        if gene_order is None:
            gene_order = genes
        sample_counts[meta["sample_id"]] = counts
        metadata_rows.append(meta)

    if gene_order is None or not sample_counts:
        raise SystemExit("No count files matched the sample sheet.")

    sample_ids = [m["sample_id"] for m in metadata_rows]

    out_counts = Path(args.out_counts_tsv).expanduser().resolve()
    out_counts.parent.mkdir(parents=True, exist_ok=True)
    with out_counts.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "gene_symbol", *sample_ids])
        for gene_id, symbol in gene_order:
            writer.writerow([gene_id, symbol, *[sample_counts[s].get(gene_id, 0) for s in sample_ids]])

    out_meta = Path(args.out_sample_metadata_tsv).expanduser().resolve()
    out_meta.parent.mkdir(parents=True, exist_ok=True)
    with out_meta.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["sample_id", "project_id", "sample_type"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(metadata_rows)

    print(f"Wrote {len(sample_ids)} samples x {len(gene_order)} genes -> {out_counts}", flush=True)
    print(f"Wrote sample metadata -> {out_meta}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
