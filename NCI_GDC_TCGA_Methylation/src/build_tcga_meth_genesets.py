#!/usr/bin/env python3
"""Master loop for the NCI_GDC_TCGA_Methylation library.

Mirrors the other NCI_GDC_TCGA_* master loops. For each selected tumor type (with matched
normals) x model, dispatch to run_methylation_diff_model.py, which chains the DIG CLI
(workflows methylation_diff_prepare -> convert methylation_cpg_diff -> provenance build).
geneset-extractor-dev stays a thin wrapper; DIG owns the workflow logic.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from meth_selection_io import (
    default_model_list_path, default_model_manifest_path, default_out_root,
    default_tumor_type_list_path, load_model_rows, load_tumor_type_rows,
    model_group_for, relative_or_absolute_path, repo_root, resolve_requested_ids, row_map,
)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Build TCGA DNA methylation tumor-vs-normal gene sets.")
    p.add_argument("--models", default="all")
    p.add_argument("--models_file")
    p.add_argument("--tumor_types", default="all")
    p.add_argument("--tumor_types_file")
    p.add_argument("--model_list", default=str(default_model_list_path()))
    p.add_argument("--tumor_type_list", default=str(default_tumor_type_list_path()))
    p.add_argument("--model_manifest", default=str(default_model_manifest_path()))
    p.add_argument("--python_bin", default=sys.executable or "python3")
    p.add_argument("--beta_matrix_tsv", required=True, help="Probe-by-sample 450K beta matrix (probe_id col + sample cols).")
    p.add_argument("--sample_metadata_tsv", required=True, help="sample_id + project_id + sample_type metadata TSV.")
    p.add_argument("--gtf", required=True, help="GTF for CpG->gene linking + biotype filtering.")
    p.add_argument("--probe_manifest_tsv", help="Optional 450K probe manifest (probe_id, chrom, pos); else DIG auto-resolves.")
    p.add_argument("--dig_dir", required=True)
    p.add_argument("--out_root", default=str(default_out_root()))
    p.add_argument("--overwrite", action="store_true")
    return p


def run_command(command):
    print("$ " + " ".join(command), flush=True)
    subprocess.run(command, check=True)


def dir_nonempty(path: Path) -> bool:
    return path.exists() and path.is_dir() and any(path.iterdir())


def require_existing_file(path_text, label):
    path = relative_or_absolute_path(path_text)
    if not path.exists() or not path.is_file():
        raise SystemExit(f"Missing {label}: {path}")
    return path


def has_solid_tissue_normal(row):
    return str(row.get("has_solid_tissue_normal", "")).strip().lower() in {"true", "1", "yes"}


def main() -> int:
    args = build_parser().parse_args()
    model_rows = load_model_rows(Path(args.model_list))
    tumor_rows = load_tumor_type_rows(Path(args.tumor_type_list))
    selected_models = resolve_requested_ids(csv_text=args.models, file_path=args.models_file, rows=model_rows, key_field="model_id")
    selected_tumor_types = resolve_requested_ids(csv_text=args.tumor_types, file_path=args.tumor_types_file, rows=tumor_rows, key_field="tumor_type_id")
    tumor_by_id = row_map(tumor_rows, "tumor_type_id")

    out_root = Path(args.out_root).resolve()
    outputs_root = out_root / "genesets"
    src_root = repo_root() / "geneset-extractor-dev" / "NCI_GDC_TCGA_Methylation" / "src"

    beta = require_existing_file(args.beta_matrix_tsv, "beta matrix TSV")
    meta = require_existing_file(args.sample_metadata_tsv, "sample metadata TSV")
    gtf = require_existing_file(args.gtf, "GTF")
    manifest = require_existing_file(args.probe_manifest_tsv, "probe manifest") if args.probe_manifest_tsv else None
    dig_dir = relative_or_absolute_path(args.dig_dir).resolve()
    if not dig_dir.is_dir():
        raise SystemExit(f"Missing dig-gene-set-extractors directory: {dig_dir}")
    model_manifest = require_existing_file(args.model_manifest, "model manifest")

    meth_models = [m for m in selected_models if model_group_for(m) == "methylation_diff"]

    conflicts = []
    for tt in selected_tumor_types:
        if not has_solid_tissue_normal(tumor_by_id[tt]):
            continue
        for m in meth_models:
            mo = outputs_root / tt / "models" / m
            if dir_nonempty(mo):
                conflicts.append(str(mo))
    if conflicts and not args.overwrite:
        raise SystemExit("Output already exists (re-run with --overwrite):\n" + "\n".join(conflicts))

    for tt in selected_tumor_types:
        row = tumor_by_id[tt]
        if not has_solid_tissue_normal(row):
            for m in meth_models:
                print(f"  skip {tt}/{m}: project has no matched solid tissue normal", flush=True)
            continue
        project_id = str(row.get("project_id", "")).strip()
        tumor_type_label = str(row.get("tumor_type_label", "")).strip()
        models_root = outputs_root / tt / "models"
        for m in meth_models:
            if args.overwrite and (models_root / m).exists():
                shutil.rmtree(models_root / m)
            cmd = [
                str(Path(args.python_bin).resolve()), str(src_root / "run_methylation_diff_model.py"),
                "--model_id", m, "--tumor_type_id", tt, "--tumor_type_label", tumor_type_label, "--project_id", project_id,
                "--beta_matrix_tsv", str(beta), "--sample_metadata_tsv", str(meta), "--gtf", str(gtf),
                "--run_root", str(models_root), "--python_bin", str(Path(args.python_bin).resolve()),
                "--dig_dir", str(dig_dir), "--model_manifest", str(model_manifest),
            ]
            if manifest is not None:
                cmd += ["--probe_manifest_tsv", str(manifest)]
            run_command(cmd)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
