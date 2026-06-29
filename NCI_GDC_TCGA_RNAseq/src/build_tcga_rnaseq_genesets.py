#!/usr/bin/env python3
"""Master loop for the NCI_GDC_TCGA_RNAseq library.

Mirrors GTEx/src/build_genesets.py. For each selected tumor type x model, dispatch
to run_tumor_vs_rest_model.py, which calls the DIG CLI
(workflows rna_de_prepare --comparison_mode group_vs_rest) then
(convert rna_deg_multi). geneset-extractor-dev stays a thin wrapper; DIG owns the
workflow logic.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from tcga_rnaseq_selection_io import (
    default_model_list_path,
    default_model_manifest_path,
    default_out_root,
    default_tumor_type_list_path,
    load_model_rows,
    load_tumor_type_rows,
    model_group_for,
    relative_or_absolute_path,
    repo_root,
    resolve_requested_ids,
    row_map,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build TCGA RNA-seq tumor-type gene sets.")
    parser.add_argument("--models", default="all")
    parser.add_argument("--models_file")
    parser.add_argument("--tumor_types", default="all")
    parser.add_argument("--tumor_types_file")
    parser.add_argument("--model_list", default=str(default_model_list_path()))
    parser.add_argument("--tumor_type_list", default=str(default_tumor_type_list_path()))
    parser.add_argument("--model_manifest", default=str(default_model_manifest_path()))
    parser.add_argument("--python_bin", default=sys.executable or "python3")
    parser.add_argument(
        "--counts_tsv",
        required=True,
        help="Merged gene_by_sample counts TSV across all selected projects (see build_tcga_inputs.py).",
    )
    parser.add_argument(
        "--sample_metadata_tsv",
        required=True,
        help="Sample metadata TSV with sample_id + project_id (group column) columns.",
    )
    parser.add_argument("--gtf", help="GTF for biotype filtering (required if any model sets require_gtf).")
    parser.add_argument("--provenance_mirror_local_prefix")
    parser.add_argument("--provenance_mirror_remote_prefix")
    parser.add_argument("--dig_dir", required=True)
    parser.add_argument("--out_root", default=str(default_out_root()))
    parser.add_argument("--overwrite", action="store_true")
    return parser


def run_command(command: list[str]) -> None:
    print("$ " + " ".join(command), flush=True)
    subprocess.run(command, check=True)


def dir_nonempty(path: Path) -> bool:
    return path.exists() and path.is_dir() and any(path.iterdir())


def overwrite_dir(path: Path) -> None:
    if path.exists():
        shutil.rmtree(path)


def model_requires_gtf(row: dict[str, str]) -> bool:
    return str(row.get("require_gtf", "")).strip().lower() in {"true", "1", "yes"}


def require_existing_file(path_text: str, label: str) -> Path:
    path = relative_or_absolute_path(path_text)
    if not path.exists() or not path.is_file():
        raise SystemExit(f"Missing {label}: {path}")
    return path


def main() -> int:
    args = build_parser().parse_args()
    model_rows = load_model_rows(Path(args.model_list))
    tumor_rows = load_tumor_type_rows(Path(args.tumor_type_list))
    selected_models = resolve_requested_ids(
        csv_text=args.models, file_path=args.models_file, rows=model_rows, key_field="model_id"
    )
    selected_tumor_types = resolve_requested_ids(
        csv_text=args.tumor_types, file_path=args.tumor_types_file, rows=tumor_rows, key_field="tumor_type_id"
    )
    tumor_by_id = row_map(tumor_rows, "tumor_type_id")
    model_by_id = row_map(model_rows, "model_id")

    out_root = Path(args.out_root).resolve()
    outputs_root = out_root / "genesets"
    src_root = repo_root() / "geneset-extractor-dev" / "NCI_GDC_TCGA_RNAseq" / "src"

    counts_tsv = require_existing_file(args.counts_tsv, "counts TSV")
    sample_metadata_tsv = require_existing_file(args.sample_metadata_tsv, "sample metadata TSV")

    dig_dir = relative_or_absolute_path(args.dig_dir).resolve()
    if not dig_dir.is_dir():
        raise SystemExit(f"Missing dig-gene-set-extractors directory: {dig_dir}")

    # Only tumor_vs_rest is implemented for the first model.
    tumor_vs_rest_models = [m for m in selected_models if model_group_for(m) == "tumor_vs_rest"]
    unsupported = [m for m in selected_models if model_group_for(m) != "tumor_vs_rest"]
    if unsupported:
        raise SystemExit(f"Models not implemented yet: {', '.join(unsupported)}")

    # gtf requirement check, library-wide (mirror GTEx).
    model_list_gtf_required = [str(r["model_id"]).strip() for r in model_rows if model_requires_gtf(r)]
    resolved_gtf: Path | None = require_existing_file(args.gtf, "GTF") if args.gtf else None
    if model_list_gtf_required and resolved_gtf is None:
        raise SystemExit(
            "model_list contains models requiring --gtf but none was provided: "
            + ", ".join(model_list_gtf_required)
        )

    model_manifest = require_existing_file(args.model_manifest, "model manifest")

    # Pre-flight conflict check.
    conflicts: list[str] = []
    for tumor_type_id in selected_tumor_types:
        for model_id in tumor_vs_rest_models:
            model_out = outputs_root / tumor_type_id / "models" / model_id
            if dir_nonempty(model_out):
                conflicts.append(str(model_out))
    if conflicts and not args.overwrite:
        raise SystemExit(
            "Output already exists (re-run with --overwrite):\n" + "\n".join(conflicts)
        )

    for tumor_type_id in selected_tumor_types:
        tumor_row = tumor_by_id[tumor_type_id]
        project_id = str(tumor_row.get("project_id", "")).strip()
        tumor_type_label = str(tumor_row.get("tumor_type_label", "")).strip()
        if not project_id or not tumor_type_label:
            raise SystemExit(f"Missing project_id/tumor_type_label for {tumor_type_id}")
        models_root = outputs_root / tumor_type_id / "models"
        for model_id in tumor_vs_rest_models:
            if args.overwrite:
                overwrite_dir(models_root / model_id)
            needs_gtf = model_requires_gtf(model_by_id[model_id])
            cmd = [
                str(Path(args.python_bin).resolve()),
                str(src_root / "run_tumor_vs_rest_model.py"),
                "--model_id", model_id,
                "--tumor_type_id", tumor_type_id,
                "--tumor_type_label", tumor_type_label,
                "--project_id", project_id,
                "--counts_tsv", str(counts_tsv),
                "--sample_metadata_tsv", str(sample_metadata_tsv),
                "--group_column", "project_id",
                "--run_root", str(models_root),
                "--python_bin", str(Path(args.python_bin).resolve()),
                "--dig_dir", str(dig_dir),
                "--model_manifest", str(model_manifest),
            ]
            if resolved_gtf is not None and needs_gtf:
                cmd += ["--gtf", str(resolved_gtf)]
            if args.provenance_mirror_local_prefix:
                cmd += ["--provenance_mirror_local_prefix", args.provenance_mirror_local_prefix]
            if args.provenance_mirror_remote_prefix:
                cmd += ["--provenance_mirror_remote_prefix", args.provenance_mirror_remote_prefix]
            run_command(cmd)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
