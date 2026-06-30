#!/usr/bin/env python3
"""Master loop for NCI_GDC_TCGA_ATAC. For each cancer type x model dispatch to
run_atac_accessibility_model.py (convert atac_bulk_matrix, cancer-type-vs-rest).
Thin wrapper; DIG owns the workflow. No DIG prepare step (self-contained converter)."""
from __future__ import annotations
import argparse, subprocess, shutil, sys
from pathlib import Path
from atac_selection_io import (default_model_list_path, default_model_manifest_path, default_out_root,
    default_tumor_type_list_path, load_model_rows, load_tumor_type_rows, model_group_for,
    relative_or_absolute_path, repo_root, resolve_requested_ids, row_map)

def build_parser():
    p=argparse.ArgumentParser(description="Build TCGA ATAC accessibility gene sets.")
    p.add_argument("--models", default="all"); p.add_argument("--models_file")
    p.add_argument("--tumor_types", default="all"); p.add_argument("--tumor_types_file")
    p.add_argument("--model_list", default=str(default_model_list_path()))
    p.add_argument("--tumor_type_list", default=str(default_tumor_type_list_path()))
    p.add_argument("--model_manifest", default=str(default_model_manifest_path()))
    p.add_argument("--python_bin", default=sys.executable or "python3")
    p.add_argument("--peak_matrix_tsv", required=True, help="Pan-cancer peak-by-sample accessibility matrix.")
    p.add_argument("--peak_bed", required=True, help="Peak coordinates BED (row-aligned to matrix).")
    p.add_argument("--sample_metadata_tsv", required=True, help="sample_id + cancer_type metadata TSV.")
    p.add_argument("--gtf", required=True)
    p.add_argument("--dig_dir", required=True); p.add_argument("--out_root", default=str(default_out_root()))
    p.add_argument("--overwrite", action="store_true")
    return p

def run_command(c): print("$ "+" ".join(c), flush=True); subprocess.run(c, check=True)
def dir_nonempty(p): return p.exists() and p.is_dir() and any(p.iterdir())
def require_file(t,l):
    p=relative_or_absolute_path(t)
    if not p.exists() or not p.is_file(): raise SystemExit(f"Missing {l}: {p}")
    return p

def main()->int:
    a=build_parser().parse_args()
    model_rows=load_model_rows(Path(a.model_list)); tumor_rows=load_tumor_type_rows(Path(a.tumor_type_list))
    sel_models=resolve_requested_ids(csv_text=a.models,file_path=a.models_file,rows=model_rows,key_field="model_id")
    sel_types=resolve_requested_ids(csv_text=a.tumor_types,file_path=a.tumor_types_file,rows=tumor_rows,key_field="tumor_type_id")
    tumor_by_id=row_map(tumor_rows,"tumor_type_id")
    out_root=Path(a.out_root).resolve(); outputs_root=out_root/"genesets"
    src_root=repo_root()/"geneset-extractor-dev"/"NCI_GDC_TCGA_ATAC"/"src"
    peak_matrix=require_file(a.peak_matrix_tsv,"peak matrix"); peak_bed=require_file(a.peak_bed,"peak bed")
    sample_meta=require_file(a.sample_metadata_tsv,"sample metadata"); gtf=require_file(a.gtf,"GTF")
    dig_dir=relative_or_absolute_path(a.dig_dir).resolve()
    if not dig_dir.is_dir(): raise SystemExit(f"Missing dig dir: {dig_dir}")
    model_manifest=require_file(a.model_manifest,"model manifest")
    models=[m for m in sel_models if model_group_for(m)=="atac_accessibility"]
    conflicts=[str(outputs_root/tt/"models"/m) for tt in sel_types for m in models if dir_nonempty(outputs_root/tt/"models"/m)]
    if conflicts and not a.overwrite: raise SystemExit("Output exists (use --overwrite):\n"+"\n".join(conflicts))
    for tt in sel_types:
        row=tumor_by_id[tt]; project_id=str(row.get("project_id","")).strip(); label=str(row.get("tumor_type_label","")).strip()
        if not project_id or not label: raise SystemExit(f"Missing project_id/label for {tt}")
        models_root=outputs_root/tt/"models"
        for m in models:
            if a.overwrite and (models_root/m).exists(): shutil.rmtree(models_root/m)
            run_command([str(Path(a.python_bin).resolve()), str(src_root/"run_atac_accessibility_model.py"),
                "--model_id",m,"--tumor_type_id",tt,"--tumor_type_label",label,"--project_id",project_id,
                "--peak_matrix_tsv",str(peak_matrix),"--peak_bed",str(peak_bed),"--sample_metadata_tsv",str(sample_meta),
                "--gtf",str(gtf),"--run_root",str(models_root),"--python_bin",str(Path(a.python_bin).resolve()),
                "--dig_dir",str(dig_dir),"--model_manifest",str(model_manifest)])
    return 0

if __name__=="__main__": raise SystemExit(main())
