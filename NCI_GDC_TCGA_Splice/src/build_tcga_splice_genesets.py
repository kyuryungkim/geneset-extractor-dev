#!/usr/bin/env python3
"""Master loop for NCI_GDC_TCGA_Splice. For each normal-bearing tumor type x model,
dispatch to run_splice_tumor_vs_normal_model.py (workflows splice_prepare_public ->
convert splice_event_matrix --tool_family tcga_spliceseq). Thin wrapper; DIG owns logic;
no DIG change. Per-project raw TCGA SpliceSeq PSI files live under --psi_dir/<project_id>.psi.tsv."""
from __future__ import annotations
import argparse, subprocess, shutil, sys
from pathlib import Path
from splice_selection_io import (default_model_list_path, default_model_manifest_path, default_out_root,
    default_tumor_type_list_path, load_model_rows, load_tumor_type_rows, model_group_for,
    relative_or_absolute_path, repo_root, resolve_requested_ids, row_map)

def build_parser():
    p=argparse.ArgumentParser(description="Build TCGA SpliceSeq tumor-vs-normal gene sets.")
    p.add_argument("--models",default="all"); p.add_argument("--models_file")
    p.add_argument("--tumor_types",default="all"); p.add_argument("--tumor_types_file")
    p.add_argument("--model_list",default=str(default_model_list_path()))
    p.add_argument("--tumor_type_list",default=str(default_tumor_type_list_path()))
    p.add_argument("--model_manifest",default=str(default_model_manifest_path()))
    p.add_argument("--python_bin",default=sys.executable or "python3")
    p.add_argument("--psi_dir",required=True,help="Dir with per-project TCGA SpliceSeq PSI files named <project_id>.psi.tsv")
    p.add_argument("--dig_dir",required=True); p.add_argument("--out_root",default=str(default_out_root()))
    p.add_argument("--overwrite",action="store_true")
    return p
def run_command(c): print("$ "+" ".join(c),flush=True); subprocess.run(c,check=True)
def dir_nonempty(p): return p.exists() and p.is_dir() and any(p.iterdir())
def has_normal(row): return str(row.get("has_solid_tissue_normal","")).strip().lower() in {"true","1","yes"}
def main()->int:
    a=build_parser().parse_args()
    model_rows=load_model_rows(Path(a.model_list)); tumor_rows=load_tumor_type_rows(Path(a.tumor_type_list))
    sel_models=resolve_requested_ids(csv_text=a.models,file_path=a.models_file,rows=model_rows,key_field="model_id")
    sel_types=resolve_requested_ids(csv_text=a.tumor_types,file_path=a.tumor_types_file,rows=tumor_rows,key_field="tumor_type_id")
    tumor_by_id=row_map(tumor_rows,"tumor_type_id")
    out_root=Path(a.out_root).resolve(); outputs_root=out_root/"genesets"
    src_root=repo_root()/"geneset-extractor-dev"/"NCI_GDC_TCGA_Splice"/"src"
    psi_dir=relative_or_absolute_path(a.psi_dir).resolve()
    if not psi_dir.is_dir(): raise SystemExit(f"Missing psi_dir: {psi_dir}")
    dig_dir=relative_or_absolute_path(a.dig_dir).resolve()
    if not dig_dir.is_dir(): raise SystemExit(f"Missing dig dir: {dig_dir}")
    model_manifest=relative_or_absolute_path(a.model_manifest).resolve()
    models=[m for m in sel_models if model_group_for(m)=="splice_tumor_vs_normal"]
    conflicts=[str(outputs_root/tt/"models"/m) for tt in sel_types if has_normal(tumor_by_id[tt]) for m in models if dir_nonempty(outputs_root/tt/"models"/m)]
    if conflicts and not a.overwrite: raise SystemExit("Output exists (use --overwrite):\n"+"\n".join(conflicts))
    for tt in sel_types:
        row=tumor_by_id[tt]
        if not has_normal(row):
            for m in models: print(f"  skip {tt}/{m}: no matched solid tissue normal",flush=True)
            continue
        project_id=str(row.get("project_id","")).strip(); label=str(row.get("tumor_type_label","")).strip()
        psi_tsv=psi_dir/f"{project_id}.psi.tsv"
        if not psi_tsv.is_file():
            for m in models: print(f"  skip {tt}/{m}: missing PSI file {psi_tsv}",flush=True)
            continue
        models_root=outputs_root/tt/"models"
        for m in models:
            if a.overwrite and (models_root/m).exists(): shutil.rmtree(models_root/m)
            run_command([str(Path(a.python_bin).resolve()), str(src_root/"run_splice_tumor_vs_normal_model.py"),
                "--model_id",m,"--tumor_type_id",tt,"--tumor_type_label",label,"--project_id",project_id,
                "--psi_tsv",str(psi_tsv),"--run_root",str(models_root),"--python_bin",str(Path(a.python_bin).resolve()),
                "--dig_dir",str(dig_dir),"--model_manifest",str(model_manifest)])
    return 0
if __name__=="__main__": raise SystemExit(main())
