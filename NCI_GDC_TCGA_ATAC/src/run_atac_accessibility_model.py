#!/usr/bin/env python3
"""Run one TCGA ATAC accessibility model (cancer-type-vs-rest) and emit extractor outputs.

Single-step, self-contained converter (no DIG prepare). The wrapper derives a binary
'cohort' condition column (focal cancer type vs pooled rest) from sample metadata, then
calls `convert atac_bulk_matrix --condition_a <PROJECT> --condition_b rest`, which computes
the log2FC contrast internally and emits OPEN/CLOSE accessibility gene sets.
"""
from __future__ import annotations
import argparse, csv, json, os, shlex, subprocess, sys
from pathlib import Path
from atac_selection_io import default_model_manifest_path

REST = "rest"

def parse_args():
    p=argparse.ArgumentParser(description="Run one TCGA ATAC accessibility model.")
    p.add_argument("--model_id",required=True); p.add_argument("--tumor_type_id",required=True)
    p.add_argument("--tumor_type_label",required=True); p.add_argument("--project_id",required=True)
    p.add_argument("--peak_matrix_tsv"); p.add_argument("--peak_bed"); p.add_argument("--sample_metadata_tsv")
    p.add_argument("--gtf"); p.add_argument("--sample_id_column",default="sample_id")
    p.add_argument("--cancer_type_column",default="cancer_type")
    p.add_argument("--run_root",required=True); p.add_argument("--python_bin",default=sys.executable or "python3")
    p.add_argument("--organism",default="human",choices=["human","mouse"])
    p.add_argument("--dig_dir",required=True); p.add_argument("--model_manifest",default=str(default_model_manifest_path()))
    p.add_argument("--write_commands_only",action="store_true"); p.add_argument("--write_model_only",action="store_true")
    return p.parse_args()

def repo_root()->Path: return Path(__file__).resolve().parents[3]
def project_token(pid): t=str(pid).strip(); return t.split("-",1)[1] if "-" in t else t
def signature_name(pid): return f"TCGA_{project_token(pid)}"
def load_settings(path):
    with Path(path).open("r",encoding="utf-8",newline="") as h: rows=list(csv.DictReader(h,delimiter="\t"))
    s={str(r.get("model_id","")).strip():{str(k):str(v) for k,v in r.items()} for r in rows if str(r.get("model_id","")).strip()}
    if not s: raise SystemExit(f"No model settings in {path}")
    return s
def shell_join(c): return " ".join(shlex.quote(x) for x in c)
def write_text(p,t): p.parent.mkdir(parents=True,exist_ok=True); p.write_text(t,encoding="utf-8",newline="\n")
def log_line(p,t):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open("a",encoding="utf-8",newline="\n") as h: h.write(t.rstrip("\n")+"\n")
def write_json(p,payload): write_text(p, json.dumps(payload,indent=2,sort_keys=True)+"\n")
def read_tsv_rows(p):
    with Path(p).open("r",encoding="utf-8",newline="") as h: return list(csv.DictReader(h,delimiter="\t"))

def build_cohort_metadata(*, sample_metadata_tsv, project_id, sample_id_column, cancer_type_column, out_path):
    rows=read_tsv_rows(sample_metadata_tsv)
    if not rows or sample_id_column not in rows[0] or cancer_type_column not in rows[0]:
        raise SystemExit(f"sample metadata needs '{sample_id_column}' and '{cancer_type_column}'")
    n_focal=0
    fns=list(rows[0].keys())+(["cohort"] if "cohort" not in rows[0] else [])
    out_path.parent.mkdir(parents=True,exist_ok=True)
    with out_path.open("w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,delimiter="\t",fieldnames=fns,lineterminator="\n"); w.writeheader()
        for r in rows:
            focal = str(r.get(cancer_type_column,"")).strip()==project_id
            r=dict(r); r["cohort"]=project_id if focal else REST
            if focal: n_focal+=1
            w.writerow(r)
    return n_focal, len(rows)-n_focal

def build_model_sidecar_payload(*, model_id, tumor_type_id, tumor_type_label, project_id, settings):
    return {"schema_version":"1","library":"NCI_GDC_TCGA_ATAC","model_id":model_id,
        "model_group":"".join(c for c in str(model_id) if c.isalpha()) or str(model_id),
        "model_label":"atac_accessibility","workflow_name":"atac_bulk_matrix","extractor_name":"atac_bulk_matrix",
        "parameters":{"contrast":"cohort_vs_rest","condition_a":project_id,"condition_b":REST,
            "contrast_metric":settings["contrast_metric"],"peak_weight_transform":settings["peak_weight_transform"],
            "link_method":settings["link_method"],"program_preset":settings["program_preset"],
            "select":settings["select"],"top_k":settings["top_k"],"genome_build":settings["genome_build"]},
        "inputs":{"tumor_type_id":tumor_type_id,"tumor_type_label":tumor_type_label,"project_id":project_id,
            "organism":"human","genome_build":settings["genome_build"]},
        "naming":{"signature_name":signature_name(project_id),"comparison_label":"cohort_vs_rest",
            "comparison_style":"atac_cohort_vs_rest","gene_set_pattern":"atac_bulk_matrix__...__direction=OPEN|CLOSE"}}

def build_extractor_cmd(*, python_bin, peak_matrix, peak_bed, cohort_meta, gtf, extractor_out, organism, project_id, settings):
    return [python_bin,"-m","geneset_extractors.cli","convert","atac_bulk_matrix",
        "--peak_matrix_tsv",str(peak_matrix),"--peak_bed",str(peak_bed),"--sample_metadata_tsv",str(cohort_meta),
        "--gtf",str(gtf),"--out_dir",str(extractor_out),"--organism",organism,"--genome_build",settings["genome_build"],
        "--sample_id_column","sample_id","--condition_column","cohort","--condition_a",project_id,"--condition_b",REST,
        "--contrast_metric",settings["contrast_metric"],"--peak_weight_transform",settings["peak_weight_transform"],
        "--normalize",settings["normalize"],"--link_method",settings["link_method"],"--program_preset",settings["program_preset"],
        "--select",settings["select"],"--top_k",settings["top_k"],"--emit_gmt","true",
        "--gmt_min_genes",settings["gmt_min_genes"],"--gmt_max_genes",settings["gmt_max_genes"],
        "--gmt_biotype_allowlist",settings["gmt_biotype_allowlist"]]

def run_command(cmd,*,cwd,env,log_path):
    log_line(log_path,f"$ {shell_join(cmd)}")
    d=subprocess.run(cmd,cwd=str(cwd),env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,check=False)
    if d.stdout: log_line(log_path,d.stdout.rstrip("\n"))
    if d.returncode!=0: raise subprocess.CalledProcessError(d.returncode,cmd)

def write_model_commands(*,model_out,model_id,extractor_cmd,dig_dir):
    write_text(model_out/"commands.md","\n".join([f"# Commands For {model_id}","","## Extractor (atac_bulk_matrix; self-contained, cancer-type-vs-rest)","",
        "```bash",f"cd {shlex.quote(str(dig_dir))}",f"PYTHONPATH={shlex.quote(str(dig_dir/'src'))} {shell_join(extractor_cmd)}","```","",
        "Note: emits OPEN/CLOSE accessibility gene sets; provenance embedded in geneset.meta.json."]))

def main()->int:
    a=parse_args(); run_root=Path(a.run_root).resolve(); dig_dir=Path(a.dig_dir).resolve()
    settings=load_settings(Path(a.model_manifest).resolve())
    if a.model_id not in settings: raise SystemExit(f"Unsupported model_id: {a.model_id}")
    s=settings[a.model_id]
    pk=dict(model_id=a.model_id,tumor_type_id=a.tumor_type_id.strip(),tumor_type_label=a.tumor_type_label.strip(),project_id=a.project_id.strip(),settings=s)
    model_out=run_root/a.model_id; workflow_out=model_out/"workflow"; extractor_out=model_out/"extractor"
    model_out.mkdir(parents=True,exist_ok=True); model_log=model_out/"run.log"
    if a.write_model_only: write_json(extractor_out/"geneset.model.json", build_model_sidecar_payload(**pk)); return 0
    if not a.peak_matrix_tsv or not a.peak_bed or not a.sample_metadata_tsv or not a.gtf:
        raise SystemExit("--peak_matrix_tsv, --peak_bed, --sample_metadata_tsv, --gtf required unless --write_model_only")
    cohort_meta=workflow_out/"sample_metadata.cohort.tsv"
    n_focal,n_rest=build_cohort_metadata(sample_metadata_tsv=Path(a.sample_metadata_tsv).resolve(),project_id=a.project_id.strip(),
        sample_id_column=a.sample_id_column,cancer_type_column=a.cancer_type_column,out_path=cohort_meta)
    log_line(model_log,f"[run_atac_accessibility_model] model_id={a.model_id} project_id={a.project_id} n_focal={n_focal} n_rest={n_rest}")
    if n_focal<1 or n_rest<1: raise SystemExit(f"{a.project_id}: need focal and rest samples (focal={n_focal}, rest={n_rest})")
    extractor_cmd=build_extractor_cmd(python_bin=a.python_bin,peak_matrix=Path(a.peak_matrix_tsv).resolve(),peak_bed=Path(a.peak_bed).resolve(),
        cohort_meta=cohort_meta,gtf=Path(a.gtf).resolve(),extractor_out=extractor_out,organism=a.organism,project_id=a.project_id.strip(),settings=s)
    write_model_commands(model_out=model_out,model_id=a.model_id,extractor_cmd=extractor_cmd,dig_dir=dig_dir)
    if a.write_commands_only: return 0
    env=os.environ.copy(); env["PYTHONPATH"]=str(dig_dir/"src")
    run_command(extractor_cmd,cwd=dig_dir,env=env,log_path=model_log)
    write_json(extractor_out/"geneset.model.json", build_model_sidecar_payload(**pk))
    return 0

if __name__=="__main__": raise SystemExit(main())
