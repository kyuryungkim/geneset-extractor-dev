#!/usr/bin/env python3
"""Run one TCGA SpliceSeq tumor-vs-normal model and emit extractor outputs.

Two-step, both steps already in DIG (no DIG change):
  workflows splice_prepare_public  (raw TCGA SpliceSeq PSI -> psi_matrix + sample_metadata,
                                    inferring tumor/adjacent_normal from _Norm suffix)
  convert splice_event_matrix --tool_family tcga_spliceseq --study_contrast condition_a_vs_b
                                    (computes the tumor-vs-normal contrast internally, Welch test)
"""
from __future__ import annotations
import argparse, csv, json, os, shlex, subprocess, sys
from pathlib import Path
from splice_selection_io import default_model_manifest_path

def parse_args():
    p=argparse.ArgumentParser(description="Run one TCGA SpliceSeq tumor-vs-normal model.")
    p.add_argument("--model_id",required=True); p.add_argument("--tumor_type_id",required=True)
    p.add_argument("--tumor_type_label",required=True); p.add_argument("--project_id",required=True)
    p.add_argument("--psi_tsv"); p.add_argument("--run_root",required=True)
    p.add_argument("--python_bin",default=sys.executable or "python3")
    p.add_argument("--organism",default="human",choices=["human","mouse"]); p.add_argument("--genome_build",default="hg38")
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
def write_json(p,payload): write_text(p,json.dumps(payload,indent=2,sort_keys=True)+"\n")

def build_model_sidecar_payload(*,model_id,tumor_type_id,tumor_type_label,project_id,settings):
    return {"schema_version":"1","library":"NCI_GDC_TCGA_Splice","model_id":model_id,
        "model_group":"".join(c for c in str(model_id) if c.isalpha()) or str(model_id),
        "model_label":"splice_tumor_vs_normal","workflow_name":"splice_prepare_public","extractor_name":"splice_event_matrix",
        "parameters":{"tool_family":settings["tool_family"],"study_contrast":settings["study_contrast"],
            "condition_a":settings["condition_a"],"condition_b":settings["condition_b"],
            "effect_metric":settings["effect_metric"],"select":settings["select"],"top_k":settings["top_k"]},
        "inputs":{"tumor_type_id":tumor_type_id,"tumor_type_label":tumor_type_label,"project_id":project_id,
            "organism":"human","genome_build":"hg38"},
        "naming":{"signature_name":signature_name(project_id),"comparison_label":"tumor_vs_normal",
            "comparison_style":"splice_tumor_vs_normal","gene_set_pattern":"<signature>__pos (increased inclusion) | __neg (decreased)"}}

def build_prepare_cmd(*,python_bin,psi_tsv,prepared_out,organism,genome_build,project_id):
    return [python_bin,"-m","geneset_extractors.cli","workflows","splice_prepare_public",
        "--input_mode","tcga_spliceseq","--psi_tsv",str(psi_tsv),"--out_dir",str(prepared_out),
        "--organism",organism,"--genome_build",genome_build,"--study_id",str(project_id),
        "--study_label",f"TCGA {project_token(project_id)} SpliceSeq"]

def build_extractor_cmd(*,python_bin,psi_matrix,sample_metadata,event_metadata,extractor_out,organism,genome_build,signature,settings):
    cmd=[python_bin,"-m","geneset_extractors.cli","convert","splice_event_matrix",
        "--psi_matrix_tsv",str(psi_matrix),"--sample_metadata_tsv",str(sample_metadata),
        "--out_dir",str(extractor_out),"--organism",organism,"--genome_build",genome_build,
        "--tool_family",settings["tool_family"],"--study_contrast",settings["study_contrast"],
        "--condition_column","condition","--condition_a",settings["condition_a"],"--condition_b",settings["condition_b"],
        "--effect_metric",settings["effect_metric"],"--select",settings["select"],"--top_k",settings["top_k"],
        "--signature_name",signature,"--emit_gmt","true","--gmt_split_signed",settings["gmt_split_signed"],
        "--gmt_min_genes",settings["gmt_min_genes"],"--gmt_max_genes",settings["gmt_max_genes"]]
    if event_metadata is not None and Path(event_metadata).exists():
        cmd+=["--event_metadata_tsv",str(event_metadata)]
    return cmd

def run_command(cmd,*,cwd,env,log_path):
    log_line(log_path,f"$ {shell_join(cmd)}")
    d=subprocess.run(cmd,cwd=str(cwd),env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,check=False)
    if d.stdout: log_line(log_path,d.stdout.rstrip("\n"))
    if d.returncode!=0: raise subprocess.CalledProcessError(d.returncode,cmd)

def write_model_commands(*,model_out,model_id,prepare_cmd,extractor_cmd,dig_dir):
    write_text(model_out/"commands.md","\n".join([f"# Commands For {model_id}","",
        "## Prepare (splice_prepare_public)","","```bash",f"cd {shlex.quote(str(dig_dir))}",
        f"PYTHONPATH={shlex.quote(str(dig_dir/'src'))} {shell_join(prepare_cmd)}","```","",
        "## Extractor (splice_event_matrix, tcga_spliceseq tumor-vs-normal)","","```bash",f"cd {shlex.quote(str(dig_dir))}",
        f"PYTHONPATH={shlex.quote(str(dig_dir/'src'))} {shell_join(extractor_cmd)}","```"]))

def main()->int:
    a=parse_args(); run_root=Path(a.run_root).resolve(); dig_dir=Path(a.dig_dir).resolve()
    settings=load_settings(Path(a.model_manifest).resolve())
    if a.model_id not in settings: raise SystemExit(f"Unsupported model_id: {a.model_id}")
    s=settings[a.model_id]
    pk=dict(model_id=a.model_id,tumor_type_id=a.tumor_type_id.strip(),tumor_type_label=a.tumor_type_label.strip(),project_id=a.project_id.strip(),settings=s)
    model_out=run_root/a.model_id; workflow_out=model_out/"workflow"; extractor_out=model_out/"extractor"
    model_out.mkdir(parents=True,exist_ok=True); model_log=model_out/"run.log"
    if a.write_model_only: write_json(extractor_out/"geneset.model.json",build_model_sidecar_payload(**pk)); return 0
    if not a.psi_tsv: raise SystemExit("--psi_tsv required unless --write_model_only")
    prepare_cmd=build_prepare_cmd(python_bin=a.python_bin,psi_tsv=Path(a.psi_tsv).resolve(),prepared_out=workflow_out,
        organism=a.organism,genome_build=a.genome_build,project_id=a.project_id.strip())
    psi_matrix=workflow_out/"psi_matrix.tsv"; sample_metadata=workflow_out/"sample_metadata.tsv"; event_metadata=workflow_out/"event_metadata.tsv"
    extractor_cmd=build_extractor_cmd(python_bin=a.python_bin,psi_matrix=psi_matrix,sample_metadata=sample_metadata,
        event_metadata=event_metadata,extractor_out=extractor_out,organism=a.organism,genome_build=a.genome_build,
        signature=signature_name(a.project_id.strip()),settings=s)
    write_model_commands(model_out=model_out,model_id=a.model_id,prepare_cmd=prepare_cmd,extractor_cmd=extractor_cmd,dig_dir=dig_dir)
    if a.write_commands_only: return 0
    env=os.environ.copy(); env["PYTHONPATH"]=str(dig_dir/"src")
    run_command(prepare_cmd,cwd=dig_dir,env=env,log_path=model_log)
    run_command(extractor_cmd,cwd=dig_dir,env=env,log_path=model_log)
    write_json(extractor_out/"geneset.model.json",build_model_sidecar_payload(**pk))
    return 0
if __name__=="__main__": raise SystemExit(main())
