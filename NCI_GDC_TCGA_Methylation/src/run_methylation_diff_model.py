#!/usr/bin/env python3
"""Run one TCGA methylation tumor-vs-normal model and emit extractor outputs.

Two-step like RNA-seq, but the workflow step is the NEW DIG workflow
`methylation_diff_prepare` (computes per-CpG delta-beta from the beta matrix), then
`convert methylation_cpg_diff` builds the gene sets, then `provenance build` rebuilds
the geneset provenance with the upstream graph. Restricts the merged beta matrix to the
focal project's samples first.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

from meth_selection_io import default_model_manifest_path

TUMOR = "Primary Tumor"
NORMAL = "Solid Tissue Normal"


def parse_args():
    p = argparse.ArgumentParser(description="Run one TCGA methylation tumor-vs-normal model.")
    p.add_argument("--model_id", required=True)
    p.add_argument("--tumor_type_id", required=True)
    p.add_argument("--tumor_type_label", required=True)
    p.add_argument("--project_id", required=True)
    p.add_argument("--beta_matrix_tsv")
    p.add_argument("--sample_metadata_tsv")
    p.add_argument("--gtf")
    p.add_argument("--probe_manifest_tsv")
    p.add_argument("--sample_id_column", default="sample_id")
    p.add_argument("--project_column", default="project_id")
    p.add_argument("--group_column", default="sample_type")
    p.add_argument("--run_root", required=True)
    p.add_argument("--python_bin", default=sys.executable or "python3")
    p.add_argument("--organism", default="human", choices=["human", "mouse"])
    p.add_argument("--genome_build", default="hg38")
    p.add_argument("--dig_dir", required=True)
    p.add_argument("--model_manifest", default=str(default_model_manifest_path()))
    p.add_argument("--write_commands_only", action="store_true")
    p.add_argument("--write_model_only", action="store_true")
    return p.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def project_token(project_id: str) -> str:
    t = str(project_id).strip()
    return t.split("-", 1)[1] if "-" in t else t


def signature_name(project_id: str) -> str:
    return f"TCGA_{project_token(project_id)}"


def load_model_settings(path: Path):
    with path.open("r", encoding="utf-8", newline="") as h:
        rows = list(csv.DictReader(h, delimiter="\t"))
    s = {str(r.get("model_id", "")).strip(): {str(k): str(v) for k, v in r.items()} for r in rows if str(r.get("model_id", "")).strip()}
    if not s:
        raise SystemExit(f"No model settings in {path}")
    return s


def shell_join(cmd):
    return " ".join(shlex.quote(c) for c in cmd)


def write_text(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def log_line(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="\n") as h:
        h.write(text.rstrip("\n") + "\n")


def write_json(path: Path, payload):
    write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def read_tsv_rows(path: Path):
    with path.open("r", encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def filter_project_beta(*, beta_matrix_tsv: Path, sample_metadata_tsv: Path, project_id: str,
                        sample_id_column: str, project_column: str, group_column: str, out_path: Path):
    meta = read_tsv_rows(sample_metadata_tsv)
    if not meta or sample_id_column not in meta[0] or project_column not in meta[0] or group_column not in meta[0]:
        raise SystemExit(f"sample metadata needs '{sample_id_column}','{project_column}','{group_column}'")
    grp = {str(r[sample_id_column]).strip(): str(r.get(group_column, "")).strip() for r in meta if str(r.get(project_column, "")).strip() == project_id}
    keep = {s for s, g in grp.items() if g in {TUMOR, NORMAL}}
    n_t = sum(1 for g in grp.values() if g == TUMOR)
    n_n = sum(1 for g in grp.values() if g == NORMAL)
    with beta_matrix_tsv.open("r", encoding="utf-8", newline="") as h:
        reader = csv.reader(h, delimiter="\t")
        header = next(reader)
        keep_idx = [0] + [i for i, c in enumerate(header) if i > 0 and c in keep]
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", encoding="utf-8", newline="") as o:
            w = csv.writer(o, delimiter="\t", lineterminator="\n")
            w.writerow([header[i] for i in keep_idx])
            for row in reader:
                if row:
                    w.writerow([row[i] for i in keep_idx])
    return n_t, n_n


def build_model_sidecar_payload(*, model_id, tumor_type_id, tumor_type_label, project_id, settings):
    return {
        "schema_version": "1", "library": "NCI_GDC_TCGA_Methylation", "model_id": model_id,
        "model_group": "".join(c for c in str(model_id) if c.isalpha()) or str(model_id),
        "model_label": "methylation_diff", "workflow_name": "methylation_diff_prepare",
        "extractor_name": "methylation_cpg_diff",
        "parameters": {
            "comparison_mode": "condition_a_vs_b", "condition_a": TUMOR, "condition_b": NORMAL,
            "array_type": settings["array_type"], "score_mode": settings["score_mode"],
            "delta_orientation": settings["delta_orientation"], "program_preset": settings["program_preset"],
            "link_method": settings["link_method"], "select": settings["select"], "top_k": settings["top_k"],
        },
        "inputs": {"tumor_type_id": tumor_type_id, "tumor_type_label": tumor_type_label,
                   "project_id": project_id, "organism": "human", "genome_build": "hg38"},
        "naming": {"signature_name": signature_name(project_id), "comparison_label": "tumor_vs_normal",
                   "comparison_style": "methylation_tumor_vs_normal",
                   "gene_set_pattern": "<signature>__pos (hypomethylation) | __neg (hypermethylation)"},
    }


def build_workflow_cmd(*, python_bin, beta_project, sample_metadata_tsv, workflow_out, organism, genome_build):
    return [python_bin, "-m", "geneset_extractors.cli", "workflows", "methylation_diff_prepare",
            "--beta_matrix_tsv", str(beta_project), "--sample_metadata_tsv", str(sample_metadata_tsv),
            "--sample_id_column", "sample_id", "--group_column", "sample_type",
            "--comparison_mode", "condition_a_vs_b", "--condition_a", TUMOR, "--condition_b", NORMAL,
            "--out_dir", str(workflow_out), "--organism", organism, "--genome_build", genome_build]


def build_extractor_cmd(*, python_bin, cpg_tsv, extractor_out, organism, genome_build, gtf, probe_manifest_tsv, settings):
    cmd = [python_bin, "-m", "geneset_extractors.cli", "convert", "methylation_cpg_diff",
           "--cpg_tsv", str(cpg_tsv), "--gtf", str(gtf), "--out_dir", str(extractor_out),
           "--organism", organism, "--genome_build", genome_build,
           "--array_type", settings["array_type"], "--probe_id_column", "probe_id",
           "--delta_column", "delta_beta",
           "--score_mode", settings["score_mode"], "--delta_orientation", settings["delta_orientation"],
           "--program_preset", settings["program_preset"], "--link_method", settings["link_method"],
           "--select", settings["select"], "--top_k", settings["top_k"],
           "--emit_gmt", "true", "--gmt_split_signed", "true",
           "--gmt_min_genes", settings["gmt_min_genes"], "--gmt_max_genes", settings["gmt_max_genes"]]
    if settings.get("gmt_biotype_allowlist", "NA") not in {"NA", ""}:
        cmd += ["--gmt_biotype_allowlist", settings["gmt_biotype_allowlist"]]
    if probe_manifest_tsv:
        cmd += ["--probe_manifest_tsv", str(probe_manifest_tsv)]
    return cmd


def build_provenance_cmd(*, python_bin, metadata_json, upstream_graph, provenance_out):
    return [python_bin, "-m", "geneset_extractors.cli", "provenance", "build", str(metadata_json),
            "--out", str(provenance_out), "--upstream_provenance_graph_json", str(upstream_graph)]


def run_command(cmd, *, cwd, env, log_path):
    log_line(log_path, f"$ {shell_join(cmd)}")
    done = subprocess.run(cmd, cwd=str(cwd), env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    if done.stdout:
        log_line(log_path, done.stdout.rstrip("\n"))
    if done.returncode != 0:
        raise subprocess.CalledProcessError(done.returncode, cmd)


def write_model_commands(*, model_out, model_id, workflow_cmd, extractor_cmd, dig_dir):
    text = "\n".join([f"# Commands For {model_id}", "", "## Workflow (methylation_diff_prepare)", "",
                      "```bash", f"cd {shlex.quote(str(dig_dir))}",
                      f"PYTHONPATH={shlex.quote(str(dig_dir / 'src'))} {shell_join(workflow_cmd)}", "```", "",
                      "## Extractor (methylation_cpg_diff)", "", "```bash", f"cd {shlex.quote(str(dig_dir))}",
                      f"PYTHONPATH={shlex.quote(str(dig_dir / 'src'))} {shell_join(extractor_cmd)}", "```", "",
                      "Then `provenance build` rebuilds geneset.provenance.json with the upstream cpg_diff graph."])
    write_text(model_out / "commands.md", text)


def main() -> int:
    args = parse_args()
    run_root = Path(args.run_root).resolve()
    dig_dir = Path(args.dig_dir).resolve()
    settings = load_model_settings(Path(args.model_manifest).resolve())
    if args.model_id not in settings:
        raise SystemExit(f"Unsupported model_id: {args.model_id}")
    s = settings[args.model_id]
    payload_kwargs = dict(model_id=args.model_id, tumor_type_id=args.tumor_type_id.strip(),
                          tumor_type_label=args.tumor_type_label.strip(), project_id=args.project_id.strip(), settings=s)

    model_out = run_root / args.model_id
    workflow_out = model_out / "workflow"
    extractor_out = model_out / "extractor"
    model_out.mkdir(parents=True, exist_ok=True)
    model_log = model_out / "run.log"

    if args.write_model_only:
        write_json(extractor_out / "geneset.model.json", build_model_sidecar_payload(**payload_kwargs))
        return 0
    if not args.beta_matrix_tsv or not args.sample_metadata_tsv or not args.gtf:
        raise SystemExit("--beta_matrix_tsv, --sample_metadata_tsv, --gtf required unless --write_model_only")

    beta_project = workflow_out / "beta.project.tsv"
    n_t, n_n = filter_project_beta(beta_matrix_tsv=Path(args.beta_matrix_tsv).resolve(),
                                   sample_metadata_tsv=Path(args.sample_metadata_tsv).resolve(),
                                   project_id=args.project_id.strip(), sample_id_column=args.sample_id_column,
                                   project_column=args.project_column, group_column=args.group_column, out_path=beta_project)
    log_line(model_log, f"[run_methylation_diff_model] model_id={args.model_id} project_id={args.project_id} n_tumor={n_t} n_normal={n_n}")
    if n_t < 1 or n_n < 1:
        raise SystemExit(f"{args.project_id}: insufficient samples (tumor={n_t}, normal={n_n})")

    workflow_cmd = build_workflow_cmd(python_bin=args.python_bin, beta_project=beta_project,
                                      sample_metadata_tsv=Path(args.sample_metadata_tsv).resolve(),
                                      workflow_out=workflow_out, organism=args.organism, genome_build=args.genome_build)
    cpg_tsv = workflow_out / "cpg_diff.tsv"
    extractor_cmd = build_extractor_cmd(python_bin=args.python_bin, cpg_tsv=cpg_tsv, extractor_out=extractor_out,
                                        organism=args.organism, genome_build=args.genome_build, gtf=Path(args.gtf).resolve(),
                                        probe_manifest_tsv=(Path(args.probe_manifest_tsv).resolve() if args.probe_manifest_tsv else None),
                                        settings=s)
    write_model_commands(model_out=model_out, model_id=args.model_id, workflow_cmd=workflow_cmd, extractor_cmd=extractor_cmd, dig_dir=dig_dir)
    if args.write_commands_only:
        return 0

    env = os.environ.copy()
    env["PYTHONPATH"] = str(dig_dir / "src")
    run_command(workflow_cmd, cwd=dig_dir, env=env, log_path=model_log)
    run_command(extractor_cmd, cwd=dig_dir, env=env, log_path=model_log)
    meta_json = extractor_out / "geneset.meta.json"
    upstream = workflow_out / "cpg_diff.provenance_graph.json"
    if meta_json.exists() and upstream.exists():
        run_command(build_provenance_cmd(python_bin=args.python_bin, metadata_json=meta_json,
                                         upstream_graph=upstream, provenance_out=extractor_out / "geneset.provenance.json"),
                    cwd=dig_dir, env=env, log_path=model_log)
    write_json(extractor_out / "geneset.model.json", build_model_sidecar_payload(**payload_kwargs))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
