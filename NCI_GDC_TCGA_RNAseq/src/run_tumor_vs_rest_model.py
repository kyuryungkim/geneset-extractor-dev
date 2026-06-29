#!/usr/bin/env python3
"""Run one TCGA tumor-vs-rest model and emit extractor outputs.

Mirrors GTEx/src/run_age_binned_model.py. The only structural difference is the
workflow step: instead of the GTEx-specific gtex_age_binned workflow, we call the
generic DIG workflow `rna_de_prepare --comparison_mode group_vs_rest`, which labels
the focal project as the target group and all other samples as "rest" internally.
The extractor step is the same `convert rna_deg_multi` used by GTEx.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

from tcga_rnaseq_selection_io import default_model_manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run one TCGA tumor-vs-rest model.")
    parser.add_argument("--model_id", required=True)
    parser.add_argument("--tumor_type_id", required=True)
    parser.add_argument("--tumor_type_label", required=True)
    parser.add_argument("--project_id", required=True, help="Focal TCGA project, e.g. TCGA-BRCA.")
    parser.add_argument("--counts_tsv")
    parser.add_argument("--sample_metadata_tsv")
    parser.add_argument("--group_column", default="project_id")
    parser.add_argument("--run_root", required=True)
    parser.add_argument("--python_bin", default=sys.executable or "python3")
    parser.add_argument("--organism", default="human", choices=["human", "mouse"])
    parser.add_argument("--genome_build", default="hg38")
    parser.add_argument("--gtf")
    parser.add_argument("--provenance_mirror_local_prefix")
    parser.add_argument("--provenance_mirror_remote_prefix")
    parser.add_argument("--dig_dir", required=True)
    parser.add_argument("--model_manifest", default=str(default_model_manifest_path()))
    parser.add_argument("--write_commands_only", action="store_true")
    parser.add_argument("--write_model_only", action="store_true")
    return parser.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def resolve_input_path(path_value: str | None, *, base_dir: Path) -> str | None:
    if not path_value:
        return None
    path = Path(path_value)
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    return str(path)


def project_token(project_id: str) -> str:
    # TCGA-BRCA -> BRCA
    token = str(project_id).strip()
    return token.split("-", 1)[1] if "-" in token else token


def signature_name(project_id: str) -> str:
    return f"TCGA_{project_token(project_id)}"


def load_model_settings(manifest_path: Path) -> dict[str, dict[str, str]]:
    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    settings = {str(r.get("model_id", "")).strip(): {str(k): str(v) for k, v in r.items()} for r in rows if str(r.get("model_id", "")).strip()}
    if not settings:
        raise SystemExit(f"No model settings found in {manifest_path}")
    return settings


def shell_join(cmd: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in cmd)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def log_line(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(text.rstrip("\n") + "\n")


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_json(path: Path, payload: dict[str, object]) -> None:
    write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def add_constant_comparison_label(deg_tsv: Path, label_column: str, label_value: str) -> None:
    """Add a constant display-label column to deg_long.tsv (cosmetic; for clean GMT
    set names). Mirrors GTEx, whose workflow emits a gmt_comparison_label column."""
    if not deg_tsv.exists():
        return
    rows = read_tsv_rows(deg_tsv)
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    if label_column in fieldnames:
        return
    fieldnames.append(label_column)
    for row in rows:
        row[label_column] = label_value
    with deg_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def model_group(model_id: str) -> str:
    return "".join(ch for ch in str(model_id) if ch.isalpha()) or str(model_id)


def build_model_sidecar_payload(
    *, model_id: str, tumor_type_id: str, tumor_type_label: str, project_id: str, settings: dict[str, str], comparison_label: str | None = None
) -> dict[str, object]:
    return {
        "schema_version": "1",
        "library": "NCI_GDC_TCGA_RNAseq",
        "model_id": model_id,
        "model_group": model_group(model_id),
        "model_label": "tumor_vs_rest",
        "workflow_name": "rna_de_prepare",
        "extractor_name": "rna_deg_multi",
        "parameters": {
            "comparison_mode": "group_vs_rest",
            "de_mode": settings["workflow_de_mode"],
            "balance_groups": settings["workflow_balance_groups"],
            "balance_seed": settings["workflow_balance_seed"],
            "gene_filter_scope": settings["workflow_gene_filter_scope"],
            "backend": settings["workflow_backend"],
            "covariates": settings["workflow_covariates"],
            "postprocess_mode": settings["extractor_postprocess_mode"],
            "score_mode": settings["extractor_score_mode"],
            "select": settings["extractor_select"],
        },
        "inputs": {
            "tumor_type_id": tumor_type_id,
            "tumor_type_label": tumor_type_label,
            "project_id": project_id,
            "organism": "human",
            "genome_build": "hg38",
        },
        "naming": {
            "signature_name": signature_name(project_id),
            "comparison_label": comparison_label or "tumor_vs_rest",
            "comparison_style": "tumor_vs_rest",
            "gene_set_pattern": "TCGA_<tumor_type>_tumor_vs_rest_up|dn",
        },
    }


def write_grouped_model_sidecars(*, extractor_out: Path, **payload_kwargs) -> None:
    manifest_path = extractor_out / "manifest.tsv"
    if not manifest_path.exists():
        return
    for row in read_tsv_rows(manifest_path):
        meta_rel = str(row.get("meta_path", "")).strip()
        if not meta_rel:
            continue
        sidecar_path = (extractor_out / meta_rel).with_name("geneset.model.json")
        write_json(sidecar_path, build_model_sidecar_payload(comparison_label=str(row.get("label", "")).strip(), **payload_kwargs))


def build_workflow_cmd(*, python_bin, workflow_out, organism, genome_build, project_id, counts_tsv, sample_metadata_tsv, group_column, settings, provenance_mirror_local_prefix, provenance_mirror_remote_prefix) -> list[str]:
    cmd = [
        python_bin, "-m", "geneset_extractors.cli", "workflows", "rna_de_prepare",
        "--modality", "bulk",
        "--counts_tsv", str(counts_tsv),
        "--matrix_orientation", "gene_by_sample",
        "--feature_id_column", "gene_id",
        "--matrix_gene_symbol_column", "gene_symbol",
        "--sample_id_column", "sample_id",
        "--sample_metadata_tsv", str(sample_metadata_tsv),
        "--group_column", group_column,
        "--comparison_mode", "group_vs_rest",
        "--condition_a", project_id,
        "--de_mode", settings["workflow_de_mode"],
        "--balance_groups", settings["workflow_balance_groups"],
        "--balance_seed", settings["workflow_balance_seed"],
        "--gene_filter_scope", settings["workflow_gene_filter_scope"],
        "--backend", settings["workflow_backend"],
        "--out_dir", str(workflow_out),
        "--organism", organism,
        "--genome_build", genome_build,
    ]
    if settings["workflow_covariates"] not in {"none", "NA", ""}:
        cmd += ["--covariates", settings["workflow_covariates"]]
    if provenance_mirror_local_prefix:
        cmd += ["--provenance_mirror_local_prefix", provenance_mirror_local_prefix]
    if provenance_mirror_remote_prefix:
        cmd += ["--provenance_mirror_remote_prefix", provenance_mirror_remote_prefix]
    return cmd


def build_extractor_cmd(*, python_bin, deg_tsv, extractor_out, organism, genome_build, project_id, settings, gtf_path, provenance_mirror_local_prefix, provenance_mirror_remote_prefix) -> list[str]:
    cmd = [
        python_bin, "-m", "geneset_extractors.cli", "convert", "rna_deg_multi",
        "--deg_tsv", str(deg_tsv),
        "--comparison_column", "comparison_id",
        "--comparison_name_column", "gmt_comparison_label",
        "--out_dir", str(extractor_out),
        "--organism", organism,
        "--genome_build", genome_build,
        "--signature_name", signature_name(project_id),
        "--postprocess_mode", settings["extractor_postprocess_mode"],
        "--score_mode", settings["extractor_score_mode"],
        "--select", settings["extractor_select"],
        "--normalize", "within_set_l1",
        "--emit_full", "true",
        "--emit_gmt", "true",
        "--gmt_split_signed", "true",
        "--gmt_name_separator", "_",
        "--gmt_signed_labels", "up_dn",
        "--gmt_require_symbol", settings["extractor_gmt_require_symbol"],
        "--emit_small_gene_sets", settings["extractor_emit_small_gene_sets"],
    ]
    if settings["extractor_disable_default_excludes"] == "true":
        cmd.append("--disable_default_excludes")
    for flag_name, key in [
        ("--padj_max", "extractor_padj_max"),
        ("--pvalue_max", "extractor_pvalue_max"),
        ("--min_abs_logfc", "extractor_min_abs_logfc"),
        ("--top_k", "extractor_top_k"),
        ("--min_score", "extractor_min_score"),
        ("--gmt_source", "extractor_gmt_source"),
        ("--gmt_topk_list", "extractor_gmt_topk_list"),
        ("--gmt_min_genes", "extractor_gmt_min_genes"),
        ("--gmt_max_genes", "extractor_gmt_max_genes"),
    ]:
        value = settings[key]
        if value != "NA":
            cmd += [flag_name, value]
    if settings["extractor_gmt_biotype_allowlist"] not in {"NA", ""}:
        cmd += ["--gmt_biotype_allowlist", settings["extractor_gmt_biotype_allowlist"]]
    if gtf_path:
        cmd += ["--gtf", gtf_path]
    if provenance_mirror_local_prefix:
        cmd += ["--provenance_mirror_local_prefix", provenance_mirror_local_prefix]
    if provenance_mirror_remote_prefix:
        cmd += ["--provenance_mirror_remote_prefix", provenance_mirror_remote_prefix]
    return cmd


def build_provenance_rebuild_cmd(*, python_bin, metadata_json, upstream_provenance_graph_json, provenance_out, provenance_mirror_local_prefix, provenance_mirror_remote_prefix) -> list[str]:
    cmd = [
        python_bin, "-m", "geneset_extractors.cli", "provenance", "build",
        str(metadata_json), "--out", str(provenance_out),
        "--upstream_provenance_graph_json", str(upstream_provenance_graph_json),
    ]
    if provenance_mirror_local_prefix:
        cmd += ["--provenance_mirror_local_prefix", provenance_mirror_local_prefix]
    if provenance_mirror_remote_prefix:
        cmd += ["--provenance_mirror_remote_prefix", provenance_mirror_remote_prefix]
    return cmd


def rebuild_grouped_provenance(*, python_bin, extractor_out, upstream_provenance_graph_json, provenance_mirror_local_prefix, provenance_mirror_remote_prefix, cwd, env, log_path) -> None:
    manifest_path = extractor_out / "manifest.tsv"
    if not manifest_path.exists():
        return
    for row in read_tsv_rows(manifest_path):
        meta_rel = str(row.get("meta_path", "")).strip()
        prov_rel = str(row.get("provenance_path", "")).strip()
        if not meta_rel or not prov_rel:
            continue
        run_command(
            build_provenance_rebuild_cmd(
                python_bin=python_bin,
                metadata_json=extractor_out / meta_rel,
                upstream_provenance_graph_json=upstream_provenance_graph_json,
                provenance_out=extractor_out / prov_rel,
                provenance_mirror_local_prefix=provenance_mirror_local_prefix,
                provenance_mirror_remote_prefix=provenance_mirror_remote_prefix,
            ),
            cwd=cwd, env=env, log_path=log_path,
        )


def run_command(cmd: list[str], *, cwd: Path, env: dict[str, str], log_path: Path) -> None:
    log_line(log_path, f"$ {shell_join(cmd)}")
    completed = subprocess.run(cmd, cwd=str(cwd), env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    if completed.stdout:
        log_line(log_path, completed.stdout.rstrip("\n"))
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(completed.returncode, cmd)


def write_model_commands(*, model_out, model_id, workflow_cmd, extractor_cmd, dig_dir) -> None:
    text = "\n".join([
        f"# Commands For {model_id}", "",
        "## Workflow", "", "```bash",
        f"cd {shlex.quote(str(dig_dir))}",
        f"PYTHONPATH={shlex.quote(str(dig_dir / 'src'))} {shell_join(workflow_cmd)}",
        "```", "",
        "## Extractor", "", "```bash",
        f"cd {shlex.quote(str(dig_dir))}",
        f"PYTHONPATH={shlex.quote(str(dig_dir / 'src'))} {shell_join(extractor_cmd)}",
        "```", "",
        "Note: `rna_deg_multi` writes grouped extractor outputs; this wrapper rebuilds per-group provenance from `manifest.tsv`.",
    ])
    write_text(model_out / "commands.md", text)


def main() -> int:
    args = parse_args()
    repo = repo_root()
    run_root = Path(args.run_root).resolve()
    dig_dir = Path(args.dig_dir).resolve()
    manifest_path = Path(args.model_manifest).resolve()
    resolved_gtf = resolve_input_path(args.gtf, base_dir=repo)
    model_settings = load_model_settings(manifest_path)
    if args.model_id not in model_settings:
        raise SystemExit(f"Unsupported model_id: {args.model_id}")
    settings = model_settings[args.model_id]

    payload_kwargs = dict(
        model_id=args.model_id,
        tumor_type_id=args.tumor_type_id.strip(),
        tumor_type_label=args.tumor_type_label.strip(),
        project_id=args.project_id.strip(),
        settings=settings,
    )

    model_out = run_root / args.model_id
    workflow_out = model_out / "workflow"
    extractor_out = model_out / "extractor"
    model_out.mkdir(parents=True, exist_ok=True)
    model_log = model_out / "run.log"

    if args.write_model_only:
        write_json(extractor_out / "geneset.model.json", build_model_sidecar_payload(**payload_kwargs))
        return 0

    if not args.counts_tsv or not args.sample_metadata_tsv:
        raise SystemExit("--counts_tsv and --sample_metadata_tsv are required unless --write_model_only is used")

    counts_tsv = Path(args.counts_tsv).resolve()
    sample_metadata_tsv = Path(args.sample_metadata_tsv).resolve()

    workflow_cmd = build_workflow_cmd(
        python_bin=args.python_bin, workflow_out=workflow_out, organism=args.organism, genome_build=args.genome_build,
        project_id=args.project_id.strip(), counts_tsv=counts_tsv, sample_metadata_tsv=sample_metadata_tsv,
        group_column=args.group_column, settings=settings,
        provenance_mirror_local_prefix=args.provenance_mirror_local_prefix,
        provenance_mirror_remote_prefix=args.provenance_mirror_remote_prefix,
    )
    deg_tsv_for_extractor = workflow_out / "deg_long.tsv"
    extractor_cmd = build_extractor_cmd(
        python_bin=args.python_bin, deg_tsv=deg_tsv_for_extractor, extractor_out=extractor_out,
        organism=args.organism, genome_build=args.genome_build, project_id=args.project_id.strip(),
        settings=settings, gtf_path=resolved_gtf,
        provenance_mirror_local_prefix=args.provenance_mirror_local_prefix,
        provenance_mirror_remote_prefix=args.provenance_mirror_remote_prefix,
    )
    write_model_commands(model_out=model_out, model_id=args.model_id, workflow_cmd=workflow_cmd, extractor_cmd=extractor_cmd, dig_dir=dig_dir)
    if args.write_commands_only:
        return 0

    env = os.environ.copy()
    env["PYTHONPATH"] = str(dig_dir / "src")
    log_line(model_log, f"[run_tumor_vs_rest_model] model_id={args.model_id} project_id={args.project_id}")

    run_command(workflow_cmd, cwd=dig_dir, env=env, log_path=model_log)
    # Add a clean constant display label so GMT sets read TCGA_<TUMOR>_tumor_vs_rest_up|dn.
    add_constant_comparison_label(deg_tsv_for_extractor, "gmt_comparison_label", "tumor_vs_rest")
    run_command(extractor_cmd, cwd=dig_dir, env=env, log_path=model_log)
    write_grouped_model_sidecars(extractor_out=extractor_out, **payload_kwargs)
    rebuild_grouped_provenance(
        python_bin=args.python_bin, extractor_out=extractor_out,
        upstream_provenance_graph_json=workflow_out / "deg_long.provenance_graph.json",
        provenance_mirror_local_prefix=args.provenance_mirror_local_prefix,
        provenance_mirror_remote_prefix=args.provenance_mirror_remote_prefix,
        cwd=dig_dir, env=env, log_path=model_log,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
