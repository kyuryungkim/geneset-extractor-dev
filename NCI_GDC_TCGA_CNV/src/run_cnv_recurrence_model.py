#!/usr/bin/env python3
"""Run one TCGA CNV recurrence model and emit extractor outputs.

Mirrors NCI_GDC_TCGA_RNAseq runners but CNV is a single-step convert (no upstream
workflow, no provenance rebuild). Restricts the merged cohort segments to the focal
project's samples, then calls `convert cnv_gene_extractor --emit_cohort_sets true`,
which emits per-sample amp/del programs plus cohort recurrent_amp/recurrent_del sets.
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

from cnv_selection_io import default_model_manifest_path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run one TCGA CNV recurrence model.")
    p.add_argument("--model_id", required=True)
    p.add_argument("--tumor_type_id", required=True)
    p.add_argument("--tumor_type_label", required=True)
    p.add_argument("--project_id", required=True)
    p.add_argument("--segments_tsv")
    p.add_argument("--sample_metadata_tsv")
    p.add_argument("--sample_id_column", default="sample_id")
    p.add_argument("--project_column", default="project_id")
    p.add_argument("--run_root", required=True)
    p.add_argument("--python_bin", default=sys.executable or "python3")
    p.add_argument("--organism", default="human", choices=["human", "mouse"])
    p.add_argument("--genome_build", default="hg38")
    p.add_argument("--gtf")
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


def load_model_settings(manifest_path: Path) -> dict[str, dict[str, str]]:
    with manifest_path.open("r", encoding="utf-8", newline="") as h:
        rows = list(csv.DictReader(h, delimiter="\t"))
    s = {str(r.get("model_id", "")).strip(): {str(k): str(v) for k, v in r.items()} for r in rows if str(r.get("model_id", "")).strip()}
    if not s:
        raise SystemExit(f"No model settings found in {manifest_path}")
    return s


def shell_join(cmd: list[str]) -> str:
    return " ".join(shlex.quote(p) for p in cmd)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8", newline="\n")


def log_line(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="\n") as h:
        h.write(text.rstrip("\n") + "\n")


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write_json(path: Path, payload: dict[str, object]) -> None:
    write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def filter_project_segments(*, segments_tsv: Path, sample_metadata_tsv: Path, project_id: str,
                            sample_id_column: str, project_column: str, out_path: Path) -> int:
    meta = read_tsv_rows(sample_metadata_tsv)
    if not meta or sample_id_column not in meta[0] or project_column not in meta[0]:
        raise SystemExit(f"sample metadata must contain '{sample_id_column}' and '{project_column}'")
    keep_samples = {str(r[sample_id_column]).strip() for r in meta if str(r.get(project_column, "")).strip() == project_id}
    if not keep_samples:
        raise SystemExit(f"{project_id}: no samples found in metadata")
    seg = read_tsv_rows(segments_tsv)
    if not seg:
        raise SystemExit(f"Empty segments file: {segments_tsv}")
    if sample_id_column not in seg[0]:
        raise SystemExit(f"segments TSV must contain a '{sample_id_column}' column")
    kept = [r for r in seg if str(r.get(sample_id_column, "")).strip() in keep_samples]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(h, delimiter="\t", fieldnames=list(seg[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(kept)
    return len({str(r.get(sample_id_column, "")).strip() for r in kept})


def build_model_sidecar_payload(*, model_id, tumor_type_id, tumor_type_label, project_id, settings, label=None) -> dict[str, object]:
    return {
        "schema_version": "1",
        "library": "NCI_GDC_TCGA_CNV",
        "model_id": model_id,
        "model_group": "".join(c for c in str(model_id) if c.isalpha()) or str(model_id),
        "model_label": "cnv_recurrence",
        "workflow_name": "cnv_gene_extractor",
        "extractor_name": "cnv_gene_extractor",
        "parameters": {
            "segments_format": settings["segments_format"],
            "program_preset": settings["program_preset"],
            "min_abs_amplitude": settings["min_abs_amplitude"],
            "select": settings["select"],
            "top_k": settings["top_k"],
            "emit_cohort_sets": settings["emit_cohort_sets"],
            "cohort_score_threshold": settings["cohort_score_threshold"],
            "cohort_min_fraction": settings["cohort_min_fraction"],
            "cohort_min_samples": settings["cohort_min_samples"],
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
            "comparison_label": label or "cnv_recurrence",
            "comparison_style": "cnv_recurrence",
            "gene_set_pattern": "cnv_gene_extractor__sample=<sample>__program=amp|del ; cohort: recurrent_amp|recurrent_del",
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
        sidecar = (extractor_out / meta_rel).with_name("geneset.model.json")
        label = str(row.get("label", row.get("program", ""))).strip()
        write_json(sidecar, build_model_sidecar_payload(label=label, **payload_kwargs))


def build_extractor_cmd(*, python_bin, segments_tsv, extractor_out, organism, genome_build, gtf, sample_id_column, settings) -> list[str]:
    cmd = [
        python_bin, "-m", "geneset_extractors.cli", "convert", "cnv_gene_extractor",
        "--segments_tsv", str(segments_tsv),
        "--gtf", str(gtf),
        "--out_dir", str(extractor_out),
        "--organism", organism,
        "--genome_build", genome_build,
        "--segments_format", settings["segments_format"],
        "--sample_id_column", sample_id_column,
        "--program_preset", settings["program_preset"],
        "--min_abs_amplitude", settings["min_abs_amplitude"],
        "--select", settings["select"],
        "--top_k", settings["top_k"],
        "--emit_cohort_sets", settings["emit_cohort_sets"],
        "--cohort_score_threshold", settings["cohort_score_threshold"],
        "--cohort_min_fraction", settings["cohort_min_fraction"],
        "--cohort_min_samples", settings["cohort_min_samples"],
        "--gmt_min_genes", settings["gmt_min_genes"],
        "--gmt_max_genes", settings["gmt_max_genes"],
    ]
    if settings.get("gmt_biotype_allowlist", "NA") not in {"NA", ""}:
        cmd += ["--gmt_biotype_allowlist", settings["gmt_biotype_allowlist"]]
    return cmd


def run_command(cmd: list[str], *, cwd: Path, env: dict[str, str], log_path: Path) -> None:
    log_line(log_path, f"$ {shell_join(cmd)}")
    done = subprocess.run(cmd, cwd=str(cwd), env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    if done.stdout:
        log_line(log_path, done.stdout.rstrip("\n"))
    if done.returncode != 0:
        raise subprocess.CalledProcessError(done.returncode, cmd)


def write_model_commands(*, model_out, model_id, extractor_cmd, dig_dir) -> None:
    text = "\n".join([
        f"# Commands For {model_id}", "", "## Extractor (single-step; CNV has no upstream prepare)", "",
        "```bash", f"cd {shlex.quote(str(dig_dir))}",
        f"PYTHONPATH={shlex.quote(str(dig_dir / 'src'))} {shell_join(extractor_cmd)}", "```", "",
        "Note: cnv_gene_extractor emits per-sample amp/del programs + cohort recurrent_amp/recurrent_del; "
        "provenance is embedded in geneset.meta.json (no separate provenance build).",
    ])
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

    if not args.segments_tsv or not args.sample_metadata_tsv or not args.gtf:
        raise SystemExit("--segments_tsv, --sample_metadata_tsv, --gtf are required unless --write_model_only")

    project_segments = workflow_out / "segments.project.tsv"
    n_samples = filter_project_segments(
        segments_tsv=Path(args.segments_tsv).resolve(), sample_metadata_tsv=Path(args.sample_metadata_tsv).resolve(),
        project_id=args.project_id.strip(), sample_id_column=args.sample_id_column, project_column=args.project_column,
        out_path=project_segments,
    )
    log_line(model_log, f"[run_cnv_recurrence_model] model_id={args.model_id} project_id={args.project_id} n_samples={n_samples}")
    if n_samples < 1:
        raise SystemExit(f"{args.project_id}: no CNV samples after filtering")

    extractor_cmd = build_extractor_cmd(
        python_bin=args.python_bin, segments_tsv=project_segments, extractor_out=extractor_out,
        organism=args.organism, genome_build=args.genome_build, gtf=Path(args.gtf).resolve(),
        sample_id_column=args.sample_id_column, settings=s,
    )
    write_model_commands(model_out=model_out, model_id=args.model_id, extractor_cmd=extractor_cmd, dig_dir=dig_dir)
    if args.write_commands_only:
        return 0

    env = os.environ.copy()
    env["PYTHONPATH"] = str(dig_dir / "src")
    run_command(extractor_cmd, cwd=dig_dir, env=env, log_path=model_log)
    write_grouped_model_sidecars(extractor_out=extractor_out, **payload_kwargs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
