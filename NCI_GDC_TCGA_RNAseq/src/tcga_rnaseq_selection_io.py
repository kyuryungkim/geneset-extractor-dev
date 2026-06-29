"""Config IO + selection helpers for the NCI_GDC_TCGA_RNAseq library.

Mirrors geneset-extractor-dev/GTEx/src/selection_io.py. The partition dimension
here is the TCGA tumor type (33 projects) rather than GTEx tissue.
"""
from __future__ import annotations

import csv
import os
import re
from pathlib import Path


def repo_root() -> Path:
    # repo/geneset-extractor-dev/NCI_GDC_TCGA_RNAseq/src/<this file>
    return Path(__file__).resolve().parents[3]


def library_root() -> Path:
    return repo_root() / "geneset-extractor-dev" / "NCI_GDC_TCGA_RNAseq"


def config_root() -> Path:
    override = os.environ.get("TCGA_RNASEQ_CONFIG_ROOT", "").strip()
    if override:
        return Path(override).expanduser().resolve()
    return library_root() / "config"


def planning_root() -> Path:
    return library_root() / "planning"


def default_out_root() -> Path:
    override = os.environ.get("TCGA_RNASEQ_OUT_ROOT", "").strip()
    if override:
        return Path(override).expanduser().resolve()
    return Path.cwd() / "tcga_rnaseq_outputs"


def default_genesets_root() -> Path:
    return default_out_root() / "genesets"


def default_model_list_path() -> Path:
    override = os.environ.get("TCGA_RNASEQ_MODEL_LIST", "").strip()
    if override:
        return Path(override).expanduser().resolve()
    return config_root() / "model_list.tsv"


def default_tumor_type_list_path() -> Path:
    override = os.environ.get("TCGA_RNASEQ_TUMOR_TYPE_LIST", "").strip()
    if override:
        return Path(override).expanduser().resolve()
    return config_root() / "tumor_type_list.tsv"


def default_model_manifest_path() -> Path:
    override = os.environ.get("TCGA_RNASEQ_MODEL_MANIFEST", "").strip()
    if override:
        return Path(override).expanduser().resolve()
    return config_root() / "model_manifest.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_model_rows(path: Path | None = None) -> list[dict[str, str]]:
    return read_tsv(path or default_model_list_path())


def load_tumor_type_rows(path: Path | None = None) -> list[dict[str, str]]:
    return read_tsv(path or default_tumor_type_list_path())


def ordered_ids(rows: list[dict[str, str]], key_field: str) -> list[str]:
    return [str(row.get(key_field, "")).strip() for row in rows if str(row.get(key_field, "")).strip()]


def enabled_ids(rows: list[dict[str, str]], key_field: str) -> list[str]:
    enabled: list[str] = []
    for row in rows:
        item_id = str(row.get(key_field, "")).strip()
        if not item_id:
            continue
        if str(row.get("enabled", "true")).strip().lower() == "true":
            enabled.append(item_id)
    return enabled


def parse_id_file(path: Path, key_field: str) -> list[str]:
    text = path.read_text(encoding="utf-8")
    rows: list[str] = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        token = line.split("\t", 1)[0].strip()
        if not token or token == key_field:
            continue
        rows.append(token)
    return rows


def resolve_requested_ids(
    *,
    csv_text: str | None,
    file_path: str | None,
    rows: list[dict[str, str]],
    key_field: str,
) -> list[str]:
    all_ids = ordered_ids(rows, key_field)
    default_ids = enabled_ids(rows, key_field)
    if csv_text and file_path:
        raise SystemExit(f"Use only one of --{key_field}s or --{key_field}s_file")
    if file_path:
        requested = parse_id_file(Path(file_path), key_field)
    elif csv_text and csv_text.strip().lower() != "all":
        requested = [item.strip() for item in csv_text.split(",") if item.strip()]
    else:
        requested = list(default_ids)
    if not requested:
        raise SystemExit(f"No {key_field} values selected")
    known = set(all_ids)
    unknown = [item for item in requested if item not in known]
    if unknown:
        raise SystemExit(f"Unknown {key_field} values: {', '.join(unknown)}")
    order_index = {item_id: index for index, item_id in enumerate(all_ids)}
    return sorted(set(requested), key=lambda item_id: order_index[item_id])


def row_map(rows: list[dict[str, str]], key_field: str) -> dict[str, dict[str, str]]:
    return {str(row.get(key_field, "")).strip(): row for row in rows if str(row.get(key_field, "")).strip()}


def model_group_for(model_id: str) -> str:
    if model_id.startswith("TR"):
        return "tumor_vs_rest"
    if model_id.startswith("TN"):
        return "tumor_vs_normal"
    raise SystemExit(f"Unsupported model prefix for {model_id}")


def relative_or_absolute_path(path_text: str) -> Path:
    path = Path(path_text)
    if path.is_absolute():
        return path
    return repo_root() / path
