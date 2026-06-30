"""Config IO + selection helpers for the NCI_GDC_TCGA_Methylation library.

Mirrors the other NCI_GDC_TCGA_* selection_io modules. Partition = TCGA tumor type.
"""
from __future__ import annotations

import csv
import os
from pathlib import Path


def repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def library_root() -> Path:
    return repo_root() / "geneset-extractor-dev" / "NCI_GDC_TCGA_Methylation"


def config_root() -> Path:
    o = os.environ.get("TCGA_METH_CONFIG_ROOT", "").strip()
    return Path(o).expanduser().resolve() if o else library_root() / "config"


def default_out_root() -> Path:
    o = os.environ.get("TCGA_METH_OUT_ROOT", "").strip()
    return Path(o).expanduser().resolve() if o else Path.cwd() / "tcga_meth_outputs"


def default_model_list_path() -> Path:
    o = os.environ.get("TCGA_METH_MODEL_LIST", "").strip()
    return Path(o).expanduser().resolve() if o else config_root() / "model_list.tsv"


def default_tumor_type_list_path() -> Path:
    o = os.environ.get("TCGA_METH_TUMOR_TYPE_LIST", "").strip()
    return Path(o).expanduser().resolve() if o else config_root() / "tumor_type_list.tsv"


def default_model_manifest_path() -> Path:
    o = os.environ.get("TCGA_METH_MODEL_MANIFEST", "").strip()
    return Path(o).expanduser().resolve() if o else config_root() / "model_manifest.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def load_model_rows(path: Path | None = None) -> list[dict[str, str]]:
    return read_tsv(path or default_model_list_path())


def load_tumor_type_rows(path: Path | None = None) -> list[dict[str, str]]:
    return read_tsv(path or default_tumor_type_list_path())


def ordered_ids(rows, key_field):
    return [str(r.get(key_field, "")).strip() for r in rows if str(r.get(key_field, "")).strip()]


def enabled_ids(rows, key_field):
    out = []
    for r in rows:
        i = str(r.get(key_field, "")).strip()
        if i and str(r.get("enabled", "true")).strip().lower() == "true":
            out.append(i)
    return out


def parse_id_file(path, key_field):
    rows = []
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        tok = line.split("\t", 1)[0].strip()
        if tok and tok != key_field:
            rows.append(tok)
    return rows


def resolve_requested_ids(*, csv_text, file_path, rows, key_field):
    all_ids = ordered_ids(rows, key_field)
    default_ids = enabled_ids(rows, key_field)
    if csv_text and file_path:
        raise SystemExit(f"Use only one of --{key_field}s or --{key_field}s_file")
    if file_path:
        requested = parse_id_file(file_path, key_field)
    elif csv_text and csv_text.strip().lower() != "all":
        requested = [i.strip() for i in csv_text.split(",") if i.strip()]
    else:
        requested = list(default_ids)
    if not requested:
        raise SystemExit(f"No {key_field} values selected")
    known = set(all_ids)
    unknown = [i for i in requested if i not in known]
    if unknown:
        raise SystemExit(f"Unknown {key_field} values: {', '.join(unknown)}")
    order = {i: n for n, i in enumerate(all_ids)}
    return sorted(set(requested), key=lambda i: order[i])


def row_map(rows, key_field):
    return {str(r.get(key_field, "")).strip(): r for r in rows if str(r.get(key_field, "")).strip()}


def model_group_for(model_id: str) -> str:
    if model_id.startswith("MD"):
        return "methylation_diff"
    raise SystemExit(f"Unsupported model prefix for {model_id}")


def relative_or_absolute_path(path_text: str) -> Path:
    p = Path(path_text)
    return p if p.is_absolute() else repo_root() / p
