"""Config IO + selection helpers for the NCI_GDC_TCGA_ATAC library."""
from __future__ import annotations
import csv, os
from pathlib import Path

def repo_root() -> Path: return Path(__file__).resolve().parents[3]
def library_root() -> Path: return repo_root() / "geneset-extractor-dev" / "NCI_GDC_TCGA_ATAC"
def config_root() -> Path:
    o = os.environ.get("TCGA_ATAC_CONFIG_ROOT", "").strip()
    return Path(o).expanduser().resolve() if o else library_root() / "config"
def default_out_root() -> Path:
    o = os.environ.get("TCGA_ATAC_OUT_ROOT", "").strip()
    return Path(o).expanduser().resolve() if o else Path.cwd() / "tcga_atac_outputs"
def default_model_list_path() -> Path:
    o = os.environ.get("TCGA_ATAC_MODEL_LIST", "").strip()
    return Path(o).expanduser().resolve() if o else config_root() / "model_list.tsv"
def default_tumor_type_list_path() -> Path:
    o = os.environ.get("TCGA_ATAC_TUMOR_TYPE_LIST", "").strip()
    return Path(o).expanduser().resolve() if o else config_root() / "tumor_type_list.tsv"
def default_model_manifest_path() -> Path:
    o = os.environ.get("TCGA_ATAC_MODEL_MANIFEST", "").strip()
    return Path(o).expanduser().resolve() if o else config_root() / "model_manifest.tsv"
def read_tsv(path: Path):
    with Path(path).open("r", encoding="utf-8", newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))
def load_model_rows(path=None): return read_tsv(path or default_model_list_path())
def load_tumor_type_rows(path=None): return read_tsv(path or default_tumor_type_list_path())
def ordered_ids(rows, key): return [str(r.get(key,"")).strip() for r in rows if str(r.get(key,"")).strip()]
def enabled_ids(rows, key):
    out=[]
    for r in rows:
        i=str(r.get(key,"")).strip()
        if i and str(r.get("enabled","true")).strip().lower()=="true": out.append(i)
    return out
def parse_id_file(path, key):
    rows=[]
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line=raw.strip()
        if not line or line.startswith("#"): continue
        tok=line.split("\t",1)[0].strip()
        if tok and tok!=key: rows.append(tok)
    return rows
def resolve_requested_ids(*, csv_text, file_path, rows, key_field):
    all_ids=ordered_ids(rows,key_field); default_ids=enabled_ids(rows,key_field)
    if csv_text and file_path: raise SystemExit(f"Use only one of --{key_field}s or --{key_field}s_file")
    if file_path: requested=parse_id_file(file_path,key_field)
    elif csv_text and csv_text.strip().lower()!="all": requested=[i.strip() for i in csv_text.split(",") if i.strip()]
    else: requested=list(default_ids)
    if not requested: raise SystemExit(f"No {key_field} values selected")
    known=set(all_ids); unknown=[i for i in requested if i not in known]
    if unknown: raise SystemExit(f"Unknown {key_field} values: {', '.join(unknown)}")
    order={i:n for n,i in enumerate(all_ids)}
    return sorted(set(requested), key=lambda i: order[i])
def row_map(rows, key): return {str(r.get(key,"")).strip(): r for r in rows if str(r.get(key,"")).strip()}
def model_group_for(model_id: str) -> str:
    if model_id.startswith("AC"): return "atac_accessibility"
    raise SystemExit(f"Unsupported model prefix for {model_id}")
def relative_or_absolute_path(path_text: str) -> Path:
    p=Path(path_text); return p if p.is_absolute() else repo_root()/p
