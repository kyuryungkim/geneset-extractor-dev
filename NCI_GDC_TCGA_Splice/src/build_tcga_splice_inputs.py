#!/usr/bin/env python3
"""Stage per-project TCGA SpliceSeq PSI files for the splice library.

TCGA SpliceSeq (MD Anderson) distributes per-tumor-type PSI matrices: rows = splice events
(symbol, exons, splice_type, from_exon, to_exon, ...), columns = samples (TCGA barcodes; normal
samples carry a _Norm suffix, which splice_prepare_public uses to infer adjacent_normal). This
helper copies/renames downloaded per-type files into <out_dir>/<project_id>.psi.tsv, the layout
the master loop expects (--psi_dir). Downloaded inputs live OUTSIDE the committed tree.

Provide a mapping TSV (--map_tsv) with columns: project_id, psi_file (path to the raw per-type file).
"""
from __future__ import annotations
import argparse, csv, shutil
from pathlib import Path
def parse_args():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--map_tsv",required=True,help="TSV: project_id, psi_file")
    p.add_argument("--out_dir",required=True)
    return p.parse_args()
def main()->int:
    a=parse_args(); out=Path(a.out_dir).expanduser().resolve(); out.mkdir(parents=True,exist_ok=True); n=0
    with open(a.map_tsv,encoding="utf-8",newline="") as h:
        for r in csv.DictReader(h,delimiter="\t"):
            pid=str(r.get("project_id","")).strip(); src=str(r.get("psi_file","")).strip()
            if not pid or not src: continue
            sp=Path(src).expanduser().resolve()
            if not sp.is_file(): print(f"  WARN: missing {sp}; skipping {pid}",flush=True); continue
            shutil.copyfile(sp, out/f"{pid}.psi.tsv"); n+=1
    print(f"Staged {n} per-project PSI files -> {out}",flush=True)
    return 0
if __name__=="__main__": raise SystemExit(main())
