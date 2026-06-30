#!/usr/bin/env python3
"""Build the peak matrix + peak BED + sample metadata for TCGA ATAC models from the
open-access Corces 2018 (GDC ATACseq-AWG) pan-cancer matrices.

Inputs: the pan-cancer normalized peak-by-sample count matrix (peak id/coords col + sample
cols), the pan-cancer peak BED (562,709 peaks, hg19), and the GDC sample lookup table
(barcode -> cancer type). Emits:
  - peak_matrix_tsv   peak_id + one column per sample_id (row-aligned to peak_bed)
  - peak_bed          chrom, start, end, peak_id  (row-aligned)
  - sample_metadata_tsv  sample_id, cancer_type
Downloaded inputs live OUTSIDE the committed tree (inputs/NCI_GDC_TCGA_ATAC/). Thin data-prep.
"""
from __future__ import annotations
import argparse, csv
from pathlib import Path

def parse_args():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--corces_matrix_tsv",required=True,help="Pan-cancer normalized peak-by-sample matrix (TXT).")
    p.add_argument("--corces_peak_bed",required=True,help="Pan-cancer peak set BED.")
    p.add_argument("--sample_lookup_tsv",required=True,help="GDC ATAC sample lookup (barcode -> cancer type).")
    p.add_argument("--out_peak_matrix_tsv",required=True); p.add_argument("--out_peak_bed",required=True)
    p.add_argument("--out_sample_metadata_tsv",required=True)
    p.add_argument("--lookup_sample_column",default="sample_id"); p.add_argument("--lookup_cancer_type_column",default="cancer_type")
    return p.parse_args()

def main()->int:
    a=parse_args()
    # sample lookup -> {sample_id: cancer_type}
    with open(a.sample_lookup_tsv,encoding="utf-8",newline="") as h:
        look={}
        for r in csv.DictReader(h,delimiter="\t"):
            sid=str(r.get(a.lookup_sample_column,"")).strip(); ct=str(r.get(a.lookup_cancer_type_column,"")).strip()
            if sid and ct: look[sid]=ct
    if not look: raise SystemExit("Empty sample lookup.")
    # pass-through matrix (assume first column is peak id/coords; keep sample cols present in lookup)
    out_mat=Path(a.out_peak_matrix_tsv).expanduser().resolve(); out_mat.parent.mkdir(parents=True,exist_ok=True)
    kept_samples=[]
    with open(a.corces_matrix_tsv,encoding="utf-8",newline="") as h, out_mat.open("w",encoding="utf-8",newline="") as o:
        reader=csv.reader(h,delimiter="\t"); header=next(reader); w=csv.writer(o,delimiter="\t",lineterminator="\n")
        peak_idx=0
        keep_idx=[peak_idx]+[i for i,c in enumerate(header) if i!=peak_idx and c in look]
        kept_samples=[header[i] for i in keep_idx[1:]]
        w.writerow(["peak_id"]+kept_samples)
        for row in reader:
            if row: w.writerow([row[i] for i in keep_idx])
    if not kept_samples: raise SystemExit("No matrix sample columns matched the lookup.")
    # peak bed pass-through (ensure peak_id column)
    out_bed=Path(a.out_peak_bed).expanduser().resolve()
    with open(a.corces_peak_bed,encoding="utf-8",newline="") as h, out_bed.open("w",encoding="utf-8",newline="") as o:
        w=csv.writer(o,delimiter="\t",lineterminator="\n")
        for i,line in enumerate(h):
            parts=line.rstrip("\n").split("\t")
            if len(parts)<3: continue
            chrom,start,end=parts[0],parts[1],parts[2]
            pid=parts[3] if len(parts)>3 else f"peak{i+1}"
            w.writerow([chrom,start,end,pid])
    out_meta=Path(a.out_sample_metadata_tsv).expanduser().resolve()
    with out_meta.open("w",encoding="utf-8",newline="") as o:
        w=csv.writer(o,delimiter="\t",lineterminator="\n"); w.writerow(["sample_id","cancer_type"])
        for s in kept_samples: w.writerow([s,look[s]])
    print(f"Wrote {len(kept_samples)} samples -> {out_mat}",flush=True)
    return 0

if __name__=="__main__": raise SystemExit(main())
