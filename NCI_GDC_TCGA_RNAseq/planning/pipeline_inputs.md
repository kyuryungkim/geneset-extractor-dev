# NCI_GDC_TCGA_RNAseq — pipeline inputs & run flow

Mirrors the GTEx library. `geneset-extractor-dev` is a thin wrapper; all workflow logic
lives in `dig-gene-set-extractors` (DIG). Reference pattern = GTEx.

## Model
- **TR1** (`tumor_vs_rest`): for one focal TCGA project, contrast its **Primary Tumor**
  samples against the pooled Primary Tumor samples of all other TCGA projects.
- Output: directional gene sets `TCGA_<TUMOR>_tumor_vs_rest_up` / `_dn`.
- 33 projects × TR1 ⇒ ~66 directional gene sets.

## Inputs (acquired outside the committed tree, e.g. `inputs/NCI_GDC_TCGA_RNAseq/`)
1. **GDC open-access STAR-Counts** files: `*.rna_seq.augmented_star_gene_counts.tsv`
   (GRCh38 / GENCODE v36; `unstranded` raw-count column used for DE).
2. **GDC sample sheet** mapping file → Sample ID, Project ID, Sample Type.
3. **GTF** (GENCODE) for protein-coding biotype filtering.

`build_tcga_inputs.py` merges (1)+(2) into the two files the runner consumes:
- `counts.tsv`            gene_by_sample: `gene_id`, `gene_symbol`, `<sample_id...>`
- `sample_metadata.tsv`   `sample_id`, `project_id`, `sample_type` (Primary Tumor only by default)

## DIG calls (built by `run_tumor_vs_rest_model.py`)
1. `workflows rna_de_prepare --modality bulk --comparison_mode group_vs_rest
   --group_column project_id --condition_a <PROJECT> --de_mode modern --balance_groups true ...`
   → emits `workflow/deg_long.tsv` (+ provenance graph, comparison manifests).
2. `convert rna_deg_multi --deg_tsv workflow/deg_long.tsv --comparison_column comparison_id
   --comparison_name_column gmt_comparison_label --signature_name TCGA_<TUMOR> ...`
   → emits `extractor/` gene sets, `genesets.gmt`, `geneset.meta.json`, `geneset.provenance.json`,
   plus the wrapper-written `geneset.model.json` and per-group provenance rebuild.

**Note:** `de_mode=harmonizome` cannot be used with `group_vs_rest` (it requires explicit
two-group contrasts). Group-size balancing for tumor-vs-rest uses the `balance_groups` flag
under `de_mode=modern`.

## Output contract (per tumor type × model)
```
genesets/<tumor_type_id>/models/<model_id>/
  workflow/    deg_long.tsv, deg_long.provenance_graph.json, comparison_*.tsv, prepare_summary.json
  extractor/   manifest.tsv, genesets.gmt, <comparison>/{geneset.tsv, geneset.full.tsv,
               geneset.meta.json, geneset.model.json, geneset.provenance.json, run_summary.*}
```

## How to run
Single task (local): `NCI_GDC_TCGA_RNAseq/run/build_tcga_rnaseq_genesets.sh --tumor_types tcga_brca
--models TR1 --counts_tsv <...> --sample_metadata_tsv <...> --gtf <...> --dig_dir <...> --out_root <...>`

Cluster array (one task per tumor_type × model): `run/submit_tcga_rnaseq_models_cluster_apptainer.sh
--submit` (or the non-apptainer `submit_tcga_rnaseq_models_cluster.sh`).

Metadata/provenance refresh and S3 publish use the shared
`run/refresh_model_metadata_and_provenance.sh` and `run/publish_library_to_s3.sh` (no library-specific
exceptions; verified).
