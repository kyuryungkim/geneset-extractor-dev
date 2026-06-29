#!/usr/bin/env bash
set -euo pipefail

# Submit TCGA RNA-seq tumor-vs-rest models as a qsub array, each task running inside
# Apptainer. Mirrors run/submit_gtex_models_cluster_apptainer.sh (simplified: one model
# family, one data version, partition = TCGA tumor type).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
REPO_ROOT="${REPO_ROOT:-${DEFAULT_REPO_ROOT}}"
SELF_PATH="${REPO_ROOT}/geneset-extractor-dev/run/submit_tcga_rnaseq_models_cluster_apptainer.sh"
DIG_DIR="${DIG_DIR:-${REPO_ROOT}/dig-gene-set-extractors}"
WORK_ROOT="${WORK_ROOT:-${PWD}}"

CONFIG_ROOT="${TCGA_RNASEQ_CONFIG_ROOT:-${REPO_ROOT}/geneset-extractor-dev/NCI_GDC_TCGA_RNAseq/config}"
MODEL_LIST="${TCGA_RNASEQ_MODEL_LIST:-${CONFIG_ROOT}/model_list.tsv}"
TUMOR_TYPE_LIST="${TCGA_RNASEQ_TUMOR_TYPE_LIST:-${CONFIG_ROOT}/tumor_type_list.tsv}"
MODEL_MANIFEST="${TCGA_RNASEQ_MODEL_MANIFEST:-${CONFIG_ROOT}/model_manifest.tsv}"

OUT_ROOT="${TCGA_RNASEQ_OUT_ROOT:-${WORK_ROOT}/tcga_rnaseq_all_models}"
WORKLIST="${TCGA_RNASEQ_WORKLIST:-${WORK_ROOT}/tcga_rnaseq_qsub_worklist.tsv}"
QSUB_LOG_ROOT="${QSUB_LOG_ROOT:-${WORK_ROOT}/qsub_logs_tcga_rnaseq}"
ARRAY_MEMORY="${TCGA_RNASEQ_ARRAY_MEMORY:-32G}"
ARRAY_WALLTIME="${TCGA_RNASEQ_ARRAY_WALLTIME:-24:00:00}"
JOB_NAME="${TCGA_RNASEQ_JOB_NAME:-tcga_rnaseq_models_apptainer}"

QSUB_BIN="${QSUB_BIN:-qsub}"
APPTAINER_BIN="${APPTAINER_BIN:-apptainer}"
APPTAINER_IMAGE="${APPTAINER_IMAGE:-}"
APPTAINER_EXTRA_ARGS="${APPTAINER_EXTRA_ARGS:-}"
APPTAINER_PYTHON_BIN="${APPTAINER_PYTHON_BIN:-python}"

# Data inputs (required for full runs; not for --write_model_only).
COUNTS_TSV="${TCGA_RNASEQ_COUNTS_TSV:-}"
SAMPLE_METADATA_TSV="${TCGA_RNASEQ_SAMPLE_METADATA_TSV:-}"
GTF="${TCGA_RNASEQ_GTF:-}"

MODE="build"           # build | write_model_only
FILTER_TUMOR_TYPE=""
FILTER_MODELS=""
OVERWRITE=""

usage() {
  cat <<'EOF'
Usage:
  ./geneset-extractor-dev/run/submit_tcga_rnaseq_models_cluster_apptainer.sh --submit \
      [--write_model_only] [--tumor_type_id ID] [--model_id MODEL[,MODEL...]] [--overwrite]
  ./geneset-extractor-dev/run/submit_tcga_rnaseq_models_cluster_apptainer.sh --help

Required env:
  APPTAINER_IMAGE
  For full runs also: TCGA_RNASEQ_COUNTS_TSV, TCGA_RNASEQ_SAMPLE_METADATA_TSV, TCGA_RNASEQ_GTF
Optional env:
  REPO_ROOT, DIG_DIR, WORK_ROOT, TCGA_RNASEQ_OUT_ROOT, TCGA_RNASEQ_WORKLIST,
  TCGA_RNASEQ_ARRAY_MEMORY, TCGA_RNASEQ_ARRAY_WALLTIME, QSUB_BIN

One array task per (tumor_type x model). Metadata/provenance refresh and S3 publish use the
shared scripts run/refresh_model_metadata_and_provenance.sh and run/publish_library_to_s3.sh.
EOF
}

require_file() { [[ -f "$1" ]] || { echo "Missing required file: $1" >&2; exit 1; }; }
require_var()  { [[ -n "${!1:-}" ]] || { echo "Missing required environment variable: $1" >&2; exit 1; }; }
append_bind_path() { [[ -z "${1:-}" ]] && return 0; [[ -d "$1" ]] && printf '%s\n' "$1" || dirname "$1"; }

write_worklist() {
  # Columns: task_id, tumor_type_id, project_id, tumor_type_label, model_id
  awk -F '\t' -v filter_tt="${FILTER_TUMOR_TYPE}" -v filter_models="${FILTER_MODELS}" \
      -v model_list="${MODEL_LIST}" '
    BEGIN {
      n_models = 0
      while ((getline line < model_list) > 0) {
        if (++ml == 1) continue                       # header
        split(line, f, "\t")
        if (tolower(f[3]) != "true") continue          # enabled column
        models[++n_models] = f[1]
      }
      n_want = 0
      if (length(filter_models) > 0) { split(filter_models, w, ","); for (i in w) want[w[i]] = 1; n_want = 1 }
      task = 0
    }
    NR == 1 { next }                                   # tumor_type_list header
    {
      tt = $1; proj = $2; label = $3
      if (length(filter_tt) > 0 && tt != filter_tt) next
      for (mi = 1; mi <= n_models; mi++) {
        m = models[mi]
        if (n_want && !(m in want)) continue
        printf "%d\t%s\t%s\t%s\t%s\n", ++task, tt, proj, label, m
      }
    }
  ' "${TUMOR_TYPE_LIST}" > "${WORKLIST}"
  if [[ "$(awk 'END { print NR }' "${WORKLIST}")" -le 0 ]]; then
    echo "TCGA RNA-seq filters produced an empty worklist" >&2
    exit 1
  fi
}

apptainer_bind_csv() {
  {
    printf '%s\n' "${REPO_ROOT}" "${DIG_DIR}" "${OUT_ROOT}"
    append_bind_path "${CONFIG_ROOT}"
    append_bind_path "${COUNTS_TSV}"
    append_bind_path "${SAMPLE_METADATA_TSV}"
    append_bind_path "${GTF}"
  } | awk 'NF && !seen[$0]++' | paste -sd, -
}

run_task() {
  local task_id="${SGE_TASK_ID:-${SLURM_ARRAY_TASK_ID:-${1:-1}}}"
  local row tt proj label model
  row="$(awk -F '\t' -v t="${task_id}" '$1 == t { print; exit }' "${WORKLIST}")"
  [[ -n "${row}" ]] || { echo "No worklist row for task ${task_id}" >&2; exit 1; }
  IFS=$'\t' read -r _ tt proj label model <<<"${row}"
  mkdir -p "${OUT_ROOT}"

  local bind_csv; bind_csv="$(apptainer_bind_csv)"
  local inner
  if [[ "${MODE}" == "write_model_only" ]]; then
    inner="'${APPTAINER_PYTHON_BIN}' '${REPO_ROOT}/geneset-extractor-dev/NCI_GDC_TCGA_RNAseq/src/run_tumor_vs_rest_model.py'"
    inner+=" --model_id '${model}' --tumor_type_id '${tt}' --tumor_type_label '${label}' --project_id '${proj}'"
    inner+=" --run_root '${OUT_ROOT}/genesets/${tt}/models' --dig_dir '${DIG_DIR}' --model_manifest '${MODEL_MANIFEST}' --write_model_only"
  else
    inner="'${APPTAINER_PYTHON_BIN}' '${REPO_ROOT}/geneset-extractor-dev/NCI_GDC_TCGA_RNAseq/src/build_tcga_rnaseq_genesets.py'"
    inner+=" --tumor_types '${tt}' --models '${model}' --counts_tsv '${COUNTS_TSV}' --sample_metadata_tsv '${SAMPLE_METADATA_TSV}'"
    inner+=" --dig_dir '${DIG_DIR}' --out_root '${OUT_ROOT}' --model_manifest '${MODEL_MANIFEST}'"
    inner+=" --model_list '${MODEL_LIST}' --tumor_type_list '${TUMOR_TYPE_LIST}'"
    [[ -n "${GTF}" ]] && inner+=" --gtf '${GTF}'"
    [[ -n "${OVERWRITE}" ]] && inner+=" --overwrite"
  fi

  declare -a EXEC_CMD
  EXEC_CMD=("${APPTAINER_BIN}" exec --bind "${bind_csv}")
  if [[ -n "${APPTAINER_EXTRA_ARGS}" ]]; then
    # shellcheck disable=SC2206
    EXTRA_ARGS=( ${APPTAINER_EXTRA_ARGS} ); EXEC_CMD+=("${EXTRA_ARGS[@]}")
  fi
  EXEC_CMD+=("${APPTAINER_IMAGE}" bash --noprofile --norc -c "export PYTHONPATH='${DIG_DIR}/src'; ${inner}")
  echo "[task ${task_id}] ${tt} / ${model}"
  exec "${EXEC_CMD[@]}"
}

do_submit() {
  require_var APPTAINER_IMAGE; require_file "${APPTAINER_IMAGE}"
  require_file "${MODEL_LIST}"; require_file "${TUMOR_TYPE_LIST}"; require_file "${MODEL_MANIFEST}"
  if [[ "${MODE}" == "build" ]]; then
    require_var TCGA_RNASEQ_COUNTS_TSV; require_var TCGA_RNASEQ_SAMPLE_METADATA_TSV
    require_file "${COUNTS_TSV}"; require_file "${SAMPLE_METADATA_TSV}"
    [[ -n "${GTF}" ]] && require_file "${GTF}"
  fi
  write_worklist
  local tasks; tasks="$(awk 'END { print NR }' "${WORKLIST}")"
  mkdir -p "${QSUB_LOG_ROOT}"
  echo "TCGA RNA-seq worklist: ${WORKLIST} (${tasks} tasks), mode=${MODE}"
  "${QSUB_BIN}" -t "1-${tasks}" -N "${JOB_NAME}" \
    -o "${QSUB_LOG_ROOT}" -e "${QSUB_LOG_ROOT}" -j y \
    -l "h_vmem=${ARRAY_MEMORY},h_rt=${ARRAY_WALLTIME}" \
    -v "REPO_ROOT=${REPO_ROOT},DIG_DIR=${DIG_DIR},WORK_ROOT=${WORK_ROOT},APPTAINER_IMAGE=${APPTAINER_IMAGE},TCGA_RNASEQ_OUT_ROOT=${OUT_ROOT},TCGA_RNASEQ_WORKLIST=${WORKLIST},TCGA_RNASEQ_COUNTS_TSV=${COUNTS_TSV},TCGA_RNASEQ_SAMPLE_METADATA_TSV=${SAMPLE_METADATA_TSV},TCGA_RNASEQ_GTF=${GTF},TCGA_TASK_MODE=${MODE},TCGA_TASK_OVERWRITE=${OVERWRITE}" \
    "${SELF_PATH}" --run_task
}

# ---- arg parsing ----
[[ $# -ge 1 ]] || { usage >&2; exit 1; }
ACTION=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help|help) usage; exit 0 ;;
    --submit) ACTION="submit"; shift ;;
    --run_task) ACTION="run_task"; shift ;;
    --write_model_only) MODE="write_model_only"; shift ;;
    --overwrite) OVERWRITE="1"; shift ;;
    --tumor_type_id) FILTER_TUMOR_TYPE="$2"; shift 2 ;;
    --model_id) FILTER_MODELS="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

# Allow env-passed task mode/overwrite (set by do_submit for array tasks).
MODE="${TCGA_TASK_MODE:-${MODE}}"
OVERWRITE="${TCGA_TASK_OVERWRITE:-${OVERWRITE}}"

case "${ACTION}" in
  submit) do_submit ;;
  run_task) run_task ;;
  *) echo "Specify --submit" >&2; usage >&2; exit 1 ;;
esac
