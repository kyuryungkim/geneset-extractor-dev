#!/usr/bin/env bash
set -euo pipefail
# Submit TCGA methylation tumor-vs-normal models as a qsub array (host Python).
# Mirrors the RNAseq/CNV submit scripts. One family (MD*); only projects with matched normals.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
REPO_ROOT="${REPO_ROOT:-${DEFAULT_REPO_ROOT}}"
SELF_PATH="${REPO_ROOT}/geneset-extractor-dev/run/submit_tcga_meth_models_cluster.sh"
DIG_DIR="${DIG_DIR:-${REPO_ROOT}/dig-gene-set-extractors}"
WORK_ROOT="${WORK_ROOT:-${PWD}}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
CONFIG_ROOT="${TCGA_METH_CONFIG_ROOT:-${REPO_ROOT}/geneset-extractor-dev/NCI_GDC_TCGA_Methylation/config}"
MODEL_LIST="${TCGA_METH_MODEL_LIST:-${CONFIG_ROOT}/model_list.tsv}"
TUMOR_TYPE_LIST="${TCGA_METH_TUMOR_TYPE_LIST:-${CONFIG_ROOT}/tumor_type_list.tsv}"
MODEL_MANIFEST="${TCGA_METH_MODEL_MANIFEST:-${CONFIG_ROOT}/model_manifest.tsv}"
OUT_ROOT="${TCGA_METH_OUT_ROOT:-${WORK_ROOT}/tcga_meth_all_models}"
WORKLIST="${TCGA_METH_WORKLIST:-${WORK_ROOT}/tcga_meth_qsub_worklist.tsv}"
QSUB_LOG_ROOT="${QSUB_LOG_ROOT:-${WORK_ROOT}/qsub_logs_tcga_meth}"
ARRAY_MEMORY="${TCGA_METH_ARRAY_MEMORY:-32G}"; ARRAY_WALLTIME="${TCGA_METH_ARRAY_WALLTIME:-24:00:00}"
JOB_NAME="${TCGA_METH_JOB_NAME:-tcga_meth_models}"; QSUB_BIN="${QSUB_BIN:-qsub}"
APPTAINER_BIN="${APPTAINER_BIN:-apptainer}"; APPTAINER_IMAGE="${APPTAINER_IMAGE:-}"
APPTAINER_EXTRA_ARGS="${APPTAINER_EXTRA_ARGS:-}"; APPTAINER_PYTHON_BIN="${APPTAINER_PYTHON_BIN:-python}"
BETA_MATRIX_TSV="${TCGA_METH_BETA_MATRIX_TSV:-}"; SAMPLE_METADATA_TSV="${TCGA_METH_SAMPLE_METADATA_TSV:-}"
GTF="${TCGA_METH_GTF:-}"; PROBE_MANIFEST_TSV="${TCGA_METH_PROBE_MANIFEST_TSV:-}"
MODE="build"; FILTER_TUMOR_TYPE=""; FILTER_MODELS=""; OVERWRITE=""
usage(){ echo "Usage: $0 --submit [--write_model_only] [--tumor_type_id ID] [--model_id M[,M]] [--overwrite]"; }
require_file(){ [[ -f "$1" ]] || { echo "Missing required file: $1" >&2; exit 1; }; }
require_var(){ [[ -n "${!1:-}" ]] || { echo "Missing required environment variable: $1" >&2; exit 1; }; }
append_bind_path(){ [[ -z "${1:-}" ]] && return 0; [[ -d "$1" ]] && printf '%s\n' "$1" || dirname "$1"; }
write_worklist(){
  awk -F '\t' -v filter_tt="${FILTER_TUMOR_TYPE}" -v filter_models="${FILTER_MODELS}" -v model_list="${MODEL_LIST}" '
    BEGIN { n_models=0
      while ((getline line < model_list) > 0) { if (++ml==1) continue; split(line,f,"\t"); if (tolower(f[3])!="true") continue; models[++n_models]=f[1] }
      if (length(filter_models)>0){ split(filter_models,w,","); for(i in w) want[w[i]]=1; n_want=1 } task=0 }
    NR==1 { next }
    { tt=$1; proj=$2; label=$3; has_normal=tolower($4)
      if (length(filter_tt)>0 && tt!=filter_tt) next
      if (has_normal!="true") next
      for (mi=1; mi<=n_models; mi++){ m=models[mi]; if (n_want && !(m in want)) continue; printf "%d\t%s\t%s\t%s\t%s\n", ++task, tt, proj, label, m } }
  ' "${TUMOR_TYPE_LIST}" > "${WORKLIST}"
  [[ "$(awk 'END{print NR}' "${WORKLIST}")" -gt 0 ]] || { echo "Empty worklist" >&2; exit 1; }
}
run_task(){
  local task_id="${SGE_TASK_ID:-${SLURM_ARRAY_TASK_ID:-${1:-1}}}"
  local row tt proj label model
  row="$(awk -F '\t' -v t="${task_id}" '$1==t{print; exit}' "${WORKLIST}")"
  [[ -n "${row}" ]] || { echo "No worklist row for task ${task_id}" >&2; exit 1; }
  IFS=$'\t' read -r _ tt proj label model <<<"${row}"
  mkdir -p "${OUT_ROOT}"
  export PYTHONPATH="${DIG_DIR}/src"
  echo "[task ${task_id}] ${tt} / ${model}"
  if [[ "${MODE}" == "write_model_only" ]]; then
    exec "${PYTHON_BIN}" "${REPO_ROOT}/geneset-extractor-dev/NCI_GDC_TCGA_Methylation/src/run_methylation_diff_model.py" \
      --model_id "${model}" --tumor_type_id "${tt}" --tumor_type_label "${label}" --project_id "${proj}" \
      --run_root "${OUT_ROOT}/genesets/${tt}/models" --dig_dir "${DIG_DIR}" --model_manifest "${MODEL_MANIFEST}" --write_model_only
  fi
  local -a cmd
  cmd=("${PYTHON_BIN}" "${REPO_ROOT}/geneset-extractor-dev/NCI_GDC_TCGA_Methylation/src/build_tcga_meth_genesets.py"
       --tumor_types "${tt}" --models "${model}" --beta_matrix_tsv "${BETA_MATRIX_TSV}"
       --sample_metadata_tsv "${SAMPLE_METADATA_TSV}" --gtf "${GTF}" --dig_dir "${DIG_DIR}"
       --out_root "${OUT_ROOT}" --model_manifest "${MODEL_MANIFEST}" --model_list "${MODEL_LIST}" --tumor_type_list "${TUMOR_TYPE_LIST}")
  [[ -n "${PROBE_MANIFEST_TSV}" ]] && cmd+=(--probe_manifest_tsv "${PROBE_MANIFEST_TSV}")
  [[ -n "${OVERWRITE}" ]] && cmd+=(--overwrite)
  exec "${cmd[@]}"
}
do_submit(){

  require_file "${MODEL_LIST}"; require_file "${TUMOR_TYPE_LIST}"; require_file "${MODEL_MANIFEST}"
  if [[ "${MODE}" == "build" ]]; then
    require_var TCGA_METH_BETA_MATRIX_TSV; require_var TCGA_METH_SAMPLE_METADATA_TSV; require_var TCGA_METH_GTF
    require_file "${BETA_MATRIX_TSV}"; require_file "${SAMPLE_METADATA_TSV}"; require_file "${GTF}"
  fi
  write_worklist
  local tasks; tasks="$(awk 'END{print NR}' "${WORKLIST}")"; mkdir -p "${QSUB_LOG_ROOT}"
  echo "TCGA methylation worklist: ${WORKLIST} (${tasks} tasks), mode=${MODE}"
  "${QSUB_BIN}" -t "1-${tasks}" -N "${JOB_NAME}" -o "${QSUB_LOG_ROOT}" -e "${QSUB_LOG_ROOT}" -j y \
    -l "h_vmem=${ARRAY_MEMORY},h_rt=${ARRAY_WALLTIME}" \
    -v "REPO_ROOT=${REPO_ROOT},DIG_DIR=${DIG_DIR},WORK_ROOT=${WORK_ROOT},PYTHON_BIN=${PYTHON_BIN},APPTAINER_IMAGE=${APPTAINER_IMAGE:-},TCGA_METH_OUT_ROOT=${OUT_ROOT},TCGA_METH_WORKLIST=${WORKLIST},TCGA_METH_BETA_MATRIX_TSV=${BETA_MATRIX_TSV},TCGA_METH_SAMPLE_METADATA_TSV=${SAMPLE_METADATA_TSV},TCGA_METH_GTF=${GTF},TCGA_METH_PROBE_MANIFEST_TSV=${PROBE_MANIFEST_TSV},TCGA_TASK_MODE=${MODE},TCGA_TASK_OVERWRITE=${OVERWRITE}" \
    "${SELF_PATH}" --run_task
}
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
MODE="${TCGA_TASK_MODE:-${MODE}}"; OVERWRITE="${TCGA_TASK_OVERWRITE:-${OVERWRITE}}"
case "${ACTION}" in submit) do_submit ;; run_task) run_task ;; *) echo "Specify --submit" >&2; usage >&2; exit 1 ;; esac
