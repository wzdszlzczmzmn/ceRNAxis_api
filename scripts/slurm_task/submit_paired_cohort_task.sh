#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=compute

set -Eeuo pipefail

# Usage:
# sbatch submit_paired_cohort_task.sh \
#   <uuid> \
#   <dataset> \
#   <mrna_file> \
#   <mirna_file> \
#   <lncrna_file> \
#   <meta_file> \
#   <outdir> \
#   <logfc_cutoff_mrna> \
#   <padj_cutoff_mrna> \
#   <logfc_cutoff_mirna> \
#   <padj_cutoff_mirna> \
#   <logfc_cutoff_lncrna> \
#   <padj_cutoff_lncrna> \
#   <deg_method> \
#   <map_info_csv>

if [ $# -lt 15 ]; then
    echo "Error: Missing arguments."
    echo "Usage: sbatch $0 <uuid> <dataset> <mrna_file> <mirna_file> <lncrna_file> <meta_file> <outdir> <logfc_cutoff_mrna> <padj_cutoff_mrna> <logfc_cutoff_mirna> <padj_cutoff_mirna> <logfc_cutoff_lncrna> <padj_cutoff_lncrna> <deg_method> <map_info_csv>"
    exit 1
fi

uuid="$1"
dataset="$2"
mrna_file="$3"
mirna_file="$4"
lncrna_file="$5"
meta_file="$6"
outdir="$7"
logfc_cutoff_mrna="$8"
padj_cutoff_mrna="$9"
logfc_cutoff_mirna="${10}"
padj_cutoff_mirna="${11}"
logfc_cutoff_lncrna="${12}"
padj_cutoff_lncrna="${13}"
deg_method="${14}"
map_info_csv="${15}"

# Fixed parameters for Module2.
expr_sample_col="sample_id"
meta_sample_col="sample_id"
group_col="c_group"
case_label="case"
control_label="control"

# 根据实际部署路径修改
script_wdr="/home/platform/workspace/ceRNAixDB"

status_file="${outdir}/status.txt"

mkdir -p "${outdir}"

write_status() {
    local task_status="$1"
    local finished_time

    finished_time=$(date +"%Y-%m-%d %H:%M:%S")

    echo "${finished_time}" > "${status_file}"
    echo "${task_status}" >> "${status_file}"
}

fail_task() {
    local message="$1"

    echo "Error: ${message}"
    write_status "fail"
    exit 1
}

echo "========================================"
echo "Starting paired cohort task"
echo "Start time: $(date +"%Y-%m-%d %H:%M:%S")"
echo "UUID: ${uuid}"
echo "dataset: ${dataset}"
echo "mrna_file: ${mrna_file}"
echo "mirna_file: ${mirna_file}"
echo "lncrna_file: ${lncrna_file}"
echo "meta_file: ${meta_file}"
echo "outdir: ${outdir}"
echo "expr_sample_col: ${expr_sample_col}"
echo "meta_sample_col: ${meta_sample_col}"
echo "group_col: ${group_col}"
echo "case_label: ${case_label}"
echo "control_label: ${control_label}"
echo "logfc_cutoff_mrna: ${logfc_cutoff_mrna}"
echo "padj_cutoff_mrna: ${padj_cutoff_mrna}"
echo "logfc_cutoff_mirna: ${logfc_cutoff_mirna}"
echo "padj_cutoff_mirna: ${padj_cutoff_mirna}"
echo "logfc_cutoff_lncrna: ${logfc_cutoff_lncrna}"
echo "padj_cutoff_lncrna: ${padj_cutoff_lncrna}"
echo "deg_method: ${deg_method}"
echo "map_info_csv: ${map_info_csv}"
echo "script_wdr: ${script_wdr}"
echo "========================================"

if [ ! -d "${script_wdr}" ]; then
    fail_task "script_wdr does not exist: ${script_wdr}"
fi

if [ ! -f "${script_wdr}/run/run_module2_all.sh" ]; then
    fail_task "target script does not exist: ${script_wdr}/run/run_module2_all.sh"
fi

if [ ! -f "${mrna_file}" ]; then
    fail_task "mrna_file does not exist: ${mrna_file}"
fi

if [ ! -f "${mirna_file}" ]; then
    fail_task "mirna_file does not exist: ${mirna_file}"
fi

if [ ! -f "${lncrna_file}" ]; then
    fail_task "lncrna_file does not exist: ${lncrna_file}"
fi

if [ ! -f "${meta_file}" ]; then
    fail_task "meta_file does not exist: ${meta_file}"
fi

if [ ! -f "${map_info_csv}" ]; then
    fail_task "map_info_csv does not exist: ${map_info_csv}"
fi

if [ "${deg_method}" != "limma" ] && [ "${deg_method}" != "deseq2" ]; then
    fail_task "Invalid deg_method: ${deg_method}. Allowed values: limma, deseq2."
fi

echo "Running run_module2_all.sh..."

bash "${script_wdr}/run/run_module2_all.sh" \
    "${dataset}" \
    "${mrna_file}" \
    "${mirna_file}" \
    "${lncrna_file}" \
    "${meta_file}" \
    "${outdir}" \
    "${expr_sample_col}" \
    "${meta_sample_col}" \
    "${group_col}" \
    "${case_label}" \
    "${control_label}" \
    "${logfc_cutoff_mrna}" \
    "${padj_cutoff_mrna}" \
    "${logfc_cutoff_mirna}" \
    "${padj_cutoff_mirna}" \
    "${logfc_cutoff_lncrna}" \
    "${padj_cutoff_lncrna}" \
    "${deg_method}" \
    "${map_info_csv}"

script_exit_code=$?
finished_time=$(date +"%Y-%m-%d %H:%M:%S")

if [ "${script_exit_code}" -ne 0 ]; then
    echo "Task failed at ${finished_time}"
    write_status "fail"
else
    echo "Task completed successfully at ${finished_time}"
    write_status "success"
fi

exit "${script_exit_code}"