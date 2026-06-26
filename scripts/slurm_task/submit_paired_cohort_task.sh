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
#   <circrna_file> \
#   <meta_file> \
#   <outdir> \
#   <logfc_cutoff_mrna> \
#   <padj_cutoff_mrna> \
#   <logfc_cutoff_mirna> \
#   <padj_cutoff_mirna> \
#   <logfc_cutoff_lncrna> \
#   <padj_cutoff_lncrna> \
#   <logfc_cutoff_circrna> \
#   <padj_cutoff_circrna> \
#   <deg_method> \
#   <map_info_csv> \
#   <cancer_type> \
#   <use_padj>
#
# Notes:
#   - mrna_file, mirna_file and meta_file are required.
#   - At least one of lncrna_file or circrna_file is required.
#   - cancer_type can be an empty string.
#   - use_padj must be TRUE or FALSE.

if [ $# -lt 20 ]; then
    echo "Error: Missing arguments."
    echo "Usage: sbatch $0 <uuid> <dataset> <mrna_file> <mirna_file> <lncrna_file> <circrna_file> <meta_file> <outdir> <logfc_cutoff_mrna> <padj_cutoff_mrna> <logfc_cutoff_mirna> <padj_cutoff_mirna> <logfc_cutoff_lncrna> <padj_cutoff_lncrna> <logfc_cutoff_circrna> <padj_cutoff_circrna> <deg_method> <map_info_csv> <cancer_type> <use_padj>"
    exit 1
fi

uuid="$1"
dataset="$2"
mrna_file="$3"
mirna_file="$4"
lncrna_file="$5"
circrna_file="$6"
meta_file="$7"
outdir="$8"

logfc_cutoff_mrna="$9"
padj_cutoff_mrna="${10}"
logfc_cutoff_mirna="${11}"
padj_cutoff_mirna="${12}"
logfc_cutoff_lncrna="${13}"
padj_cutoff_lncrna="${14}"
logfc_cutoff_circrna="${15}"
padj_cutoff_circrna="${16}"
deg_method="${17}"
map_info_csv="${18}"
cancer_type="${19}"
use_padj="${20}"

# Fixed parameters for Module2.
expr_sample_col="sample_id"
meta_sample_col="sample_id"
group_col="c_group"
case_label="case"
control_label="control"

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
echo "circrna_file: ${circrna_file}"
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
echo "logfc_cutoff_circrna: ${logfc_cutoff_circrna}"
echo "padj_cutoff_circrna: ${padj_cutoff_circrna}"
echo "deg_method: ${deg_method}"
echo "map_info_csv: ${map_info_csv}"
echo "cancer_type: ${cancer_type}"
echo "use_padj: ${use_padj}"
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

has_lncrna_file="FALSE"
has_circrna_file="FALSE"

if [ -n "${lncrna_file}" ] && [ -f "${lncrna_file}" ]; then
    has_lncrna_file="TRUE"
fi

if [ -n "${circrna_file}" ] && [ -f "${circrna_file}" ]; then
    has_circrna_file="TRUE"
fi

if [ "${has_lncrna_file}" != "TRUE" ] && [ "${has_circrna_file}" != "TRUE" ]; then
    fail_task "At least one of lncrna_file or circrna_file must exist."
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

if [ "${use_padj}" != "TRUE" ] && [ "${use_padj}" != "FALSE" ]; then
    fail_task "Invalid use_padj: ${use_padj}. Allowed values: TRUE, FALSE."
fi

echo "Input file availability:"
echo "has_lncrna_file: ${has_lncrna_file}"
echo "has_circrna_file: ${has_circrna_file}"

echo "Running run_module2_all.sh..."

bash "${script_wdr}/run/run_module2_all.sh" \
    "${dataset}" \
    "${mrna_file}" \
    "${mirna_file}" \
    "${lncrna_file}" \
    "${circrna_file}" \
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
    "${logfc_cutoff_circrna}" \
    "${padj_cutoff_circrna}" \
    "${deg_method}" \
    "${map_info_csv}" \
    "${cancer_type}" \
    "${use_padj}"

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