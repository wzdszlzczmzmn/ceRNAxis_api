#!/bin/bash -l

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=compute

set -Eeuo pipefail


# Usage:
# sbatch submit_scst_hybrid_reference_task.sh \
#   <uuid> \
#   <dataset> \
#   <data_type> \
#   <exp_file> \
#   <meta_file> \
#   <tcga_type> \
#   <lncrna_type> \
#   <outdir> \
#   <group_col> \
#   <logfc_cutoff_mrna> \
#   <padj_cutoff_mrna> \
#   <map_info_csv> \
#   <use_padj>


if [ "$#" -ne 13 ]; then
    echo "Error: Invalid number of arguments."
    echo
    echo "Usage:"
    echo "  sbatch $0 \\"
    echo "    <uuid> \\"
    echo "    <dataset> \\"
    echo "    <data_type> \\"
    echo "    <exp_file> \\"
    echo "    <meta_file> \\"
    echo "    <tcga_type> \\"
    echo "    <lncrna_type> \\"
    echo "    <outdir> \\"
    echo "    <group_col> \\"
    echo "    <logfc_cutoff_mrna> \\"
    echo "    <padj_cutoff_mrna> \\"
    echo "    <map_info_csv> \\"
    echo "    <use_padj>"
    exit 1
fi


uuid="$1"
dataset="$2"
data_type="$3"
exp_file="$4"
meta_file="$5"
tcga_type="$6"
lncrna_type="$7"
outdir="$8"
group_col="$9"
logfc_cutoff_mrna="${10}"
padj_cutoff_mrna="${11}"
map_info_csv="${12}"
use_padj="${13}"


# Adjust this path to the actual deployment directory.
script_wdr="/home/platform/workspace/ceRNAixDB"

target_script="${script_wdr}/run/run_sc_st.sh"
status_file="${outdir}/status.txt"


mkdir -p "${outdir}"


write_status() {
    local task_status="$1"
    local finished_time

    finished_time=$(date +"%Y-%m-%d %H:%M:%S")

    {
        echo "${finished_time}"
        echo "${task_status}"
    } > "${status_file}"
}


fail_task() {
    local message="$1"

    echo "Error: ${message}" >&2
    write_status "fail"
    exit 1
}


handle_unexpected_error() {
    local exit_code=$?
    local line_number="$1"

    echo "Unexpected error at line ${line_number}, exit code ${exit_code}." >&2
    write_status "fail"

    exit "${exit_code}"
}


trap 'handle_unexpected_error ${LINENO}' ERR


echo "========================================"
echo "Starting SC/ST Hybrid Reference task"
echo "Start time: $(date +"%Y-%m-%d %H:%M:%S")"
echo "UUID: ${uuid}"
echo "dataset: ${dataset}"
echo "data_type: ${data_type}"
echo "exp_file: ${exp_file}"
echo "meta_file: ${meta_file}"
echo "tcga_type: ${tcga_type}"
echo "lncrna_type: ${lncrna_type}"
echo "outdir: ${outdir}"
echo "group_col: ${group_col}"
echo "logfc_cutoff_mrna: ${logfc_cutoff_mrna}"
echo "padj_cutoff_mrna: ${padj_cutoff_mrna}"
echo "map_info_csv: ${map_info_csv}"
echo "use_padj: ${use_padj}"
echo "script_wdr: ${script_wdr}"
echo "target_script: ${target_script}"
echo "========================================"


if [ ! -d "${script_wdr}" ]; then
    fail_task "script_wdr does not exist: ${script_wdr}"
fi


if [ ! -f "${target_script}" ]; then
    fail_task "target script does not exist: ${target_script}"
fi


if [ ! -f "${exp_file}" ]; then
    fail_task "exp_file does not exist: ${exp_file}"
fi


if [ ! -f "${meta_file}" ]; then
    fail_task "meta_file does not exist: ${meta_file}"
fi


if [ ! -f "${map_info_csv}" ]; then
    fail_task "map_info_csv does not exist: ${map_info_csv}"
fi


case "${data_type}" in
    sc)
        expected_id_column="cell_id"
        ;;
    st)
        expected_id_column="spot_id"
        ;;
    *)
        fail_task (
            "Invalid data_type: ${data_type}. "
            "Allowed values: sc, st."
        )
        ;;
esac


if [ "${use_padj}" != "TRUE" ] &&
   [ "${use_padj}" != "FALSE" ]; then
    fail_task (
        "Invalid use_padj: ${use_padj}. "
        "Allowed values: TRUE, FALSE."
    )
fi


case "${exp_file,,}" in
    *.parquet)
        ;;
    *)
        fail_task (
            "Invalid exp_file extension: ${exp_file}. "
            "Expected a .parquet file."
        )
        ;;
esac


case "${meta_file,,}" in
    *.csv)
        ;;
    *)
        fail_task (
            "Invalid meta_file extension: ${meta_file}. "
            "Expected a .csv file."
        )
        ;;
esac


if [ -z "${dataset}" ]; then
    fail_task "dataset cannot be empty."
fi


if [ -z "${group_col}" ]; then
    fail_task "group_col cannot be empty."
fi


echo "Expected identifier column/index: ${expected_id_column}"
echo "Running run_sc_st.sh..."


set +e

bash "${target_script}" \
    "${dataset}" \
    "${exp_file}" \
    "${meta_file}" \
    "${tcga_type}" \
    "${lncrna_type}" \
    "${outdir}" \
    "${group_col}" \
    "${logfc_cutoff_mrna}" \
    "${padj_cutoff_mrna}" \
    "${map_info_csv}" \
    "${use_padj}"

script_exit_code=$?

set -e


finished_time=$(date +"%Y-%m-%d %H:%M:%S")


if [ "${script_exit_code}" -ne 0 ]; then
    echo (
        "SC/ST Hybrid Reference task failed at "
        "${finished_time}, exit code: ${script_exit_code}"
    ) >&2

    write_status "fail"
else
    echo (
        "SC/ST Hybrid Reference task completed successfully "
        "at ${finished_time}"
    )

    write_status "success"
fi


exit "${script_exit_code}"