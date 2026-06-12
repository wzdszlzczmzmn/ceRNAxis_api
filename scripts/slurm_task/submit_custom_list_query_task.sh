#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --partition=compute

set -Eeuo pipefail

# Usage:
# sbatch submit_custom_list_query_task.sh \
#   <uuid> \
#   <map_info_csv> \
#   <cerna_axis_csv> \
#   <outdir> \
#   <out_prefix>

if [ $# -lt 5 ]; then
    echo "Error: Missing arguments."
    echo "Usage: sbatch $0 <uuid> <map_info_csv> <cerna_axis_csv> <outdir> <out_prefix>"
    exit 1
fi

uuid="$1"
map_info_csv="$2"
cerna_axis_csv="$3"
outdir="$4"
out_prefix="$5"

# 根据实际部署路径修改
script_wdr="/home/platform/workspace/ceRNAixDB"

status_file="${outdir}/status.txt"

mkdir -p "${outdir}"

write_status() {
    local status="$1"
    local finished_time

    finished_time=$(date +"%Y-%m-%d %H:%M:%S")

    echo "${finished_time}" > "${status_file}"
    echo "${status}" >> "${status_file}"
}

echo "========================================"
echo "Starting custom list query task"
echo "Start time: $(date +"%Y-%m-%d %H:%M:%S")"
echo "UUID: ${uuid}"
echo "map_info_csv: ${map_info_csv}"
echo "cerna_axis_csv: ${cerna_axis_csv}"
echo "outdir: ${outdir}"
echo "out_prefix: ${out_prefix}"
echo "script_wdr: ${script_wdr}"
echo "========================================"

if [ ! -d "${script_wdr}" ]; then
    echo "Error: script_wdr does not exist: ${script_wdr}"
    write_status "fail"
    exit 1
fi

if [ ! -f "${script_wdr}/run/run_map_immune_gene_axis.sh" ]; then
    echo "Error: target script does not exist: ${script_wdr}/run/run_map_immune_gene_axis.sh"
    write_status "fail"
    exit 1
fi

if [ ! -f "${map_info_csv}" ]; then
    echo "Error: map_info_csv does not exist: ${map_info_csv}"
    write_status "fail"
    exit 1
fi

if [ ! -f "${cerna_axis_csv}" ]; then
    echo "Error: cerna_axis_csv does not exist: ${cerna_axis_csv}"
    write_status "fail"
    exit 1
fi

echo "Running run_map_immune_gene_axis.sh..."

bash "${script_wdr}/run/run_map_immune_gene_axis.sh" \
    "${map_info_csv}" \
    "${cerna_axis_csv}" \
    "${outdir}" \
    "${out_prefix}"

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