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
#   <miRNA_str> \
#   <mRNA_str> \
#   <mRNA_str_up> \
#   <mRNA_str_down> \
#   <lncRNA_str> \
#   <circRNA_str> \
#   <outdir> \
#   <out_prefix> \
#   <cancer_type> \
#   <has_mRNA_direction>

if [ $# -lt 11 ]; then
    echo "Error: Missing arguments."
    echo "Usage: sbatch $0 <uuid> <miRNA_str> <mRNA_str> <mRNA_str_up> <mRNA_str_down> <lncRNA_str> <circRNA_str> <outdir> <out_prefix> <cancer_type> <has_mRNA_direction>"
    exit 1
fi

uuid="$1"
miRNA_str="$2"
mRNA_str="$3"
mRNA_str_up="$4"
mRNA_str_down="$5"
lncRNA_str="$6"
circRNA_str="$7"
outdir="$8"
out_prefix="$9"
cancer_type="${10}"
has_mRNA_direction="${11}"

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
echo "miRNA count string length: ${#miRNA_str}"
echo "mRNA count string length: ${#mRNA_str}"
echo "mRNA_str_up length: ${#mRNA_str_up}"
echo "mRNA_str_down length: ${#mRNA_str_down}"
echo "lncRNA count string length: ${#lncRNA_str}"
echo "circRNA count string length: ${#circRNA_str}"
echo "outdir: ${outdir}"
echo "out_prefix: ${out_prefix}"
echo "cancer_type: ${cancer_type}"
echo "has_mRNA_direction: ${has_mRNA_direction}"
echo "script_wdr: ${script_wdr}"
echo "========================================"

if [ ! -d "${script_wdr}" ]; then
    echo "Error: script_wdr does not exist: ${script_wdr}"
    write_status "fail"
    exit 1
fi

if [ ! -f "${script_wdr}/run/run_module1.sh" ]; then
    echo "Error: target script does not exist: ${script_wdr}/run/run_module1.sh"
    write_status "fail"
    exit 1
fi

if [ -z "${outdir}" ]; then
    echo "Error: outdir is empty."
    write_status "fail"
    exit 1
fi

if [ -z "${out_prefix}" ]; then
    echo "Error: out_prefix is empty."
    write_status "fail"
    exit 1
fi

if [ -z "${cancer_type}" ]; then
    echo "Error: cancer_type is empty."
    write_status "fail"
    exit 1
fi

if [ "${has_mRNA_direction}" != "True" ] && [ "${has_mRNA_direction}" != "False" ]; then
    echo "Error: has_mRNA_direction must be True or False."
    write_status "fail"
    exit 1
fi

echo "Running run_module1.sh..."

bash "${script_wdr}/run/run_module1.sh" \
    "${miRNA_str}" \
    "${mRNA_str}" \
    "${mRNA_str_up}" \
    "${mRNA_str_down}" \
    "${lncRNA_str}" \
    "${circRNA_str}" \
    "${outdir}" \
    "${out_prefix}" \
    "${cancer_type}" \
    "${has_mRNA_direction}"

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