import subprocess

from analysis.models import CustomListQueryTask, PairedCohortTask
from analysis.utils.custom_list_query_task_utils import validate_task_mirna_axis_file, get_task_output_dir, \
    get_task_out_prefix
from analysis.utils.immune_annotation_path_utils import validate_immune_annotation_file
from analysis.utils.paired_cohort_task_utils import validate_paired_cohort_input_files, \
    get_paired_cohort_task_output_dir
from analysis.utils.slurm_path_utils import validate_slurm_script, get_pipeline_stdout_file_path, \
    get_pipeline_stderr_file_path, get_slurm_job_name

CUSTOM_LIST_QUERY_SLURM_SCRIPT = "submit_custom_list_query_task.sh"
PAIRED_COHORT_SLURM_SCRIPT = "submit_paired_cohort_task.sh"


def sbatch_custom_list_query_task(task_uuid) -> dict:
    """
    Submit CustomListQueryTask to SLURM.

    Expected SLURM script arguments:
        1. uuid
        2. map_info_csv
        3. cerna_axis_csv
        4. outdir
        5. out_prefix
    """

    try:
        task = CustomListQueryTask.objects.get(uuid=task_uuid)
    except CustomListQueryTask.DoesNotExist:
        return {
            "success": False,
            "msg": f"CustomListQueryTask not found: {task_uuid}",
            "stdout": "",
            "stderr": "",
        }

    try:
        script_path = validate_slurm_script(
            CUSTOM_LIST_QUERY_SLURM_SCRIPT
        )

        map_info_file = validate_immune_annotation_file(
            task.map_info
        )

        cerna_axis_file = validate_task_mirna_axis_file(
            task
        )

        output_dir = get_task_output_dir(task)
        output_dir.mkdir(parents=True, exist_ok=True)

        stdout_file = get_pipeline_stdout_file_path(output_dir)
        stderr_file = get_pipeline_stderr_file_path(output_dir)

        out_prefix = get_task_out_prefix(task)

    except Exception as e:
        return {
            "success": False,
            "msg": str(e),
            "stdout": "",
            "stderr": "",
        }

    command = [
        "sbatch",
        f"--job-name={get_slurm_job_name(task.uuid)}",
        f"--output={stdout_file}",
        f"--error={stderr_file}",
        str(script_path),
        str(task.uuid),
        str(map_info_file),
        str(cerna_axis_file),
        str(output_dir),
        out_prefix,
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        return {
            "success": False,
            "msg": "Failed to submit SLURM job.",
            "stdout": result.stdout,
            "stderr": result.stderr,
            "command": command,
        }

    return {
        "success": True,
        "msg": "SLURM job submitted successfully.",
        "stdout": result.stdout,
        "stderr": result.stderr,
        "command": command,
    }


def sbatch_paired_cohort_task(task_uuid) -> dict:
    """
    Submit PairedCohortTask to SLURM.

    Expected SLURM script arguments:
        1. uuid
        2. dataset
        3. mrna_file
        4. mirna_file
        5. lncrna_file
        6. meta_file
        7. outdir
        8. logfc_cutoff_mrna
        9. padj_cutoff_mrna
        10. logfc_cutoff_mirna
        11. padj_cutoff_mirna
        12. logfc_cutoff_lncrna
        13. padj_cutoff_lncrna
        14. deg_method
        15. map_info_csv
    """

    try:
        task = PairedCohortTask.objects.get(uuid=task_uuid)
    except PairedCohortTask.DoesNotExist:
        return {
            "success": False,
            "msg": f"PairedCohortTask not found: {task_uuid}",
            "stdout": "",
            "stderr": "",
        }

    try:
        script_path = validate_slurm_script(
            PAIRED_COHORT_SLURM_SCRIPT
        )

        map_info_file = validate_immune_annotation_file(
            task.map_info
        )

        input_files = validate_paired_cohort_input_files(
            task
        )

        output_dir = get_paired_cohort_task_output_dir(task)
        output_dir.mkdir(parents=True, exist_ok=True)

        stdout_file = get_pipeline_stdout_file_path(output_dir)
        stderr_file = get_pipeline_stderr_file_path(output_dir)

        dataset = str(task.task_name).strip()

        if not dataset:
            raise ValueError("Missing task_name for paired cohort task.")

        if task.deg_method not in ["limma", "deseq2"]:
            raise ValueError(
                "Invalid deg_method. Allowed values are: limma, deseq2."
            )

    except Exception as e:
        return {
            "success": False,
            "msg": str(e),
            "stdout": "",
            "stderr": "",
        }

    command = [
        "sbatch",
        f"--job-name={get_slurm_job_name(task.uuid)}",
        f"--output={stdout_file}",
        f"--error={stderr_file}",
        str(script_path),
        str(task.uuid),
        dataset,
        str(input_files["mrna_file"]),
        str(input_files["mirna_file"]),
        str(input_files["lncrna_file"]),
        str(input_files["circrna_file"]),
        str(input_files["meta_file"]),
        str(output_dir),
        str(task.logfc_cutoff_mrna),
        str(task.padj_cutoff_mrna),
        str(task.logfc_cutoff_mirna),
        str(task.padj_cutoff_mirna),
        str(task.logfc_cutoff_lncrna),
        str(task.padj_cutoff_lncrna),
        str(task.logfc_cutoff_circrna),
        str(task.padj_cutoff_circrna),
        task.deg_method,
        str(map_info_file),
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        return {
            "success": False,
            "msg": "Failed to submit SLURM job.",
            "stdout": result.stdout,
            "stderr": result.stderr,
            "command": command,
        }

    return {
        "success": True,
        "msg": "SLURM job submitted successfully.",
        "stdout": result.stdout,
        "stderr": result.stderr,
        "command": command,
    }
