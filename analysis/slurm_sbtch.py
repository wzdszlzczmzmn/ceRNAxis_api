import subprocess

from analysis.models import CustomListQueryTask
from analysis.utils.custom_list_query_task_utils import validate_task_mirna_axis_file, get_task_output_dir, \
    get_task_out_prefix
from analysis.utils.immune_annotation_path_utils import validate_immune_annotation_file
from analysis.utils.slurm_path_utils import validate_slurm_script, get_pipeline_stdout_file_path, \
    get_pipeline_stderr_file_path, get_slurm_job_name

CUSTOM_LIST_QUERY_SLURM_SCRIPT = "submit_custom_list_query_task.sh"


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
