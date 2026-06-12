import subprocess

from analysis.utils.slurm_path_utils import get_slurm_script_path


def normalize_slurm_job_name(task_uuid: str) -> str:
    return str(task_uuid).replace("-", "_")


def squeue_by_job_name(task_uuid: str) -> str | bool:
    script_path = get_slurm_script_path("task_query.sh")

    result = subprocess.run(
        [
            str(script_path),
            normalize_slurm_job_name(task_uuid),
        ],
        capture_output=True,
        text=True,
    )

    if result.returncode == 0:
        return result.stdout.strip()

    return False
