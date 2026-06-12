from pathlib import Path

from django.conf import settings


class SlurmPathError(ValueError):
    pass


PIPELINE_STDOUT_FILENAME = "Pipeline.out"
PIPELINE_STDERR_FILENAME = "Pipeline.err"


def validate_safe_script_name(script_name: str) -> None:
    if not script_name:
        raise SlurmPathError("Missing script_name.")

    if "/" in script_name or "\\" in script_name or ".." in script_name:
        raise SlurmPathError("Invalid script_name.")


def get_slurm_script_home() -> Path:
    return Path(settings.SLURM_SCRIPT_HOME).resolve()


def get_slurm_script_path(script_name: str) -> Path:
    validate_safe_script_name(script_name)

    script_home = get_slurm_script_home()
    script_path = (script_home / script_name).resolve()

    if not str(script_path).startswith(str(script_home)):
        raise SlurmPathError("Invalid SLURM script path.")

    return script_path


def validate_slurm_script(script_name: str) -> Path:
    script_path = get_slurm_script_path(script_name)

    if not script_path.exists() or not script_path.is_file():
        raise FileNotFoundError(f"SLURM script not found: {script_name}")

    return script_path


def get_slurm_job_name(task_uuid) -> str:
    return str(task_uuid).replace("-", "_")


def get_pipeline_stdout_file_path(output_dir) -> Path:
    return Path(output_dir).resolve() / PIPELINE_STDOUT_FILENAME


def get_pipeline_stderr_file_path(output_dir) -> Path:
    return Path(output_dir).resolve() / PIPELINE_STDERR_FILENAME
