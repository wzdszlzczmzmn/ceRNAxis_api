from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from django.utils import timezone

from analysis.slurm_squeue import squeue_by_job_name


STATUS_FILENAME = "status.txt"


@dataclass
class TaskStatusResult:
    position: int = 0


class TaskStatusSyncError(Exception):
    pass


def sync_task_status(task) -> TaskStatusResult:
    """
    Sync task status using Slurm queue and status.txt.

    Returns:
        TaskStatusResult(position=...)
    """

    if task.status in [
        task.Status.Success,
        task.Status.Failed,
    ]:
        return TaskStatusResult(position=0)

    task_status = squeue_by_job_name(str(task.uuid))

    if task_status == "empty":
        return sync_task_status_from_status_file(task)

    if task_status == "R":
        task.status = task.Status.Running
        task.save(update_fields=["status"])
        return TaskStatusResult(position=0)

    if isinstance(task_status, str) and task_status.startswith("PD"):
        task.status = task.Status.Pending
        task.save(update_fields=["status"])

        return TaskStatusResult(
            position=parse_pending_position(task_status)
        )

    raise TaskStatusSyncError(
        f"Unknown task status from Slurm: {task_status}"
    )


def sync_task_status_from_status_file(task) -> TaskStatusResult:
    status_file_path = (
        Path(task.get_output_dir_absolute_path())
        / STATUS_FILENAME
    )

    if not status_file_path.exists() or not status_file_path.is_file():
        raise TaskStatusSyncError(
            "Status file not found. The task may be in an inconsistent state."
        )

    finish_time, raw_status = read_status_file(status_file_path)

    if raw_status == "success":
        task.status = task.Status.Success
    else:
        task.status = task.Status.Failed

    task.finish_time = finish_time
    task.save(update_fields=["status", "finish_time"])

    return TaskStatusResult(position=0)


def read_status_file(status_file_path: Path):
    with status_file_path.open("r", encoding="utf-8") as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]

    if len(lines) < 2:
        raise TaskStatusSyncError("Invalid status.txt format.")

    finish_time_str = lines[0]
    raw_status = lines[1].lower()

    try:
        naive_finish_time = datetime.strptime(
            finish_time_str,
            "%Y-%m-%d %H:%M:%S",
        )
    except ValueError as e:
        raise TaskStatusSyncError(
            f"Invalid finish_time in status.txt: {finish_time_str}"
        ) from e

    finish_time = timezone.make_aware(
        naive_finish_time,
        timezone.get_current_timezone(),
    )

    return finish_time, raw_status


def parse_pending_position(task_status: str) -> int:
    """
    Example:
        "PD 3" -> 3
        "PD" -> 0
    """

    parts = task_status.split()

    if len(parts) < 2:
        return 0

    try:
        return int(parts[1])
    except ValueError:
        return 0
