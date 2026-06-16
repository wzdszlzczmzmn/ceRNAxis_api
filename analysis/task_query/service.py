import uuid as uuid_lib

from .registry import (
    find_task_by_uuid,
    TaskNotFoundError,
    MultipleTaskMatchedError,
)
from .status_sync import (
    sync_task_status,
    TaskStatusSyncError,
)


class InvalidTaskUUIDError(ValueError):
    pass


def validate_task_uuid(task_uuid: str) -> None:
    try:
        uuid_lib.UUID(str(task_uuid))
    except ValueError as e:
        raise InvalidTaskUUIDError("Invalid Task UUID format.") from e


def query_task_by_uuid(task_uuid: str) -> dict:
    validate_task_uuid(task_uuid)

    task, config = find_task_by_uuid(task_uuid)

    sync_result = sync_task_status(task)

    return {
        "task_type": config.task_type,
        "data": config.formatter(
            task,
            position=sync_result.position,
        ),
    }
