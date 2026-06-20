import uuid as uuid_lib

from .registry import (
    get_query_config_by_task_type,
)
from .status_sync import (
    sync_task_status,
)
from ..task_common.registry import find_task_by_uuid


class InvalidTaskUUIDError(ValueError):
    pass


def validate_task_uuid(task_uuid: str) -> None:
    try:
        uuid_lib.UUID(str(task_uuid))
    except ValueError as e:
        raise InvalidTaskUUIDError("Invalid Task UUID format.") from e


def query_task_by_uuid(task_uuid: str) -> dict:
    validate_task_uuid(task_uuid)

    task, task_model_config = find_task_by_uuid(task_uuid)

    sync_result = sync_task_status(task)

    query_config = get_query_config_by_task_type(
        task_model_config.task_type
    )

    return {
        "task_type": task_model_config.task_type,
        "data": query_config.formatter(
            task,
            position=sync_result.position,
        ),
    }
