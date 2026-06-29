import uuid as uuid_lib


class WorkflowNetworkViewError(Exception):
    status_code = 400

    def __init__(self, msg: str):
        self.msg = msg
        super().__init__(msg)


class MissingTaskUUIDError(WorkflowNetworkViewError):
    status_code = 400


class InvalidTaskUUIDError(WorkflowNetworkViewError):
    status_code = 400


class TaskNotFoundError(WorkflowNetworkViewError):
    status_code = 404


class TaskNotSuccessError(WorkflowNetworkViewError):
    status_code = 400


def get_required_task_uuid(request) -> str:
    task_uuid = str(request.query_params.get("taskUUID", "")).strip()

    if not task_uuid:
        raise MissingTaskUUIDError("Missing required parameter: taskUUID.")

    try:
        uuid_lib.UUID(task_uuid)
    except ValueError as e:
        raise InvalidTaskUUIDError("Invalid taskUUID format.") from e

    return task_uuid


def get_task_or_raise(model_class, task_uuid: str, task_label: str):
    try:
        return model_class.objects.get(uuid=task_uuid)
    except model_class.DoesNotExist as e:
        raise TaskNotFoundError(
            f"{task_label} with UUID {task_uuid} not found."
        ) from e


def require_success_task(task, task_label: str) -> None:
    if task.status != task.Status.Success:
        raise TaskNotSuccessError(
            f"{task_label} is not completed successfully. "
            f"Current status: {task.get_status_display()}."
        )


def attach_task_metadata(result: dict, task, task_type: str, extra: dict | None = None) -> dict:
    result["task_type"] = task_type
    result["task_uuid"] = str(task.uuid)
    result["task_name"] = task.task_name

    if hasattr(task, "map_info"):
        result["map_info"] = task.map_info

    if extra:
        result.update(extra)

    return result
