from dataclasses import dataclass
from typing import Callable, Type

from django.db import models

from analysis.models import CustomListQueryTask

from .formatters import (
    format_custom_list_query_task,
)


@dataclass(frozen=True)
class TaskTypeConfig:
    task_type: str
    model: Type[models.Model]
    formatter: Callable


TASK_REGISTRY = [
    TaskTypeConfig(
        task_type="CustomListQueryTask",
        model=CustomListQueryTask,
        formatter=format_custom_list_query_task,
    ),
]


class TaskNotFoundError(Exception):
    pass


class MultipleTaskMatchedError(Exception):
    pass


def find_task_by_uuid(task_uuid: str):
    matched = []

    for config in TASK_REGISTRY:
        try:
            task = config.model.objects.get(uuid=task_uuid)
            matched.append((task, config))
        except config.model.DoesNotExist:
            continue

    if not matched:
        raise TaskNotFoundError(
            f"Task with UUID {task_uuid} not found."
        )

    if len(matched) > 1:
        raise MultipleTaskMatchedError(
            f"Multiple task records found for UUID {task_uuid}."
        )

    return matched[0]
