from dataclasses import dataclass
from typing import Type

from django.db import models

from analysis.models import CustomListQueryTask, PairedCohortTask, HybridReferenceTask, SCSTHybridReferenceTask


@dataclass(frozen=True)
class TaskModelConfig:
    task_type: str
    model: Type[models.Model]


TASK_MODEL_REGISTRY = [
    TaskModelConfig(
        task_type="CustomListQueryTask",
        model=CustomListQueryTask,
    ),
    TaskModelConfig(
        task_type="PairedCohortTask",
        model=PairedCohortTask,
    ),
    TaskModelConfig(
        task_type="HybridReferenceTask",
        model=HybridReferenceTask,
    ),
    TaskModelConfig(
        task_type="SCSTHybridReferenceTask",
        model=SCSTHybridReferenceTask,
    ),
]


class TaskNotFoundError(Exception):
    pass


class MultipleTaskMatchedError(Exception):
    pass


def find_task_by_uuid(task_uuid: str):
    matched = []

    for config in TASK_MODEL_REGISTRY:
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
