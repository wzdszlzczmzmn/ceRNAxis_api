from dataclasses import dataclass
from typing import Callable

from .formatters import (
    format_custom_list_query_task,
    format_paired_cohort_task, format_hybrid_reference_task, format_scst_hybrid_reference_task,
)


@dataclass(frozen=True)
class TaskQueryConfig:
    task_type: str
    formatter: Callable


TASK_QUERY_REGISTRY = [
    TaskQueryConfig(
        task_type="CustomListQueryTask",
        formatter=format_custom_list_query_task,
    ),
    TaskQueryConfig(
        task_type="PairedCohortTask",
        formatter=format_paired_cohort_task,
    ),
    TaskQueryConfig(
        task_type="HybridReferenceTask",
        formatter=format_hybrid_reference_task,
    ),
    TaskQueryConfig(
        task_type="SCSTHybridReferenceTask",
        formatter=format_scst_hybrid_reference_task,
    ),
]


class TaskQueryConfigNotFoundError(Exception):
    pass


def get_query_config_by_task_type(task_type: str) -> TaskQueryConfig:
    for config in TASK_QUERY_REGISTRY:
        if config.task_type == task_type:
            return config

    raise TaskQueryConfigNotFoundError(
        f"Task query config not found for task_type: {task_type}."
    )
