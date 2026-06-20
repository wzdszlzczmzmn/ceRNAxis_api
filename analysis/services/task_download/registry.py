from dataclasses import dataclass
from typing import Callable

from django.db import models

from analysis.models import CustomListQueryTask, PairedCohortTask

from .resolvers import (
    build_custom_list_query_archive_name,
    build_paired_cohort_archive_name,
    resolve_custom_list_query_result_files,
    resolve_paired_cohort_result_files,
)


@dataclass(frozen=True)
class TaskDownloadConfig:
    task_type: str
    model: type[models.Model]
    result_file_resolver: Callable
    archive_name_builder: Callable


TASK_DOWNLOAD_REGISTRY = [
    TaskDownloadConfig(
        task_type="CustomListQueryTask",
        model=CustomListQueryTask,
        result_file_resolver=resolve_custom_list_query_result_files,
        archive_name_builder=build_custom_list_query_archive_name,
    ),
    TaskDownloadConfig(
        task_type="PairedCohortTask",
        model=PairedCohortTask,
        result_file_resolver=resolve_paired_cohort_result_files,
        archive_name_builder=build_paired_cohort_archive_name,
    ),
]


class TaskDownloadConfigNotFoundError(Exception):
    pass


def get_download_config_by_task_type(task_type: str) -> TaskDownloadConfig:
    for config in TASK_DOWNLOAD_REGISTRY:
        if config.task_type == task_type:
            return config

    raise TaskDownloadConfigNotFoundError(
        f"Task download config not found for task_type: {task_type}."
    )
