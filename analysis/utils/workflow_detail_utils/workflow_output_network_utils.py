from pathlib import Path

from analysis.utils.paired_cohort_task_utils import (
    validate_safe_name,
    PairedCohortTaskPathError,
)


WORKFLOW_CERNA_AXIS_FILENAME_SUFFIX = "_ceRNA_axis.csv"
WORKFLOW_IMMUNE_AXIS_FILENAME_SUFFIX = "_map_immune_axis.csv"


SCST_CERNA_AXIS_FILENAME_TEMPLATE = (
    "{task_name}_ceRNA_{group_value}_axis.csv"
)


SCST_IMMUNE_AXIS_FILENAME_TEMPLATE = (
    "{task_name}_map_immune_axis_{group_value}.csv"
)


class WorkflowOutputNetworkPathError(ValueError):
    pass


def get_workflow_output_dir(task) -> Path:
    return Path(task.get_output_dir_absolute_path()).resolve()


def get_workflow_output_file_path(
    *,
    task,
    filename_suffix: str,
    error_message: str,
) -> Path:
    task_name = str(task.task_name).strip()

    validate_safe_name(task_name, "task_name")

    output_dir = get_workflow_output_dir(task)

    file_path = (
        output_dir / f"{task_name}{filename_suffix}"
    ).resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise WorkflowOutputNetworkPathError(error_message)

    return file_path


def get_workflow_cerna_axis_file_path(task) -> Path:
    return get_workflow_output_file_path(
        task=task,
        filename_suffix=WORKFLOW_CERNA_AXIS_FILENAME_SUFFIX,
        error_message="Invalid workflow ceRNA axis file path.",
    )


def get_workflow_immune_axis_file_path(task) -> Path:
    return get_workflow_output_file_path(
        task=task,
        filename_suffix=WORKFLOW_IMMUNE_AXIS_FILENAME_SUFFIX,
        error_message="Invalid workflow immune axis file path.",
    )


def validate_scst_group_value(
    group_value: str,
) -> str:

    group_value = str(
        group_value or ""
    ).strip()

    if not group_value:
        raise WorkflowOutputNetworkPathError(
            "Missing group value."
        )

    if (
        "/" in group_value
        or "\\" in group_value
        or ".." in group_value
    ):
        raise WorkflowOutputNetworkPathError(
            "Invalid group value."
        )

    return group_value


def get_scst_workflow_cerna_axis_file_path(
    task,
    group_value,
) -> Path:

    group_value = validate_scst_group_value(
        group_value
    )

    output_dir = get_workflow_output_dir(task)

    filename = (
        SCST_CERNA_AXIS_FILENAME_TEMPLATE.format(
            task_name=task.task_name,
            group_value=group_value,
        )
    )

    file_path = (
        output_dir / filename
    ).resolve()

    if not str(file_path).startswith(
        str(output_dir)
    ):
        raise WorkflowOutputNetworkPathError(
            "Invalid SC/ST ceRNA axis path."
        )

    return file_path


def get_scst_workflow_immune_axis_file_path(
    task,
    group_value,
) -> Path:

    group_value = validate_scst_group_value(
        group_value
    )

    output_dir = get_workflow_output_dir(task)

    filename = (
        SCST_IMMUNE_AXIS_FILENAME_TEMPLATE.format(
            task_name=task.task_name,
            group_value=group_value,
        )
    )

    file_path = (
        output_dir / filename
    ).resolve()

    if not str(file_path).startswith(
        str(output_dir)
    ):
        raise WorkflowOutputNetworkPathError(
            "Invalid SC/ST immune axis path."
        )

    return file_path
