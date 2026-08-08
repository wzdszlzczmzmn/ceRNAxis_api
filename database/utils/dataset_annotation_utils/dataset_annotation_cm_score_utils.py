from pathlib import Path

import pandas as pd

from analysis.utils.workflow_detail_utils.workflow_cm_score_utils import (
    build_cm_score_options_response_common,
    build_cm_score_result_response_common,
    get_available_cm_score_items_from_dir,
    get_cm_score_file_path_from_dir,
    read_cm_score_file_by_path,
)


def build_dataset_annotation_cm_score_options_response(
    *,
    base_response: dict,
    cm_results_dir: Path,
) -> dict:
    """
    Build Dataset Annotation CM-score item options.

    Missing/non-directory CM-results paths intentionally produce an
    empty option list, matching Workflow CM Score Options behavior.
    """
    cm_results_dir = Path(
        cm_results_dir
    ).resolve()

    items = get_available_cm_score_items_from_dir(
        cm_results_dir=cm_results_dir,
    )

    return build_cm_score_options_response_common(
        base_response={
            **base_response,
            "cm_results_dir": cm_results_dir.name,
        },
        items=items,
    )


def read_dataset_annotation_cm_score_file(
    *,
    cm_results_dir: Path,
    item_value: str,
) -> tuple[Path, pd.DataFrame]:
    """
    Resolve, read and validate one Dataset Annotation CM-score CSV.
    """
    file_path = get_cm_score_file_path_from_dir(
        cm_results_dir=cm_results_dir,
        item_value=item_value,
    )

    return read_cm_score_file_by_path(
        file_path=file_path,
    )


def build_dataset_annotation_cm_score_result_response(
    *,
    base_response: dict,
    item_value: str,
    file_path: Path,
    dataframe: pd.DataFrame,
) -> dict:
    """
    Build the same CM-score data contract used by Workflow.
    """
    return build_cm_score_result_response_common(
        base_response=base_response,
        item_value=item_value,
        file_path=file_path,
        dataframe=dataframe,
    )
