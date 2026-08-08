from collections import OrderedDict
from pathlib import Path
import csv
import anndata as ad
import pandas as pd

from analysis.utils.hybrid_reference_task_utils import (
    SCSTHybridReferenceTaskInputError,
    SCSTHybridReferenceTaskPathError,
    get_scst_expected_id_column,
    get_scst_hybrid_reference_input_file_path,
    validate_scst_meta_csv_schema,
)


class WorkflowVizInfoInputError(ValueError):
    """
    Raised when workflow visualization information cannot be
    extracted from task input files.
    """


class WorkflowVizInfoPathError(ValueError):
    """
    Raised when a workflow visualization information file path
    is invalid.
    """


def get_scst_hybrid_reference_meta_file_path(task) -> Path:
    """
    Return the resolved meta.csv path for an SC/ST Hybrid Reference task.

    Expected location:
        {task_workspace}/input/meta.csv

    The actual filename is resolved through the existing
    SC/ST Hybrid Reference input-file path utility.
    """
    try:
        file_path = get_scst_hybrid_reference_input_file_path(
            task=task,
            field_name="meta_file",
        )
    except SCSTHybridReferenceTaskPathError as exc:
        raise WorkflowVizInfoPathError(
            str(exc)
        ) from exc

    return file_path


def validate_scst_hybrid_reference_meta_file(task) -> Path:
    """
    Validate that the SC/ST task metadata file exists and is a file.
    """
    file_path = get_scst_hybrid_reference_meta_file_path(task)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            "SC/ST Hybrid Reference metadata file not found: "
            f"{file_path.name}"
        )

    return file_path


def normalize_group_value(value) -> str:
    """
    Normalize a value read from the metadata group column.

    Current behavior:
    - None becomes an empty string.
    - Other values are converted to string.
    - Leading and trailing whitespace is removed.
    """
    if value is None:
        return ""

    return str(value).strip()


def read_scst_group_value_counts(
    *,
    file_path: Path,
    group_col: str,
) -> tuple[OrderedDict[str, int], int]:
    """
    Read group values from an SC/ST metadata CSV.

    Returns:
        (
            OrderedDict({
                group_value: row_count,
                ...
            }),
            total_data_row_count,
        )

    Group ordering follows the first occurrence in meta.csv.
    """
    group_col = str(group_col or "").strip()

    if not group_col:
        raise WorkflowVizInfoInputError(
            "SC/ST Hybrid Reference task does not define group_col."
        )

    group_counts: OrderedDict[str, int] = OrderedDict()
    row_count = 0

    try:
        with file_path.open(
            "r",
            newline="",
            encoding="utf-8-sig",
        ) as file_obj:
            reader = csv.DictReader(file_obj)

            if not reader.fieldnames:
                raise WorkflowVizInfoInputError(
                    "Metadata file is empty or missing a header."
                )

            normalized_fieldnames = [
                str(column).strip()
                for column in reader.fieldnames
            ]

            reader.fieldnames = normalized_fieldnames

            if group_col not in normalized_fieldnames:
                raise WorkflowVizInfoInputError(
                    "Metadata file is missing group column: "
                    f"'{group_col}'."
                )

            for row_number, row in enumerate(
                reader,
                start=2,
            ):
                row_count += 1

                group_value = normalize_group_value(
                    row.get(group_col)
                )

                if not group_value:
                    raise WorkflowVizInfoInputError(
                        "Metadata file contains an empty group value "
                        f"in column '{group_col}' at row {row_number}."
                    )

                if group_value not in group_counts:
                    group_counts[group_value] = 0

                group_counts[group_value] += 1

    except UnicodeDecodeError as exc:
        raise WorkflowVizInfoInputError(
            "Metadata file must be UTF-8 encoded: "
            f"{file_path.name}."
        ) from exc

    except csv.Error as exc:
        raise WorkflowVizInfoInputError(
            f"Invalid metadata CSV file: {file_path.name}. "
            f"{str(exc)}"
        ) from exc

    if row_count == 0:
        raise WorkflowVizInfoInputError(
            "Metadata file has no data rows."
        )

    if not group_counts:
        raise WorkflowVizInfoInputError(
            f"Metadata group column '{group_col}' has no values."
        )

    return group_counts, row_count


def get_scst_hybrid_reference_group_info(task) -> dict:
    """
    Extract group information from an SC/ST Hybrid Reference task's
    metadata file.

    The task's group_col determines which metadata column is read.

    Returned structure:
        {
            "meta_file": "meta.csv",
            "id_column": "cell_id",
            "group_col": "Celltype",
            "group_count": 3,
            "sample_count": 100,
            "group_values": [
                "T cell",
                "B cell",
                "Myeloid",
            ],
            "group_options": [
                {
                    "value": "T cell",
                    "label": "T cell",
                    "count": 40,
                },
                ...
            ],
        }
    """
    group_col = str(
        getattr(task, "group_col", "") or ""
    ).strip()

    if not group_col:
        raise WorkflowVizInfoInputError(
            "SC/ST Hybrid Reference task does not define group_col."
        )

    data_type = str(
        getattr(task, "data_type", "") or ""
    ).strip()

    try:
        expected_id_column = get_scst_expected_id_column(
            data_type
        )
    except SCSTHybridReferenceTaskInputError as exc:
        raise WorkflowVizInfoInputError(
            str(exc)
        ) from exc

    meta_file_path = validate_scst_hybrid_reference_meta_file(
        task
    )

    # Reuse the existing metadata schema validation:
    # - SC requires cell_id as the first column.
    # - ST requires spot_id as the first column.
    # - group_col must exist.
    # - group_col cannot equal the identifier column.
    try:
        validate_scst_meta_csv_schema(
            file_path=meta_file_path,
            expected_id_column=expected_id_column,
            group_col=group_col,
        )
    except SCSTHybridReferenceTaskInputError as exc:
        raise WorkflowVizInfoInputError(
            str(exc)
        ) from exc

    group_counts, sample_count = (
        read_scst_group_value_counts(
            file_path=meta_file_path,
            group_col=group_col,
        )
    )

    group_values = list(group_counts.keys())

    group_options = [
        {
            "value": group_value,
            "label": group_value,
            "count": count,
        }
        for group_value, count in group_counts.items()
    ]

    return {
        "meta_file": meta_file_path.name,
        "id_column": expected_id_column,
        "group_col": group_col,
        "group_count": len(group_values),
        "sample_count": sample_count,
        "group_values": group_values,
        "group_options": group_options,
    }


def validate_scst_hybrid_reference_group_value(
    *,
    task,
    group_value: str,
) -> str:
    group_value = str(
        group_value or ""
    ).strip()

    if not group_value:
        raise WorkflowVizInfoInputError(
            "Missing groupValue."
        )

    group_info = (
        get_scst_hybrid_reference_group_info(
            task
        )
    )

    valid_group_values = (
        group_info.get("group_values") or []
    )

    if group_value not in valid_group_values:
        raise WorkflowVizInfoInputError(
            "Invalid groupValue. "
            "Allowed values are: "
            f"{', '.join(valid_group_values)}."
        )

    return group_value


def get_scst_hybrid_reference_expression_file_path(task) -> Path:
    """
    Return the resolved expression.h5ad path for an SC/ST Hybrid Reference task.

    SC/ST is now H5AD-only, so resolve the uploaded input directly from
    the current SC/ST input-file mapping using field_name="exp_file".
    """
    try:
        return get_scst_hybrid_reference_input_file_path(
            task=task,
            field_name="exp_file",
        )
    except SCSTHybridReferenceTaskPathError as exc:
        raise WorkflowVizInfoPathError(
            str(exc)
        ) from exc


def validate_scst_hybrid_reference_expression_file(task) -> Path:
    """
    Validate that expression.h5ad exists and is a file.
    """
    file_path = get_scst_hybrid_reference_expression_file_path(task)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            "SC/ST Hybrid Reference expression file not found: "
            f"{file_path.name}"
        )

    return file_path


def read_scst_h5ad_group_value_counts(
    *,
    file_path: Path,
    group_col: str,
) -> tuple[OrderedDict[str, int], int]:
    """
    Read group values from expression.h5ad -> adata.obs[group_col].

    Returns:
        (
            OrderedDict({
                group_value: cell_count,
                ...
            }),
            total_cell_count,
        )

    Group ordering follows the first occurrence in adata.obs.
    """
    group_col = str(group_col or "").strip()

    if not group_col:
        raise WorkflowVizInfoInputError(
            "SC/ST Hybrid Reference task does not define group_col."
        )

    try:
        adata = ad.read_h5ad(
            file_path,
            backed="r",
        )
    except (OSError, ValueError) as exc:
        raise WorkflowVizInfoInputError(
            f"Unable to read H5AD file: {file_path.name}. "
            f"{str(exc)}"
        ) from exc

    try:
        if group_col not in adata.obs.columns:
            raise WorkflowVizInfoInputError(
                "H5AD obs is missing group column: "
                f"'{group_col}'."
            )

        group_series = adata.obs[group_col]
        cell_count = adata.n_obs

        if cell_count == 0:
            raise WorkflowVizInfoInputError(
                "H5AD file has no observations/cells."
            )

        group_counts: OrderedDict[str, int] = OrderedDict()

        for obs_name, value in group_series.items():
            if pd.isna(value):
                group_value = ""
            else:
                group_value = normalize_group_value(value)

            if not group_value:
                raise WorkflowVizInfoInputError(
                    "H5AD obs contains an empty group value "
                    f"in column '{group_col}' for observation "
                    f"'{obs_name}'."
                )

            if group_value not in group_counts:
                group_counts[group_value] = 0

            group_counts[group_value] += 1

        if not group_counts:
            raise WorkflowVizInfoInputError(
                f"H5AD obs group column '{group_col}' has no values."
            )

        return group_counts, cell_count

    finally:
        if getattr(adata, "isbacked", False):
            adata.file.close()


def get_scst_hybrid_reference_h5ad_group_info(task) -> dict:
    """
    Extract group information from expression.h5ad -> adata.obs.

    The task's group_col determines which obs column is read.

    Returned structure:
        {
            "expression_file": "expression.h5ad",
            "id_column": "cell_id",
            "group_col": "Celltype (malignancy)",
            "group_count": 3,
            "sample_count": 33043,
            "group_values": [
                "Malignant",
                "Non-malignant",
                ...
            ],
            "group_options": [
                {
                    "value": "Malignant",
                    "label": "Malignant",
                    "count": 10000,
                },
                ...
            ],
        }
    """
    group_col = str(
        getattr(task, "group_col", "") or ""
    ).strip()

    if not group_col:
        raise WorkflowVizInfoInputError(
            "SC/ST Hybrid Reference task does not define group_col."
        )

    data_type = str(
        getattr(task, "data_type", "") or ""
    ).strip()

    try:
        expected_id_column = get_scst_expected_id_column(
            data_type
        )
    except SCSTHybridReferenceTaskInputError as exc:
        raise WorkflowVizInfoInputError(
            str(exc)
        ) from exc

    expression_file_path = (
        validate_scst_hybrid_reference_expression_file(task)
    )

    group_counts, sample_count = (
        read_scst_h5ad_group_value_counts(
            file_path=expression_file_path,
            group_col=group_col,
        )
    )

    group_values = list(group_counts.keys())

    group_options = [
        {
            "value": group_value,
            "label": group_value,
            "count": count,
        }
        for group_value, count in group_counts.items()
    ]

    return {
        "expression_file": expression_file_path.name,
        "id_column": expected_id_column,
        "group_col": group_col,
        "group_count": len(group_values),
        "sample_count": sample_count,
        "group_values": group_values,
        "group_options": group_options,
    }


def validate_scst_hybrid_reference_h5ad_group_value(
    *,
    task,
    group_value: str,
) -> str:
    """
    Validate groupValue against expression.h5ad -> adata.obs[group_col].
    """
    group_value = str(
        group_value or ""
    ).strip()

    if not group_value:
        raise WorkflowVizInfoInputError(
            "Missing groupValue."
        )

    group_info = (
        get_scst_hybrid_reference_h5ad_group_info(
            task
        )
    )

    valid_group_values = (
        group_info.get("group_values") or []
    )

    if group_value not in valid_group_values:
        raise WorkflowVizInfoInputError(
            "Invalid groupValue. "
            "Allowed values are: "
            f"{', '.join(valid_group_values)}."
        )

    return group_value
