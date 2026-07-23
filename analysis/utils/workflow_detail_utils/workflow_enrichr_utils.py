import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.utils.custom_list_query_task_utils import (
    get_task_output_dir,
    validate_task_name_for_filename,
)


CUSTOM_LIST_QUERY_ENRICHR_DIRECTIONS = {
    "up",
    "down",
}

CUSTOM_LIST_QUERY_ENRICHR_REQUIRED_COLUMNS = {
    "Gene_set",
    "Term",
    "Overlap",
    "P-value",
    "Adjusted P-value",
    "Odds Ratio",
    "Combined Score",
    "Genes",
}

OVERLAP_PATTERN = re.compile(
    r"^\s*(\d+)\s*/\s*(\d+)\s*$"
)


class CustomListQueryEnrichrInputError(ValueError):
    pass


class CustomListQueryEnrichrPathError(ValueError):
    pass


def normalize_enrichr_direction(direction: str) -> str:
    normalized_direction = str(
        direction or ""
    ).strip().lower()

    if normalized_direction not in (
        CUSTOM_LIST_QUERY_ENRICHR_DIRECTIONS
    ):
        raise CustomListQueryEnrichrInputError(
            "Invalid direction. Allowed values are: up, down."
        )

    return normalized_direction


def get_custom_list_query_enrichr_filename(
    task_name: str,
    direction: str,
) -> str:
    normalized_task_name = str(
        task_name or ""
    ).strip()

    validate_task_name_for_filename(
        normalized_task_name
    )

    normalized_direction = normalize_enrichr_direction(
        direction
    )

    return (
        f"{normalized_task_name}"
        f"_mRNA_{normalized_direction}_enrichr.csv"
    )


def get_custom_list_query_enrichr_file_path(
    task,
    direction: str,
) -> Path:
    output_dir = get_task_output_dir(task)

    filename = get_custom_list_query_enrichr_filename(
        task_name=task.task_name,
        direction=direction,
    )

    file_path = (
        output_dir / filename
    ).resolve()

    try:
        file_path.relative_to(output_dir)
    except ValueError as exc:
        raise CustomListQueryEnrichrPathError(
            "Invalid Enrichr result file path."
        ) from exc

    return file_path


def validate_custom_list_query_enrichr_file(
    task,
    direction: str,
) -> Path:
    file_path = get_custom_list_query_enrichr_file_path(
        task=task,
        direction=direction,
    )

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Enrichr result file not found: {file_path.name}"
        )

    return file_path


def read_custom_list_query_enrichr_file(
    task,
    direction: str,
) -> tuple[Path, pd.DataFrame]:
    file_path = validate_custom_list_query_enrichr_file(
        task=task,
        direction=direction,
    )

    try:
        dataframe = pd.read_csv(file_path)
    except Exception as exc:
        raise CustomListQueryEnrichrInputError(
            f"Failed to read Enrichr result file: {exc}"
        ) from exc

    missing_columns = (
        CUSTOM_LIST_QUERY_ENRICHR_REQUIRED_COLUMNS
        - set(dataframe.columns)
    )

    if missing_columns:
        raise CustomListQueryEnrichrInputError(
            "Enrichr result file is missing required columns: "
            f"{', '.join(sorted(missing_columns))}."
        )

    return file_path, dataframe


def parse_overlap_value(
    value,
) -> tuple[int | None, int | None]:
    if value is None:
        return None, None

    matched = OVERLAP_PATTERN.fullmatch(
        str(value)
    )

    if matched is None:
        return None, None

    return (
        int(matched.group(1)),
        int(matched.group(2)),
    )


def parse_gene_list(value) -> list[str]:
    if value is None:
        return []

    try:
        if pd.isna(value):
            return []
    except (TypeError, ValueError):
        pass

    return [
        gene.strip()
        for gene in str(value).split(";")
        if gene.strip()
    ]


def to_optional_float(value):
    if value is None:
        return None

    try:
        numeric_value = float(value)
    except (TypeError, ValueError):
        return None

    if not math.isfinite(numeric_value):
        return None

    return numeric_value


def normalize_dataframe_value(value):
    if value is None:
        return None

    if isinstance(value, np.generic):
        value = value.item()

    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass

    if isinstance(value, float):
        if not math.isfinite(value):
            return None

    return value


def serialize_custom_list_query_enrichr_dataframe(
    dataframe: pd.DataFrame,
    direction: str,
) -> list[dict]:
    normalized_direction = normalize_enrichr_direction(
        direction
    )

    results = []

    for row_number, (_, row) in enumerate(
        dataframe.iterrows(),
        start=1,
    ):
        overlap_count, overlap_total = parse_overlap_value(
            row.get("Overlap")
        )

        genes = parse_gene_list(
            row.get("Genes")
        )

        results.append(
            {
                "id": row_number,
                "direction": normalized_direction,

                "gene_set": normalize_dataframe_value(
                    row.get("Gene_set")
                ),
                "term": normalize_dataframe_value(
                    row.get("Term")
                ),

                "overlap": normalize_dataframe_value(
                    row.get("Overlap")
                ),
                "overlap_count": overlap_count,
                "overlap_total": overlap_total,

                "p_value": to_optional_float(
                    row.get("P-value")
                ),
                "adjusted_p_value": to_optional_float(
                    row.get("Adjusted P-value")
                ),
                "old_p_value": to_optional_float(
                    row.get("Old P-value")
                ),
                "old_adjusted_p_value": to_optional_float(
                    row.get("Old Adjusted P-value")
                ),

                "odds_ratio": to_optional_float(
                    row.get("Odds Ratio")
                ),
                "combined_score": to_optional_float(
                    row.get("Combined Score")
                ),

                "genes": genes,
                "gene_count": len(genes),
            }
        )

    return results


def build_custom_list_query_enrichr_response(
    *,
    task,
    direction: str,
    file_path: Path,
    dataframe: pd.DataFrame,
) -> dict:
    normalized_direction = normalize_enrichr_direction(
        direction
    )

    results = serialize_custom_list_query_enrichr_dataframe(
        dataframe=dataframe,
        direction=normalized_direction,
    )

    return {
        "uuid": str(task.uuid),
        "task_type": "CustomListQueryTask",
        "task_name": task.task_name,

        "cancer_type": getattr(
            task,
            "cancer_type",
            "",
        ) or "",

        "has_mrna_direction": bool(
            getattr(
                task,
                "has_mrna_direction",
                False,
            )
        ),

        "direction": normalized_direction,
        "enrichr_file": file_path.name,

        "plot": {
            "chart_type": "bar",
            "orientation": "horizontal",
            "x_field": "combined_score",
            "y_field": "term",
        },

        "summary": {
            "raw_count": len(dataframe),
            "returned_count": len(results),
        },

        "results": results,
    }
