from collections import OrderedDict
from pathlib import Path
import csv

from analysis.utils.workflow_detail_utils.workflow_deg_volcano_utils import (
    SCST_WORKFLOW_DEG_SCOPE_ALL,
    SCST_WORKFLOW_DEG_SCOPES,
)
from analysis.utils.workflow_detail_utils.workflow_log2fc_background_utils import (
    SCST_HYBRID_REFERENCE_VALID_BACKGROUND_TYPES,
    WorkflowLog2FCBackgroundInputError,
    get_workflow_available_background_types,
    read_log2fc_background_file_by_path,
)

from database.utils.dataset_annotation_utils.path_utils import (
    DatasetAnnotationPathError,
)
from database.utils.dataset_annotation_utils.scst_path_utils import (
    build_scst_group_value_result_paths,
    get_scst_dataset_group_by_fields,
    get_scst_dataset_id_column,
    get_scst_dataset_meta_file_path,
    get_scst_group_by_meta_column,
    is_scst_existing_file,
    is_scst_non_empty_file_directory,
    resolve_scst_dataset_group_annotation_dir,
    should_skip_scst_group_value_for_results,
    validate_scst_data_type,
)


class SCSTDatasetAnnotationMetadataError(
    DatasetAnnotationPathError
):
    """
    Raised when configured SC/ST dataset metadata
    is missing or inconsistent.
    """


def normalize_scst_group_value(
    value,
) -> str:
    if value is None:
        return ""

    return str(value).strip()


def read_scst_dataset_group_value_counts(
    *,
    meta_file: Path,
    group_by_fields: list[str],
) -> tuple[
    dict[str, OrderedDict[str, int]],
    int,
]:
    """
    Read all configured group-by columns in one pass.

    For each configured group-by value:
    - map it to the actual metadata CSV column
    - preserve first-occurrence order
    - collect unique values
    - collect row counts

    The returned dictionary remains keyed by the configured
    group-by value, not by the physical metadata column name.

    Returns:
        (
            {
                group_by: OrderedDict({
                    group_value: count,
                    ...
                }),
                ...
            },
            total_row_count,
        )
    """
    meta_file = Path(
        meta_file
    ).resolve()

    if (
        not meta_file.exists()
        or not meta_file.is_file()
    ):
        raise FileNotFoundError(
            "SC/ST dataset metadata file not found: "
            f"{meta_file.name}"
        )

    normalized_group_by_fields = []

    for raw_group_by in group_by_fields:
        group_by = str(
            raw_group_by or ""
        ).strip()

        if (
            group_by
            and group_by
            not in normalized_group_by_fields
        ):
            normalized_group_by_fields.append(
                group_by
            )

    meta_column_by_group_by = {
        group_by: get_scst_group_by_meta_column(
            group_by
        )
        for group_by
        in normalized_group_by_fields
    }

    group_counts = {
        group_by: OrderedDict()
        for group_by
        in normalized_group_by_fields
    }

    row_count = 0

    try:
        with meta_file.open(
            "r",
            newline="",
            encoding="utf-8-sig",
        ) as file_obj:
            reader = csv.DictReader(
                file_obj
            )

            if not reader.fieldnames:
                raise (
                    SCSTDatasetAnnotationMetadataError(
                        "SC/ST dataset metadata file "
                        "is empty or missing a header."
                    )
                )

            normalized_fieldnames = [
                str(column).strip()
                for column
                in reader.fieldnames
            ]

            reader.fieldnames = (
                normalized_fieldnames
            )

            missing_meta_columns = [
                (
                    group_by,
                    meta_column_by_group_by[
                        group_by
                    ],
                )
                for group_by
                in normalized_group_by_fields
                if (
                    meta_column_by_group_by[
                        group_by
                    ]
                    not in normalized_fieldnames
                )
            ]

            if missing_meta_columns:
                missing_text = ", ".join(
                    (
                        f"{group_by} -> "
                        f"{meta_column}"
                    )
                    for (
                        group_by,
                        meta_column,
                    ) in missing_meta_columns
                )

                raise (
                    SCSTDatasetAnnotationMetadataError(
                        "SC/ST dataset metadata file "
                        "is missing mapped group-by "
                        "column(s): "
                        f"{missing_text}."
                    )
                )

            for row_number, row in enumerate(
                    reader,
                    start=2,
            ):
                row_count += 1

                for group_by in normalized_group_by_fields:
                    meta_column = (
                        meta_column_by_group_by[group_by]
                    )

                    group_value = (
                        normalize_scst_group_value(
                            row.get(meta_column)
                        )
                    )

                    # Empty group values are valid in SC/ST metadata,
                    # but they do not form an annotation group.
                    if not group_value:
                        continue

                    counts = group_counts[group_by]

                    if group_value not in counts:
                        counts[group_value] = 0

                    counts[group_value] += 1

    except UnicodeDecodeError as exc:
        raise SCSTDatasetAnnotationMetadataError(
            "SC/ST dataset metadata must be "
            "UTF-8 encoded."
        ) from exc

    except csv.Error as exc:
        raise SCSTDatasetAnnotationMetadataError(
            "Invalid SC/ST dataset metadata CSV."
        ) from exc

    if row_count == 0:
        raise SCSTDatasetAnnotationMetadataError(
            "SC/ST dataset metadata has no data rows."
        )

    return group_counts, row_count


def get_scst_dataset_available_background_types(
    *,
    background_file: Path,
) -> list[str]:
    """
    Return only supported interaction types that are actually
    present in the SC/ST background file.

    File existence alone is not enough for Log2FC correlation
    availability.
    """
    try:
        _, df = read_log2fc_background_file_by_path(
            background_file
        )
    except (
        FileNotFoundError,
        WorkflowLog2FCBackgroundInputError,
    ):
        return []

    return get_workflow_available_background_types(
        df=df,
        valid_types=(
            SCST_HYBRID_REFERENCE_VALID_BACKGROUND_TYPES
        ),
    )


def build_scst_group_value_visualization_availability(
    *,
    annotation_dir: Path,
    dataset_name: str,
    group_value: str,
) -> dict:
    """
    Build visualization availability for one SC/ST group value.

    Expected resources:
        annotation_network:
            {dataset}_ceRNA_{group_value}_axis.csv
            {dataset}_map_immune_axis_{group_value}.csv

        axis_final:
            {dataset}_ceRNA_{group_value}_axis_final.csv

        cmap:
            {dataset}_CMap_{group_value}.csv

        volcano:
            all:
                {dataset}_deg_{group_value}.csv

            intersect:
                {dataset}_mRNA_deg_{group_value}_intersect.csv

            RNA type:
                mRNA only

        log2fc_correlation:
            {dataset}_ceRNA_{group_value}_background.csv

        exp_correlation:
            {dataset}_ceRNA_corr_{group_value}.csv

        survival:
            {dataset}_survival_analysis_{group_value}.csv

        deg_pathway:
            {dataset}_mRNA_gsea_{group_value}.csv

        CMdrug:
            {dataset}_CMdrug_result_{group_value}/
            available only when at least one direct child file exists
    """
    paths = build_scst_group_value_result_paths(
        annotation_dir=annotation_dir,
        dataset_name=dataset_name,
        group_value=group_value,
    )

    cerna_axis_available = (
        is_scst_existing_file(
            paths["cerna_axis_file"]
        )
    )
    map_immune_axis_available = (
        is_scst_existing_file(
            paths["map_immune_axis_file"]
        )
    )

    deg_all_available = (
        is_scst_existing_file(
            paths["deg_file"]
        )
    )

    deg_intersect_available = (
        is_scst_existing_file(
            paths["deg_intersect_file"]
        )
    )

    # Keep the same semantics as the existing SC/ST workflow:
    # RNA-type availability is determined from the default/all
    # DEG scope. SC/ST currently supports mRNA only.
    available_deg_rna_types = (
        ["mRNA"]
        if deg_all_available
        else []
    )

    available_deg_scopes = []

    if available_deg_rna_types:
        if (
            SCST_WORKFLOW_DEG_SCOPE_ALL
            in SCST_WORKFLOW_DEG_SCOPES
            and deg_all_available
        ):
            available_deg_scopes.append(
                SCST_WORKFLOW_DEG_SCOPE_ALL
            )

        if (
            "intersect"
            in SCST_WORKFLOW_DEG_SCOPES
            and deg_intersect_available
        ):
            available_deg_scopes.append(
                "intersect"
            )

    available_background_types = (
        get_scst_dataset_available_background_types(
            background_file=(
                paths[
                    "log2fc_background_file"
                ]
            ),
        )
    )

    visualizations = {
        "annotation_network": {
            "available": bool(
                cerna_axis_available
                and map_immune_axis_available
            ),
            "network_source_task_type": (
                "SCSTHybridReferenceTask"
            ),
        },
        "axis_final": {
            "available": is_scst_existing_file(
                paths["axis_final_file"]
            ),
        },
        "cmap": {
            "available": is_scst_existing_file(
                paths["cmap_file"]
            ),
        },
        "volcano": {
            "available": bool(
                available_deg_rna_types
                and available_deg_scopes
            ),
            "default_rna_type": "mRNA",
            "default_deg_scope": (
                SCST_WORKFLOW_DEG_SCOPE_ALL
            ),
            "available_deg_rna_types": (
                available_deg_rna_types
            ),
            "available_deg_scopes": (
                available_deg_scopes
            ),
        },
        "log2fc_correlation": {
            "available": bool(
                available_background_types
            ),
            "available_background_types": (
                available_background_types
            ),
        },
        "exp_correlation": {
            "available": is_scst_existing_file(
                paths["exp_correlation_file"]
            ),
        },
        "survival": {
            "available": is_scst_existing_file(
                paths["survival_file"]
            ),
        },
        "deg_pathway": {
            "available": is_scst_existing_file(
                paths["mrna_gsea_file"]
            ),
        },
        "CMdrug": {
            "available": (
                is_scst_non_empty_file_directory(
                    paths["cmdrug_result_dir"]
                )
            ),
        },
    }

    available_visualization_count = sum(
        1
        for visualization
        in visualizations.values()
        if visualization["available"]
    )

    return {
        "available": (
            available_visualization_count > 0
        ),
        "available_visualization_count": (
            available_visualization_count
        ),
        "visualizations": visualizations,
    }


def build_scst_dataset_group_value_option(
    *,
    annotation_dir: Path,
    dataset_name: str,
    group_value: str,
    count: int,
) -> dict:
    availability = (
        build_scst_group_value_visualization_availability(
            annotation_dir=annotation_dir,
            dataset_name=dataset_name,
            group_value=group_value,
        )
    )

    return {
        "value": group_value,
        "label": group_value,
        "count": count,
        **availability,
    }


def build_scst_dataset_group_by_option(
    *,
    dataset_name: str,
    data_type: str,
    group_by: str,
    group_counts: OrderedDict[str, int],
) -> dict:
    """
    Build one group-by option.

    ``available`` now means that:
    - the group-by annotation directory exists, and
    - at least one metadata group value has at least one
      available visualization.
    """
    annotation_dir = (
        resolve_scst_dataset_group_annotation_dir(
            dataset_name=dataset_name,
            group_by=group_by,
            data_type=data_type,
        )
    )

    annotation_dir_available = (
        annotation_dir.exists()
        and annotation_dir.is_dir()
    )

    filtered_group_counts = OrderedDict(
        (
            group_value,
            count,
        )
        for (
            group_value,
            count,
        ) in group_counts.items()
        if not should_skip_scst_group_value_for_results(
            group_value
        )
    )

    group_values = list(
        filtered_group_counts.keys()
    )

    # Values containing "/" are intentionally omitted because
    # the upstream SC/ST annotation workflow does not generate
    # result files for them.
    #
    # Build the same visualization schema even when the group-by
    # annotation directory does not exist. Path-based checks will
    # naturally resolve every visualization to unavailable.
    group_value_options = [
        build_scst_dataset_group_value_option(
            annotation_dir=annotation_dir,
            dataset_name=dataset_name,
            group_value=group_value,
            count=count,
        )
        for (
            group_value,
            count,
        ) in filtered_group_counts.items()
    ]

    available_group_value_options = [
        option
        for option in group_value_options
        if option["available"]
    ]

    return {
        "value": group_by,
        "label": group_by,
        "available": bool(
            available_group_value_options
        ),
        "annotation_dir_available": (
            annotation_dir_available
        ),
        "group_value_count": len(
            group_values
        ),
        "skipped_group_value_count": (
            len(group_counts)
            - len(filtered_group_counts)
        ),
        "available_group_value_count": len(
            available_group_value_options
        ),
        "default_group_value": (
            available_group_value_options[0][
                "value"
            ]
            if available_group_value_options
            else None
        ),
        "group_values": group_values,
        "group_value_options": (
            group_value_options
        ),
    }


def build_scst_dataset_annotation_availability(
    *,
    dataset_name: str,
    data_type: str,
) -> dict:
    """
    Build SC/ST Dataset Annotation availability.

    Hierarchy:
        dataset
            -> group_by
                -> group_value
                    -> visualizations

    Overall availability is derived upward from visualization
    availability rather than directory existence alone.
    """
    data_type = validate_scst_data_type(
        data_type
    )

    group_by_fields = (
        get_scst_dataset_group_by_fields(
            dataset_name
        )
    )

    if not group_by_fields:
        return {
            "available": False,
            "id_column": (
                get_scst_dataset_id_column(
                    data_type
                )
            ),
            "sample_count": 0,
            "configured_group_by_count": 0,
            "available_group_by_count": 0,
            "default_group_by": None,
            "group_by_options": [],
        }

    meta_file = (
        get_scst_dataset_meta_file_path(
            dataset_name=dataset_name,
            data_type=data_type,
        )
    )

    (
        group_counts_by_group_by,
        sample_count,
    ) = read_scst_dataset_group_value_counts(
        meta_file=meta_file,
        group_by_fields=group_by_fields,
    )

    group_by_options = [
        build_scst_dataset_group_by_option(
            dataset_name=dataset_name,
            data_type=data_type,
            group_by=group_by,
            group_counts=(
                group_counts_by_group_by[
                    group_by
                ]
            ),
        )
        for group_by in group_by_fields
    ]

    available_group_by_options = [
        option
        for option in group_by_options
        if option["available"]
    ]

    return {
        "available": bool(
            available_group_by_options
        ),
        "id_column": (
            get_scst_dataset_id_column(
                data_type
            )
        ),
        "sample_count": sample_count,
        "configured_group_by_count": len(
            group_by_options
        ),
        "available_group_by_count": len(
            available_group_by_options
        ),
        "default_group_by": (
            available_group_by_options[0][
                "value"
            ]
            if available_group_by_options
            else None
        ),
        "group_by_options": group_by_options,
    }
