from pathlib import Path

import numpy as np
import pandas as pd


WORKFLOW_DEG_BASE_REQUIRED_COLUMNS = {
    "gene_name",
    "log2FC",
    "regulation",
}

WORKFLOW_DEG_VALID_REGULATION_GROUPS = [
    "NotSig",
    "Down",
    "Up",
]

WORKFLOW_DEG_PAIRED_COHORT_RNA_TYPES = [
    "mRNA",
    "miRNA",
    "lncRNA",
    "circRNA",
]

WORKFLOW_DEG_HYBRID_REFERENCE_RNA_TYPES = [
    "mRNA",
]

WORKFLOW_DEG_SCOPE_ALL = "all"
WORKFLOW_DEG_SCOPE_INTERSECT = "intersect"


# p-value == 0 cannot be transformed directly with -log10.  For volcano
# visualization, zero p-values are placed above the most significant finite
# p-value while the original p-value itself remains 0.
ZERO_PVALUE_NEG_LOG10_OFFSET = 1.0
ZERO_PVALUE_FALLBACK_NEG_LOG10 = 300.0

WORKFLOW_DEG_PAIRED_COHORT_SCOPES = [
    WORKFLOW_DEG_SCOPE_ALL,
]

WORKFLOW_DEG_HYBRID_REFERENCE_SCOPES = [
    WORKFLOW_DEG_SCOPE_ALL,
    WORKFLOW_DEG_SCOPE_INTERSECT,
]


SCST_WORKFLOW_DEG_RNA_TYPES = [
    "mRNA",
]

SCST_WORKFLOW_DEG_SCOPE_ALL = "all"
SCST_WORKFLOW_DEG_SCOPE_INTERSECT = "intersect"

SCST_WORKFLOW_DEG_SCOPES = [
    SCST_WORKFLOW_DEG_SCOPE_ALL,
    SCST_WORKFLOW_DEG_SCOPE_INTERSECT,
]

SCST_WORKFLOW_DEG_FILENAME_TEMPLATE = (
    "{task_name}_deg_{group_value}.csv"
)

SCST_WORKFLOW_DEG_INTERSECT_FILENAME_TEMPLATE = (
    "{task_name}_mRNA_deg_{group_value}_intersect.csv"
)


class WorkflowDEGVolcanoInputError(ValueError):
    pass


class WorkflowDEGVolcanoPathError(ValueError):
    pass


def validate_safe_name(value: str, field_name: str) -> None:
    value = str(value or "").strip()

    if not value:
        raise WorkflowDEGVolcanoPathError(
            f"Missing required parameter: {field_name}."
        )

    if "/" in value or "\\" in value or ".." in value:
        raise WorkflowDEGVolcanoPathError(
            f"Invalid {field_name} parameter."
        )


def get_workflow_task_output_dir(task) -> Path:
    return Path(task.get_output_dir_absolute_path()).resolve()


def get_workflow_deg_filename(
    *,
    task_name: str,
    deg_method: str,
    rna_type: str,
    deg_scope: str = WORKFLOW_DEG_SCOPE_ALL,
) -> str:
    validate_safe_name(task_name, "task_name")
    validate_safe_name(deg_method, "deg_method")
    validate_safe_name(rna_type, "rna_type")
    validate_safe_name(deg_scope, "deg_scope")

    if deg_scope == WORKFLOW_DEG_SCOPE_ALL:
        return f"{task_name}_{deg_method}_{rna_type}.csv"

    if deg_scope == WORKFLOW_DEG_SCOPE_INTERSECT:
        return f"{task_name}_{deg_method}_{rna_type}_intersect.csv"

    raise WorkflowDEGVolcanoPathError(
        f"Invalid DEG scope: {deg_scope}."
    )


def get_workflow_deg_file_path(
    *,
    task,
    rna_type: str,
    deg_scope: str = WORKFLOW_DEG_SCOPE_ALL,
) -> Path:
    task_name = str(task.task_name).strip()
    deg_method = str(task.deg_method).strip()
    rna_type = str(rna_type).strip()
    deg_scope = str(deg_scope or WORKFLOW_DEG_SCOPE_ALL).strip()

    output_dir = get_workflow_task_output_dir(task)

    filename = get_workflow_deg_filename(
        task_name=task_name,
        deg_method=deg_method,
        rna_type=rna_type,
        deg_scope=deg_scope,
    )

    file_path = (output_dir / filename).resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise WorkflowDEGVolcanoPathError(
            "Invalid workflow DEG file path."
        )

    return file_path


def validate_workflow_deg_file(
    *,
    task,
    rna_type: str,
    deg_scope: str = WORKFLOW_DEG_SCOPE_ALL,
) -> Path:
    file_path = get_workflow_deg_file_path(
        task=task,
        rna_type=rna_type,
        deg_scope=deg_scope,
    )

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"DEG file not found: {file_path.name}."
        )

    return file_path


def read_workflow_deg_file(
    *,
    task,
    rna_type: str,
    deg_scope: str = WORKFLOW_DEG_SCOPE_ALL,
) -> tuple[Path, pd.DataFrame]:
    file_path = validate_workflow_deg_file(
        task=task,
        rna_type=rna_type,
        deg_scope=deg_scope,
    )

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise WorkflowDEGVolcanoInputError(
            f"Failed to read DEG file: {str(e)}"
        )

    return file_path, df


def get_workflow_available_deg_rna_types(
    *,
    task,
    valid_rna_types: list[str],
    deg_scope: str = WORKFLOW_DEG_SCOPE_ALL,
) -> list[str]:
    output_dir = get_workflow_task_output_dir(task)

    if not output_dir.exists() or not output_dir.is_dir():
        return []

    available_rna_types = []

    for rna_type in valid_rna_types:
        try:
            file_path = get_workflow_deg_file_path(
                task=task,
                rna_type=rna_type,
                deg_scope=deg_scope,
            )
        except WorkflowDEGVolcanoPathError:
            continue

        if file_path.exists() and file_path.is_file():
            available_rna_types.append(rna_type)

    return available_rna_types


def get_workflow_available_deg_scopes(
    *,
    task,
    rna_type: str,
    valid_scopes: list[str],
) -> list[str]:
    output_dir = get_workflow_task_output_dir(task)

    if not output_dir.exists() or not output_dir.is_dir():
        return []

    available_scopes = []

    for deg_scope in valid_scopes:
        try:
            file_path = get_workflow_deg_file_path(
                task=task,
                rna_type=rna_type,
                deg_scope=deg_scope,
            )
        except WorkflowDEGVolcanoPathError:
            continue

        if file_path.exists() and file_path.is_file():
            available_scopes.append(deg_scope)

    return available_scopes


def normalize_workflow_deg_dataframe(
    *,
    df: pd.DataFrame,
    pvalue_source: str,
) -> tuple[pd.DataFrame, int, int, int]:
    required_columns = WORKFLOW_DEG_BASE_REQUIRED_COLUMNS | {
        pvalue_source,
    }

    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise WorkflowDEGVolcanoInputError(
            "Missing required columns: "
            f"{sorted(missing_columns)}."
        )

    normalized_df = df[
        [
            "gene_name",
            "log2FC",
            pvalue_source,
            "regulation",
        ]
    ].copy()

    normalized_df = normalized_df.rename(
        columns={
            pvalue_source: "pvalue",
        }
    )

    normalized_df = normalized_df.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    raw_count = int(normalized_df.shape[0])

    normalized_df["gene_name"] = (
        normalized_df["gene_name"]
        .astype(str)
        .str.strip()
    )

    normalized_df["regulation"] = (
        normalized_df["regulation"]
        .astype(str)
        .str.strip()
    )

    normalized_df["log2FC"] = pd.to_numeric(
        normalized_df["log2FC"],
        errors="coerce",
    )

    normalized_df["pvalue"] = pd.to_numeric(
        normalized_df["pvalue"],
        errors="coerce",
    )

    normalized_df = normalized_df.dropna(
        subset=[
            "gene_name",
            "log2FC",
            "pvalue",
            "regulation",
        ]
    )

    normalized_df = normalized_df[
        normalized_df["gene_name"] != ""
    ]

    # Keep p-value == 0.  These values commonly arise from numerical
    # underflow and should not be discarded from a volcano plot.
    # Negative values and values > 1 remain invalid.
    normalized_df = normalized_df[
        normalized_df["pvalue"] >= 0
    ]

    normalized_df = normalized_df[
        normalized_df["pvalue"] <= 1
    ]

    normalized_df = normalized_df[
        normalized_df["regulation"].isin(
            WORKFLOW_DEG_VALID_REGULATION_GROUPS
        )
    ].copy()

    cleaned_count = int(normalized_df.shape[0])
    dropped_count = raw_count - cleaned_count

    normalized_df["is_zero_pvalue"] = (
        normalized_df["pvalue"] == 0
    )

    positive_mask = normalized_df["pvalue"] > 0
    zero_mask = normalized_df["is_zero_pvalue"]

    # First calculate the real -log10(p) values for all finite positive
    # p-values.  The zero p-values are assigned a separate plotting height
    # below, without changing their original pvalue field.
    normalized_df["neg_log10_pvalue"] = np.nan
    normalized_df.loc[
        positive_mask,
        "neg_log10_pvalue",
    ] = -np.log10(
        normalized_df.loc[positive_mask, "pvalue"]
    )

    positive_pvalues = normalized_df.loc[
        positive_mask,
        "pvalue",
    ]

    if not positive_pvalues.empty:
        min_positive_pvalue = float(positive_pvalues.min())
        max_finite_neg_log10 = float(
            -np.log10(min_positive_pvalue)
        )
        zero_pvalue_neg_log10 = (
            max_finite_neg_log10
            + ZERO_PVALUE_NEG_LOG10_OFFSET
        )
    else:
        min_positive_pvalue = None
        max_finite_neg_log10 = None
        zero_pvalue_neg_log10 = (
            ZERO_PVALUE_FALLBACK_NEG_LOG10
        )

    normalized_df.loc[
        zero_mask,
        "neg_log10_pvalue",
    ] = zero_pvalue_neg_log10

    # Store zero-p-value plotting metadata separately from the original
    # statistical values.  The response builder exposes this metadata to the
    # frontend so it can draw the dedicated "extremely significant" line.
    normalized_df.attrs["zero_pvalue_plot"] = {
        "count": int(zero_mask.sum()),
        "min_positive_pvalue": min_positive_pvalue,
        "max_finite_neg_log10_pvalue": max_finite_neg_log10,
        "neg_log10_offset": ZERO_PVALUE_NEG_LOG10_OFFSET,
        "neg_log10_plot_y": float(zero_pvalue_neg_log10),
        "used_fallback": min_positive_pvalue is None,
    }

    return normalized_df, raw_count, cleaned_count, dropped_count


def build_workflow_deg_volcano_groups(
    df: pd.DataFrame,
) -> dict:
    groups = {}

    for group in WORKFLOW_DEG_VALID_REGULATION_GROUPS:
        sub_df = df[df["regulation"] == group]

        groups[group] = [
            {
                "gene_name": row["gene_name"],
                "log2FC": float(row["log2FC"]),
                "pvalue": float(row["pvalue"]),
                "neg_log10_pvalue": float(row["neg_log10_pvalue"]),
                "is_zero_pvalue": bool(row["is_zero_pvalue"]),
            }
            for _, row in sub_df.iterrows()
        ]

    return groups


def build_workflow_deg_volcano_response_data(
    *,
    task,
    task_type: str,
    rna_type: str,
    deg_scope: str,
    deg_file_name: str,
    df: pd.DataFrame,
    use_padj: bool,
) -> dict:
    base_response = {
        "uuid": str(task.uuid),
        "task_type": task_type,
        "task_name": task.task_name,
    }

    return build_deg_volcano_response_data_from_dataframe(
        df=df,
        deg_file_name=deg_file_name,
        rna_type=rna_type,
        deg_scope=deg_scope,
        deg_method=getattr(task, "deg_method", None),
        use_padj=use_padj,
        base_response=base_response,
    )


def validate_deg_file_path(file_path: Path) -> Path:
    file_path = Path(file_path).resolve()

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"DEG file not found: {file_path.name}."
        )

    return file_path


def read_deg_file_by_path(file_path: Path) -> tuple[Path, pd.DataFrame]:
    file_path = validate_deg_file_path(file_path)

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise WorkflowDEGVolcanoInputError(
            f"Failed to read DEG file: {str(e)}"
        )

    return file_path, df


def build_deg_volcano_response_data_from_dataframe(
    *,
    df: pd.DataFrame,
    deg_file_name: str,
    rna_type: str,
    deg_scope: str,
    deg_method: str,
    use_padj: bool,
    base_response: dict | None = None,
) -> dict:
    pvalue_source = "padj" if use_padj else "pvalue"
    pvalue_label = "adjusted p-value" if use_padj else "raw p-value"

    try:
        volcano_df, raw_count, cleaned_count, dropped_count = (
            normalize_workflow_deg_dataframe(
                df=df,
                pvalue_source=pvalue_source,
            )
        )
    except WorkflowDEGVolcanoInputError as e:
        raise WorkflowDEGVolcanoInputError(
            f"{str(e)} This DEG volcano view uses {pvalue_label}, "
            f"so the DEG file must contain column: {pvalue_source}."
        )

    groups = build_workflow_deg_volcano_groups(volcano_df)
    zero_pvalue_plot = volcano_df.attrs.get(
        "zero_pvalue_plot",
        {},
    )

    response_data = {
        "deg_method": deg_method,
        "rna_type": rna_type,
        "deg_scope": deg_scope,
        "deg_file": deg_file_name,

        "use_padj": use_padj,
        "pvalue_source": pvalue_source,
        "pvalue_label": pvalue_label,
        "zero_pvalue_plot": zero_pvalue_plot,

        "summary": {
            "raw_count": raw_count,
            "cleaned_count": cleaned_count,
            "dropped_count": dropped_count,
            "not_sig": len(groups["NotSig"]),
            "down": len(groups["Down"]),
            "up": len(groups["Up"]),
        },
        "groups": groups,
    }

    if base_response:
        response_data = {
            **base_response,
            **response_data,
        }

    return response_data


def validate_scst_deg_group_value(
    group_value: str,
) -> str:
    """
    Validate a cell type/group value used in an SC/ST DEG filename.

    Spaces, parentheses, plus signs and similar label characters are
    permitted. Path separators and traversal sequences are forbidden.
    """
    group_value = str(
        group_value or ""
    ).strip()

    if not group_value:
        raise WorkflowDEGVolcanoPathError(
            "Missing required parameter: group_value."
        )

    if "\x00" in group_value:
        raise WorkflowDEGVolcanoPathError(
            "Invalid group_value parameter."
        )

    if (
        "/" in group_value
        or "\\" in group_value
        or ".." in group_value
    ):
        raise WorkflowDEGVolcanoPathError(
            "Invalid group_value parameter."
        )

    return group_value


def validate_scst_deg_scope(
    deg_scope: str,
) -> str:
    """
    Validate an SC/ST DEG scope.
    """
    deg_scope = str(
        deg_scope or SCST_WORKFLOW_DEG_SCOPE_ALL
    ).strip()

    if deg_scope not in SCST_WORKFLOW_DEG_SCOPES:
        raise WorkflowDEGVolcanoPathError(
            "Invalid SC/ST DEG scope. "
            "Allowed values are: "
            f"{', '.join(SCST_WORKFLOW_DEG_SCOPES)}."
        )

    return deg_scope


def get_scst_workflow_deg_filename(
    *,
    task_name: str,
    group_value: str,
    rna_type: str = "mRNA",
    deg_scope: str = SCST_WORKFLOW_DEG_SCOPE_ALL,
) -> str:
    """
    Build an SC/ST Hybrid Reference DEG filename.

    all:
        {task_name}_deg_{group_value}.csv

    intersect:
        {dataset}_mRNA_deg_{group_value}_intersect.csv
    """
    task_name = str(
        task_name or ""
    ).strip()

    rna_type = str(
        rna_type or ""
    ).strip()

    validate_safe_name(
        task_name,
        "task_name",
    )

    validate_safe_name(
        rna_type,
        "rna_type",
    )

    group_value = validate_scst_deg_group_value(
        group_value
    )

    deg_scope = validate_scst_deg_scope(
        deg_scope
    )

    if rna_type not in SCST_WORKFLOW_DEG_RNA_TYPES:
        raise WorkflowDEGVolcanoPathError(
            f"Unsupported SC/ST DEG RNA type: {rna_type}."
        )

    if deg_scope == SCST_WORKFLOW_DEG_SCOPE_ALL:
        return SCST_WORKFLOW_DEG_FILENAME_TEMPLATE.format(
            task_name=task_name,
            group_value=group_value,
        )

    return (
        SCST_WORKFLOW_DEG_INTERSECT_FILENAME_TEMPLATE.format(
            task_name=task_name,
            group_value=group_value,
        )
    )


def get_scst_workflow_deg_file_path(
    *,
    task,
    group_value: str,
    rna_type: str = "mRNA",
    deg_scope: str = SCST_WORKFLOW_DEG_SCOPE_ALL,
) -> Path:
    """
    Return the SC/ST DEG file path for one group and scope.
    """
    output_dir = get_workflow_task_output_dir(
        task
    )

    filename = get_scst_workflow_deg_filename(
        task_name=task.task_name,
        group_value=group_value,
        rna_type=rna_type,
        deg_scope=deg_scope,
    )

    file_path = (
        output_dir / filename
    ).resolve()

    try:
        file_path.relative_to(output_dir)
    except ValueError as exc:
        raise WorkflowDEGVolcanoPathError(
            "Invalid SC/ST workflow DEG file path."
        ) from exc

    return file_path


def validate_scst_workflow_deg_file(
    *,
    task,
    group_value: str,
    rna_type: str = "mRNA",
    deg_scope: str = SCST_WORKFLOW_DEG_SCOPE_ALL,
) -> Path:
    """
    Validate that an SC/ST DEG result file exists.
    """
    file_path = get_scst_workflow_deg_file_path(
        task=task,
        group_value=group_value,
        rna_type=rna_type,
        deg_scope=deg_scope,
    )

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"SC/ST DEG file not found: {file_path.name}."
        )

    return file_path


def read_scst_workflow_deg_file(
    *,
    task,
    group_value: str,
    rna_type: str = "mRNA",
    deg_scope: str = SCST_WORKFLOW_DEG_SCOPE_ALL,
) -> tuple[Path, pd.DataFrame]:
    """
    Read one SC/ST DEG result file.
    """
    file_path = validate_scst_workflow_deg_file(
        task=task,
        group_value=group_value,
        rna_type=rna_type,
        deg_scope=deg_scope,
    )

    return read_deg_file_by_path(
        file_path
    )


def get_scst_workflow_available_deg_rna_types(
    *,
    task,
    group_value: str,
    valid_rna_types: list[str] | None = None,
    deg_scope: str = SCST_WORKFLOW_DEG_SCOPE_ALL,
) -> list[str]:
    """
    Return available DEG RNA types for one SC/ST group and scope.

    The current SC/ST workflow supports only mRNA.
    """
    valid_rna_types = (
        valid_rna_types
        or SCST_WORKFLOW_DEG_RNA_TYPES
    )

    try:
        deg_scope = validate_scst_deg_scope(
            deg_scope
        )
    except WorkflowDEGVolcanoPathError:
        return []

    output_dir = get_workflow_task_output_dir(
        task
    )

    if not output_dir.exists() or not output_dir.is_dir():
        return []

    available_rna_types = []

    for rna_type in valid_rna_types:
        try:
            file_path = get_scst_workflow_deg_file_path(
                task=task,
                group_value=group_value,
                rna_type=rna_type,
                deg_scope=deg_scope,
            )
        except WorkflowDEGVolcanoPathError:
            continue

        if file_path.exists() and file_path.is_file():
            available_rna_types.append(
                rna_type
            )

    return available_rna_types


def get_scst_workflow_available_deg_scopes(
    *,
    task,
    group_value: str,
    rna_type: str = "mRNA",
    valid_scopes: list[str] | None = None,
) -> list[str]:
    """
    Return available DEG scopes for one SC/ST group.

    Possible scopes:
        all
        intersect
    """
    valid_scopes = (
        valid_scopes
        or SCST_WORKFLOW_DEG_SCOPES
    )

    output_dir = get_workflow_task_output_dir(
        task
    )

    if not output_dir.exists() or not output_dir.is_dir():
        return []

    available_scopes = []

    for deg_scope in valid_scopes:
        try:
            file_path = get_scst_workflow_deg_file_path(
                task=task,
                group_value=group_value,
                rna_type=rna_type,
                deg_scope=deg_scope,
            )
        except WorkflowDEGVolcanoPathError:
            continue

        if file_path.exists() and file_path.is_file():
            available_scopes.append(
                deg_scope
            )

    return available_scopes


def get_scst_workflow_deg_availability(
    *,
    task,
    group_values: list[str],
    valid_rna_types: list[str] | None = None,
    valid_scopes: list[str] | None = None,
) -> list[dict]:
    """
    Return DEG availability information for every SC/ST group.

    Each group reports:
    - Available RNA types
    - Available scopes
    - File information for all and intersect scopes
    """
    valid_rna_types = (
        valid_rna_types
        or SCST_WORKFLOW_DEG_RNA_TYPES
    )

    valid_scopes = (
        valid_scopes
        or SCST_WORKFLOW_DEG_SCOPES
    )

    results = []

    for raw_group_value in group_values:
        group_value = str(
            raw_group_value or ""
        ).strip()

        if not group_value:
            continue

        available_scopes = (
            get_scst_workflow_available_deg_scopes(
                task=task,
                group_value=group_value,
                rna_type="mRNA",
                valid_scopes=valid_scopes,
            )
        )

        # available_deg_rna_types 表示至少存在一个可用 DEG scope
        # 时，对应 RNA type 可以用于 DEG 可视化。
        available_rna_types = []

        for rna_type in valid_rna_types:
            rna_type_available = any(
                rna_type in (
                    get_scst_workflow_available_deg_rna_types(
                        task=task,
                        group_value=group_value,
                        valid_rna_types=[rna_type],
                        deg_scope=deg_scope,
                    )
                )
                for deg_scope in valid_scopes
            )

            if rna_type_available:
                available_rna_types.append(
                    rna_type
                )

        deg_files = {}

        for deg_scope in valid_scopes:
            try:
                file_path = get_scst_workflow_deg_file_path(
                    task=task,
                    group_value=group_value,
                    rna_type="mRNA",
                    deg_scope=deg_scope,
                )

                file_exists = (
                    file_path.exists()
                    and file_path.is_file()
                )

                deg_files[deg_scope] = {
                    "file": (
                        file_path.name
                        if file_exists
                        else None
                    ),
                    "exists": file_exists,
                }

            except WorkflowDEGVolcanoPathError:
                deg_files[deg_scope] = {
                    "file": None,
                    "exists": False,
                }

        results.append(
            {
                "group_value": group_value,
                "deg_available": bool(
                    available_rna_types
                ),
                "available_deg_rna_types": (
                    available_rna_types
                ),
                "available_deg_scopes": (
                    available_scopes
                ),
                "deg_files": deg_files,
            }
        )

    return results
