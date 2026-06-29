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

WORKFLOW_DEG_PAIRED_COHORT_SCOPES = [
    WORKFLOW_DEG_SCOPE_ALL,
]

WORKFLOW_DEG_HYBRID_REFERENCE_SCOPES = [
    WORKFLOW_DEG_SCOPE_ALL,
    WORKFLOW_DEG_SCOPE_INTERSECT,
]


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

    normalized_df = normalized_df[
        normalized_df["pvalue"] > 0
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

    normalized_df["neg_log10_pvalue"] = -np.log10(
        normalized_df["pvalue"]
    )

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
        deg_method=task.deg_method,
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

    response_data = {
        "deg_method": deg_method,
        "rna_type": rna_type,
        "deg_scope": deg_scope,
        "deg_file": deg_file_name,

        "use_padj": use_padj,
        "pvalue_source": pvalue_source,
        "pvalue_label": pvalue_label,

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
