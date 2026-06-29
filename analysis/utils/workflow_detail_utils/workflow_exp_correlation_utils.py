from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd


WORKFLOW_EXP_CORRELATION_FILENAME_SUFFIX = "_ceRNA_corr.csv"

WORKFLOW_EXP_CORRELATION_REQUIRED_COLUMNS = {
    "gene1",
    "gene2",
    "type",
    "correlation type",
    "correlation",
    "p value",
}

WORKFLOW_EXP_CORRELATION_OPTIONAL_COLUMNS = {
    "patient",
}

BASE_EXP_CORRELATION_TYPES = [
    "miRNA-mRNA",
    "miRNA-lncRNA",
    "lncRNA-mRNA",
]

CIRCRNA_EXP_CORRELATION_TYPES = [
    "miRNA-circRNA",
    "circRNA-mRNA",
]

PAIRED_COHORT_VALID_EXP_CORRELATION_TYPES = [
    *BASE_EXP_CORRELATION_TYPES,
    *CIRCRNA_EXP_CORRELATION_TYPES,
]

HYBRID_REFERENCE_VALID_EXP_CORRELATION_TYPES = [
    *BASE_EXP_CORRELATION_TYPES,
    *CIRCRNA_EXP_CORRELATION_TYPES,
]

WORKFLOW_EXP_CORRELATION_TYPE_FILE_FIELD_MAP = {
    "miRNA-mRNA": {
        "gene1_file": "mirna_file",
        "gene2_file": "mrna_file",
    },
    "miRNA-lncRNA": {
        "gene1_file": "mirna_file",
        "gene2_file": "lncrna_file",
    },
    "lncRNA-mRNA": {
        "gene1_file": "lncrna_file",
        "gene2_file": "mrna_file",
    },
    "miRNA-circRNA": {
        "gene1_file": "mirna_file",
        "gene2_file": "circrna_file",
    },
    "circRNA-mRNA": {
        "gene1_file": "circrna_file",
        "gene2_file": "mrna_file",
    },
}

WORKFLOW_EXPR_SAMPLE_COL = "sample_id"

WORKFLOW_EXP_CORRELATION_TYPE_RNA_MAP = {
    "miRNA-mRNA": {
        "gene1_rna_type": "miRNA",
        "gene2_rna_type": "mRNA",
    },
    "miRNA-lncRNA": {
        "gene1_rna_type": "miRNA",
        "gene2_rna_type": "lncRNA",
    },
    "lncRNA-mRNA": {
        "gene1_rna_type": "lncRNA",
        "gene2_rna_type": "mRNA",
    },
    "miRNA-circRNA": {
        "gene1_rna_type": "miRNA",
        "gene2_rna_type": "circRNA",
    },
    "circRNA-mRNA": {
        "gene1_rna_type": "circRNA",
        "gene2_rna_type": "mRNA",
    },
}

WORKFLOW_EXP_CORRELATION_STATS = {
    "Pearson Correlation": ("pearson_r", "pearson_p"),
    "Spearman Correlation": ("spearman_r", "spearman_p"),
    "Kendall's tau": ("kendall_tau", "kendall_p"),
}


class WorkflowExpCorrelationInputError(ValueError):
    pass


class WorkflowExpCorrelationPathError(ValueError):
    pass


def validate_safe_name(value: str, field_name: str) -> None:
    value = str(value or "").strip()

    if not value:
        raise WorkflowExpCorrelationPathError(
            f"Missing required parameter: {field_name}."
        )

    if "/" in value or "\\" in value or ".." in value:
        raise WorkflowExpCorrelationPathError(
            f"Invalid {field_name} parameter."
        )


def get_workflow_task_output_dir(task) -> Path:
    return Path(task.get_output_dir_absolute_path()).resolve()


def get_workflow_task_input_dir(task) -> Path:
    return Path(task.get_input_dir_absolute_path()).resolve()


def get_workflow_exp_correlation_file_path(task) -> Path:
    task_name = str(task.task_name).strip()

    validate_safe_name(task_name, "task_name")

    output_dir = get_workflow_task_output_dir(task)

    file_path = (
        output_dir / f"{task_name}{WORKFLOW_EXP_CORRELATION_FILENAME_SUFFIX}"
    ).resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise WorkflowExpCorrelationPathError(
            "Invalid workflow expression correlation file path."
        )

    return file_path


def validate_workflow_exp_correlation_file(task) -> Path:
    file_path = get_workflow_exp_correlation_file_path(task)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Workflow expression correlation file not found: {file_path.name}"
        )

    return file_path


def read_workflow_exp_correlation_file(task) -> tuple[Path, pd.DataFrame]:
    file_path = validate_workflow_exp_correlation_file(task)

    return read_exp_correlation_file_by_path(file_path)


def normalize_workflow_exp_correlation_dataframe(
    df: pd.DataFrame,
) -> pd.DataFrame:
    normalized_df = df.copy()

    string_columns = [
        "patient",
        "gene1",
        "gene2",
        "type",
        "correlation type",
    ]

    for col in string_columns:
        if col in normalized_df.columns:
            normalized_df[col] = (
                normalized_df[col]
                .astype(str)
                .str.strip()
            )

    for col in ["correlation", "p value"]:
        if col in normalized_df.columns:
            normalized_df[col] = pd.to_numeric(
                normalized_df[col],
                errors="coerce",
            )

    normalized_df = normalized_df.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    return normalized_df


def get_workflow_available_exp_correlation_pairs(
    df: pd.DataFrame,
    valid_types: list[str],
) -> list[dict]:
    normalized_df = normalize_workflow_exp_correlation_dataframe(df)

    normalized_df = normalized_df[
        normalized_df["type"].isin(valid_types)
    ].copy()

    pair_df = (
        normalized_df[["gene1", "gene2", "type"]]
        .dropna()
        .drop_duplicates()
        .sort_values(["type", "gene1", "gene2"])
    )

    return [
        {
            "gene1": row["gene1"],
            "gene2": row["gene2"],
            "type": row["type"],
        }
        for _, row in pair_df.iterrows()
    ]


def get_workflow_available_exp_correlation_types(
    df: pd.DataFrame,
    valid_types: list[str],
) -> list[str]:
    normalized_df = normalize_workflow_exp_correlation_dataframe(df)

    observed_types = set(
        normalized_df["type"]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )

    return [
        type_value
        for type_value in valid_types
        if type_value in observed_types
    ]


def get_workflow_exp_file_fields_by_type(type_value: str) -> dict:
    type_value = str(type_value).strip()

    if type_value not in WORKFLOW_EXP_CORRELATION_TYPE_FILE_FIELD_MAP:
        raise WorkflowExpCorrelationInputError(
            "Invalid type. Allowed values are: "
            f"{', '.join(WORKFLOW_EXP_CORRELATION_TYPE_FILE_FIELD_MAP.keys())}."
        )

    return WORKFLOW_EXP_CORRELATION_TYPE_FILE_FIELD_MAP[type_value]


def safe_float_or_none(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None

    if pd.isna(result) or not np.isfinite(result):
        return None

    return result


def extract_workflow_exp_correlation_stats(
    cor_df: pd.DataFrame,
) -> dict:
    stats = {
        "pearson_r": None,
        "pearson_p": None,
        "spearman_r": None,
        "spearman_p": None,
        "kendall_tau": None,
        "kendall_p": None,
    }

    normalized_df = normalize_workflow_exp_correlation_dataframe(cor_df)

    for correlation_type, fields in WORKFLOW_EXP_CORRELATION_STATS.items():
        value_field, pvalue_field = fields

        sub_df = normalized_df[
            normalized_df["correlation type"] == correlation_type
        ]

        if sub_df.empty:
            continue

        row = sub_df.iloc[0]

        stats[value_field] = safe_float_or_none(
            row.get("correlation")
        )
        stats[pvalue_field] = safe_float_or_none(
            row.get("p value")
        )

    return stats


def get_selected_workflow_exp_correlation_pair_df(
    *,
    cor_df: pd.DataFrame,
    gene1: str,
    gene2: str,
    type_value: str,
) -> pd.DataFrame:
    normalized_df = normalize_workflow_exp_correlation_dataframe(cor_df)

    pair_df = normalized_df[
        (normalized_df["gene1"] == gene1)
        & (normalized_df["gene2"] == gene2)
        & (normalized_df["type"] == type_value)
    ].copy()

    if pair_df.empty:
        raise WorkflowExpCorrelationInputError(
            "Selected gene pair was not found in workflow expression "
            "correlation file."
        )

    return pair_df


def calculate_regression(plot_df: pd.DataFrame) -> dict:
    if plot_df.shape[0] < 2:
        return {
            "slope": None,
            "intercept": None,
        }

    if plot_df["gene1_expr"].nunique() < 2:
        return {
            "slope": None,
            "intercept": None,
        }

    slope, intercept = np.polyfit(
        plot_df["gene1_expr"],
        plot_df["gene2_expr"],
        1,
    )

    return {
        "slope": float(slope),
        "intercept": float(intercept),
    }


def get_workflow_exp_rna_types_by_type(type_value: str) -> dict:
    type_value = str(type_value or "").strip()

    if type_value not in WORKFLOW_EXP_CORRELATION_TYPE_RNA_MAP:
        raise WorkflowExpCorrelationInputError(
            "Invalid type. Allowed values are: "
            f"{', '.join(WORKFLOW_EXP_CORRELATION_TYPE_RNA_MAP.keys())}."
        )

    return WORKFLOW_EXP_CORRELATION_TYPE_RNA_MAP[type_value]


def build_expression_pair_points_from_files(
    *,
    gene1_expr_file: Path,
    gene2_expr_file: Path,
    gene1: str,
    gene2: str,
    sample_col: str = WORKFLOW_EXPR_SAMPLE_COL,
    sample_ids: list[str] | None = None,
) -> dict:
    for file_path in [gene1_expr_file, gene2_expr_file]:
        if not file_path.exists() or not file_path.is_file():
            raise FileNotFoundError(
                f"Expression file not found: {file_path.name}."
            )

    try:
        gene1_df = pd.read_csv(
            gene1_expr_file,
            index_col=sample_col,
        )

        gene2_df = pd.read_csv(
            gene2_expr_file,
            index_col=sample_col,
        )

        gene1_df.columns = [
            str(col).strip()
            for col in gene1_df.columns
        ]

        gene2_df.columns = [
            str(col).strip()
            for col in gene2_df.columns
        ]
    except ValueError as e:
        raise WorkflowExpCorrelationInputError(
            f"Failed to read expression file by sample column "
            f"{sample_col}: {str(e)}"
        )
    except Exception as e:
        raise WorkflowExpCorrelationInputError(
            f"Failed to read expression file: {str(e)}"
        )

    if gene1 not in gene1_df.columns:
        raise WorkflowExpCorrelationInputError(
            f"gene1 not found in expression file: {gene1}."
        )

    if gene2 not in gene2_df.columns:
        raise WorkflowExpCorrelationInputError(
            f"gene2 not found in expression file: {gene2}."
        )

    gene1_expr = gene1_df[[gene1]].copy()
    gene2_expr = gene2_df[[gene2]].copy()

    gene1_expr.index = gene1_expr.index.astype(str)
    gene2_expr.index = gene2_expr.index.astype(str)

    merged_df = pd.merge(
        gene1_expr,
        gene2_expr,
        left_index=True,
        right_index=True,
        how="inner",
    )

    if sample_ids is not None:
        sample_id_set = {
            str(sample_id)
            for sample_id in sample_ids
        }

        merged_df = merged_df[
            merged_df.index.astype(str).isin(sample_id_set)
        ].copy()

    merged_df.columns = [
        "gene1_expr",
        "gene2_expr",
    ]

    merged_df["gene1_expr"] = pd.to_numeric(
        merged_df["gene1_expr"],
        errors="coerce",
    )

    merged_df["gene2_expr"] = pd.to_numeric(
        merged_df["gene2_expr"],
        errors="coerce",
    )

    raw_count = int(merged_df.shape[0])

    plot_df = (
        merged_df
        .replace([np.inf, -np.inf], np.nan)
        .dropna(subset=["gene1_expr", "gene2_expr"])
        .copy()
    )

    cleaned_count = int(plot_df.shape[0])
    dropped_count = raw_count - cleaned_count

    regression = calculate_regression(plot_df)

    points = [
        {
            "sample_id": str(sample_id),
            "gene1_expr": float(row["gene1_expr"]),
            "gene2_expr": float(row["gene2_expr"]),
        }
        for sample_id, row in plot_df.iterrows()
    ]

    return {
        "summary": {
            "raw_count": raw_count,
            "cleaned_count": cleaned_count,
            "dropped_count": dropped_count,
            "has_points": True,
            "point_unavailable_reason": "",
        },
        "regression": regression,
        "points": points,
    }


def validate_exp_correlation_file_path(file_path: Path) -> Path:
    file_path = Path(file_path).resolve()

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Expression correlation file not found: {file_path.name}"
        )

    return file_path


def read_exp_correlation_file_by_path(
    file_path: Path,
) -> tuple[Path, pd.DataFrame]:
    file_path = validate_exp_correlation_file_path(file_path)

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise WorkflowExpCorrelationInputError(
            f"Failed to read expression correlation file: {str(e)}"
        )

    missing_columns = (
        WORKFLOW_EXP_CORRELATION_REQUIRED_COLUMNS - set(df.columns)
    )

    if missing_columns:
        raise WorkflowExpCorrelationInputError(
            "Expression correlation file is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}."
        )

    return file_path, df


def get_workflow_exp_rna_types_for_types(
    valid_types: list[str],
) -> list[str]:
    rna_types = []

    for type_value in valid_types:
        if type_value not in WORKFLOW_EXP_CORRELATION_TYPE_RNA_MAP:
            continue

        type_map = WORKFLOW_EXP_CORRELATION_TYPE_RNA_MAP[type_value]

        for key in ["gene1_rna_type", "gene2_rna_type"]:
            rna_type = type_map[key]

            if rna_type not in rna_types:
                rna_types.append(rna_type)

    return rna_types


@lru_cache(maxsize=512)
def _read_expression_gene_set_cached(
    file_path_str: str,
    sample_col: str,
    mtime_ns: int,
    size: int,
) -> frozenset[str]:
    try:
        columns = pd.read_csv(
            file_path_str,
            nrows=0,
        ).columns
    except ValueError as e:
        raise WorkflowExpCorrelationInputError(
            f"Failed to read expression file header: {str(e)}"
        )
    except Exception as e:
        raise WorkflowExpCorrelationInputError(
            f"Failed to read expression file header: {str(e)}"
        )

    normalized_columns = [
        str(col).strip()
        for col in columns
    ]

    if sample_col not in normalized_columns:
        raise WorkflowExpCorrelationInputError(
            f"Expression file is missing sample column: {sample_col}."
        )

    return frozenset(
        col
        for col in normalized_columns
        if col and col != sample_col
    )


def read_expression_gene_set_from_file(
    file_path: Path,
    sample_col: str = WORKFLOW_EXPR_SAMPLE_COL,
) -> set[str]:
    file_path = Path(file_path).resolve()

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Expression file not found: {file_path.name}."
        )

    stat = file_path.stat()

    return set(
        _read_expression_gene_set_cached(
            str(file_path),
            sample_col,
            stat.st_mtime_ns,
            stat.st_size,
        )
    )


def filter_workflow_exp_correlation_pairs_by_expression_genes(
    *,
    pairs: list[dict],
    expression_gene_sets: dict[str, set[str]],
) -> list[dict]:
    filtered_pairs = []

    for pair in pairs:
        type_value = str(pair.get("type", "")).strip()
        gene1 = str(pair.get("gene1", "")).strip()
        gene2 = str(pair.get("gene2", "")).strip()

        if not type_value or not gene1 or not gene2:
            continue

        rna_type_map = get_workflow_exp_rna_types_by_type(type_value)

        gene1_rna_type = rna_type_map["gene1_rna_type"]
        gene2_rna_type = rna_type_map["gene2_rna_type"]

        gene1_set = expression_gene_sets.get(gene1_rna_type)
        gene2_set = expression_gene_sets.get(gene2_rna_type)

        if not gene1_set or not gene2_set:
            continue

        if gene1 not in gene1_set:
            continue

        if gene2 not in gene2_set:
            continue

        filtered_pairs.append(pair)

    return filtered_pairs


def get_workflow_available_exp_correlation_types_from_pairs(
    *,
    pairs: list[dict],
    valid_types: list[str],
) -> list[str]:
    observed_types = {
        str(pair.get("type", "")).strip()
        for pair in pairs
        if pair.get("type")
    }

    return [
        type_value
        for type_value in valid_types
        if type_value in observed_types
    ]
