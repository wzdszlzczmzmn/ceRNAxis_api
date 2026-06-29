import math

import numpy as np
import pandas as pd


DEG_PATHWAY_REQUIRED_COLUMNS = {
    "Term",
    "NES",
    "FDR q-val",
}

DEG_PATHWAY_OPTIONAL_COLUMNS = [
    "Name",
    "ES",
    "NOM p-val",
    "FWER p-val",
    "Tag %",
    "Gene %",
    "Lead_genes",
]


class DEGPathwayInputError(ValueError):
    pass


def safe_float_or_none(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None

    if pd.isna(result) or not math.isfinite(result):
        return None

    return result


def get_optional_value(row, column_name: str):
    if column_name not in row.index:
        return None

    value = row.get(column_name)

    if pd.isna(value):
        return None

    return value


def validate_deg_pathway_dataframe_columns(
    df: pd.DataFrame,
    required_columns: set[str] | None = None,
) -> None:
    required_columns = required_columns or DEG_PATHWAY_REQUIRED_COLUMNS

    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise DEGPathwayInputError(
            "DEG pathway file is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}."
        )


def normalize_deg_pathway_dataframe(
    df: pd.DataFrame,
) -> pd.DataFrame:
    validate_deg_pathway_dataframe_columns(df)

    normalized_df = df.copy()

    normalized_df = normalized_df.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    normalized_df["Term"] = (
        normalized_df["Term"]
        .astype(str)
        .str.strip()
    )

    normalized_df["NES"] = pd.to_numeric(
        normalized_df["NES"],
        errors="coerce",
    )

    normalized_df["FDR q-val"] = pd.to_numeric(
        normalized_df["FDR q-val"],
        errors="coerce",
    )

    if "ES" in normalized_df.columns:
        normalized_df["ES"] = pd.to_numeric(
            normalized_df["ES"],
            errors="coerce",
        )

    if "NOM p-val" in normalized_df.columns:
        normalized_df["NOM p-val"] = pd.to_numeric(
            normalized_df["NOM p-val"],
            errors="coerce",
        )

    if "FWER p-val" in normalized_df.columns:
        normalized_df["FWER p-val"] = pd.to_numeric(
            normalized_df["FWER p-val"],
            errors="coerce",
        )

    normalized_df = normalized_df.dropna(
        subset=[
            "Term",
            "NES",
            "FDR q-val",
        ]
    )

    normalized_df = normalized_df[
        normalized_df["Term"] != ""
    ].copy()

    normalized_df = normalized_df[
        normalized_df["FDR q-val"] > 0
    ].copy()

    normalized_df = normalized_df[
        normalized_df["FDR q-val"] <= 1
    ].copy()

    normalized_df["neg_log10_fdr_qval"] = -np.log10(
        normalized_df["FDR q-val"]
    )

    normalized_df = normalized_df.sort_values(
        by="NES",
        ascending=False,
    )

    return normalized_df


def format_deg_pathway_row(row) -> dict:
    return {
        "term": row["Term"],
        "nes": safe_float_or_none(row["NES"]),
        "fdr_qval": safe_float_or_none(row["FDR q-val"]),
        "neg_log10_fdr_qval": safe_float_or_none(
            row["neg_log10_fdr_qval"]
        ),

        "name": get_optional_value(row, "Name"),
        "es": safe_float_or_none(
            get_optional_value(row, "ES")
        ),
        "nom_pval": safe_float_or_none(
            get_optional_value(row, "NOM p-val")
        ),
        "fwer_pval": safe_float_or_none(
            get_optional_value(row, "FWER p-val")
        ),
        "tag_percent": get_optional_value(row, "Tag %"),
        "gene_percent": get_optional_value(row, "Gene %"),
        "lead_genes": get_optional_value(row, "Lead_genes"),
    }


def build_deg_pathway_data_from_dataframe(
    *,
    task,
    gsea_file_name: str,
    df: pd.DataFrame,
    title: str = "DEG Pathway Enrichment",
) -> dict:
    return build_deg_pathway_data_from_dataframe_common(
        gsea_file_name=gsea_file_name,
        df=df,
        title=title,
        base_response={
            "uuid": str(task.uuid),
            "task_type": task.__class__.__name__,
            "task_name": task.task_name,
        },
    )


def build_deg_pathway_data_from_dataframe_common(
    *,
    gsea_file_name: str,
    df: pd.DataFrame,
    title: str = "DEG Pathway Enrichment",
    base_response: dict | None = None,
) -> dict:
    raw_count = int(df.shape[0])

    pathway_df = normalize_deg_pathway_dataframe(df)

    cleaned_count = int(pathway_df.shape[0])
    dropped_count = raw_count - cleaned_count

    results = [
        format_deg_pathway_row(row)
        for _, row in pathway_df.iterrows()
    ]

    response_data = {
        "gsea_file": gsea_file_name,
        "title": title,
        "x_field": "NES",
        "y_field": "Term",
        "size_field": "neg_log10_fdr_qval",
        "summary": {
            "raw_count": raw_count,
            "cleaned_count": cleaned_count,
            "dropped_count": dropped_count,
        },
        "results": results,
    }

    if base_response:
        response_data = {
            **base_response,
            **response_data,
        }

    return response_data
