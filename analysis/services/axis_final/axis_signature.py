import pandas as pd


AXIS_SIGNATURE_FIELDS = [
    "axis_type",
    "miRNA",
    "mRNA",
    "lncRNA",
    "circRNA",
]


EMPTY_AXIS_VALUES = {
    "",
    "nan",
    "none",
    "null",
    "<na>",
    "na",
    "n/a",
}


AXIS_SIGNATURE_PART_COUNT = 5


def normalize_axis_signature_value(value) -> str:
    """
    Normalize one axis component before building axis_signature.

    Rules:
    - None / NaN / pd.NA -> ""
    - trim whitespace
    - common empty-like strings -> ""
    - keep original case for RNA names
    """
    if value is None:
        return ""

    try:
        if pd.isna(value):
            return ""
    except TypeError:
        pass

    value = str(value).strip()

    if value.lower() in EMPTY_AXIS_VALUES:
        return ""

    return value


def build_axis_signature(
    *,
    axis_type,
    miRNA,
    mRNA,
    lncRNA="",
    circRNA="",
) -> str:
    """
    Build stable axis signature for cross-source matching.

    Format:
        axis_type|miRNA|mRNA|lncRNA|circRNA

    Examples:
        miRNA-mRNA-lncRNA|hsa-miR-210|ESPL1|CDCA3|
        miRNA-mRNA-circRNA|hsa-miR-210|ESPL1||hsa_circ_000001
    """
    return "|".join(
        [
            normalize_axis_signature_value(axis_type),
            normalize_axis_signature_value(miRNA),
            normalize_axis_signature_value(mRNA),
            normalize_axis_signature_value(lncRNA),
            normalize_axis_signature_value(circRNA),
        ]
    )


def build_axis_signature_from_record(record: dict) -> str:
    """
    Build axis_signature from serialized axis_final record.
    """
    return build_axis_signature(
        axis_type=record.get("axis_type"),
        miRNA=record.get("miRNA"),
        mRNA=record.get("mRNA"),
        lncRNA=record.get("lncRNA"),
        circRNA=record.get("circRNA"),
    )


def build_axis_signature_from_series(row: pd.Series) -> str:
    """
    Build axis_signature from pandas Series row.
    """
    return build_axis_signature(
        axis_type=row.get("axis_type"),
        miRNA=row.get("miRNA"),
        mRNA=row.get("mRNA"),
        lncRNA=row.get("lncRNA"),
        circRNA=row.get("circRNA"),
    )


def add_axis_signature_to_dataframe(
    df: pd.DataFrame,
    signature_column: str = "axis_signature",
) -> pd.DataFrame:
    """
    Return a copy of df with axis_signature column appended.

    This should be used by reference import scripts before writing
    DatasetAxisFinalOccurrence records.
    """
    result_df = df.copy()

    for field in AXIS_SIGNATURE_FIELDS:
        if field not in result_df.columns:
            result_df[field] = pd.NA

    result_df[signature_column] = result_df.apply(
        build_axis_signature_from_series,
        axis=1,
    )

    return result_df


def add_axis_signature_to_records(
    records: list[dict],
    signature_key: str = "axis_signature",
) -> list[dict]:
    """
    Return copied records with axis_signature appended.

    This should be used when user axis_final is read from output file
    and enriched with dataset project matches.
    """
    enriched_records = []

    for record in records:
        enriched_records.append(
            {
                **record,
                signature_key: build_axis_signature_from_record(record),
            }
        )

    return enriched_records


def parse_axis_signature(axis_signature: str) -> dict:
    """
    Parse:
        axis_type|miRNA|mRNA|lncRNA|circRNA

    Empty lncRNA/circRNA parts are preserved.
    """
    axis_signature = str(axis_signature or "").strip()

    parts = axis_signature.split("|", AXIS_SIGNATURE_PART_COUNT - 1)

    if len(parts) < AXIS_SIGNATURE_PART_COUNT:
        parts.extend(
            [""] * (AXIS_SIGNATURE_PART_COUNT - len(parts))
        )

    axis_type, mirna, mrna, lncrna, circrna = parts

    return {
        "axis_type": axis_type.strip(),
        "miRNA": mirna.strip(),
        "mRNA": mrna.strip(),
        "lncRNA": lncrna.strip(),
        "circRNA": circrna.strip(),
    }
