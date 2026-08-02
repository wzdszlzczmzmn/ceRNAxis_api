from __future__ import annotations

import hashlib
import json
import math
import unicodedata
from pathlib import Path
from typing import Any

import pandas as pd

from database.models import CanonicalAxisType


EMPTY_TEXT_VALUES = {
    "",
    "nan",
    "none",
    "null",
    "<na>",
    "na",
    "n/a",
}


class AxisImportValidationError(ValueError):
    pass


def normalize_text(value: Any) -> str:
    if value is None:
        return ""

    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass

    value = unicodedata.normalize(
        "NFKC",
        str(value),
    ).strip()

    if value.lower() in EMPTY_TEXT_VALUES:
        return ""

    return value


def normalize_optional_float(
    value: Any,
    *,
    field_name: str,
    row_index: int,
) -> float | None:
    if value is None:
        return None

    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass

    text = str(value).strip()

    if text.lower() in EMPTY_TEXT_VALUES:
        return None

    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise AxisImportValidationError(
            f"Row {row_index}: "
            f"field '{field_name}' must be numeric or empty."
        ) from exc

    if not math.isfinite(number):
        raise AxisImportValidationError(
            f"Row {row_index}: "
            f"field '{field_name}' must be finite."
        )

    return number


def derive_canonical_axis_type(
    *,
    miRNA: str,
    mRNA: str,
    lncRNA: str,
    circRNA: str,
    row_index: int,
) -> str:
    if not miRNA:
        raise AxisImportValidationError(
            f"Row {row_index}: miRNA cannot be empty."
        )

    if not mRNA:
        raise AxisImportValidationError(
            f"Row {row_index}: mRNA cannot be empty."
        )

    if lncRNA and circRNA:
        raise AxisImportValidationError(
            f"Row {row_index}: "
            "lncRNA and circRNA cannot both be non-empty."
        )

    if lncRNA:
        return CanonicalAxisType.MRNA_MIRNA_LNCRNA

    if circRNA:
        return CanonicalAxisType.MRNA_MIRNA_CIRCRNA

    return CanonicalAxisType.MRNA_MIRNA


def build_axis_signature(
    *,
    miRNA: str,
    mRNA: str,
    lncRNA: str = "",
    circRNA: str = "",
) -> str:
    return "|".join([
        miRNA,
        mRNA,
        lncRNA,
        circRNA,
    ])


def build_axis_key(axis_signature: str) -> str:
    return hashlib.sha256(
        axis_signature.encode("utf-8")
    ).hexdigest()


def calculate_file_sha256(file_path: Path) -> str:
    hasher = hashlib.sha256()

    with Path(file_path).open("rb") as handle:
        for chunk in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            hasher.update(chunk)

    return hasher.hexdigest()


def to_json_value(value: Any) -> Any:
    if value is None:
        return None

    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass

    if hasattr(value, "item"):
        try:
            value = value.item()
        except (TypeError, ValueError):
            pass

    try:
        json.dumps(value)
    except (TypeError, ValueError):
        return str(value)

    return value
