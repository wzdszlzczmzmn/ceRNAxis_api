from __future__ import annotations

import re

from django.db.models import Q
from rest_framework.exceptions import ValidationError


AXIS_PATTERN_FIELDS = (
    "miRNA",
    "mRNA",
    "lncRNA",
    "circRNA",
)

AXIS_PATTERN_PART_COUNT = len(AXIS_PATTERN_FIELDS)

MAX_AXIS_PATTERN_LENGTH = 1600
MAX_AXIS_PATTERN_PART_LENGTH = 512


class AxisRecurrentPatternError(ValueError):
    pass


def normalize_pattern_value(value) -> str:
    if value is None:
        return ""

    return str(value).strip()


def parse_axis_recurrent_pattern(
    pattern: str,
) -> dict[str, str]:
    """
    Parse:
        miRNA|mRNA|lncRNA|circRNA

    Rules:
        "*"           -> any value, including empty
        empty segment -> field must equal ""
        plain text    -> case-insensitive exact match
        text with "*" -> case-insensitive wildcard match
    """
    pattern = normalize_pattern_value(pattern)

    if not pattern:
        return {}

    if len(pattern) > MAX_AXIS_PATTERN_LENGTH:
        raise AxisRecurrentPatternError(
            "Pattern length cannot exceed "
            f"{MAX_AXIS_PATTERN_LENGTH} characters."
        )

    parts = pattern.split("|")

    if len(parts) != AXIS_PATTERN_PART_COUNT:
        raise AxisRecurrentPatternError(
            "Invalid axis pattern. Expected exactly four parts: "
            "miRNA|mRNA|lncRNA|circRNA."
        )

    parsed_pattern: dict[str, str] = {}

    for field_name, raw_value in zip(
        AXIS_PATTERN_FIELDS,
        parts,
    ):
        value = normalize_pattern_value(raw_value)

        if len(value) > MAX_AXIS_PATTERN_PART_LENGTH:
            raise AxisRecurrentPatternError(
                f"Pattern part '{field_name}' cannot exceed "
                f"{MAX_AXIS_PATTERN_PART_LENGTH} characters."
            )

        parsed_pattern[field_name] = value

    return parsed_pattern


def wildcard_value_to_regex(value: str) -> str:
    """
    Convert '*' wildcard syntax into an anchored safe regex.

    Regex metacharacters supplied by users are escaped.
    """
    escaped_parts = [
        re.escape(part)
        for part in value.split("*")
    ]

    return "^" + ".*".join(escaped_parts) + "$"


def build_axis_recurrent_pattern_query(
    pattern: str,
    *,
    field_prefix: str = "axis__",
) -> Q:
    """
    Build a Q object for a queryset whose RNA fields are located on
    CanonicalAxis.

    AxisStructureRecurrentSummary uses the default prefix ``axis__``.
    Passing ``field_prefix=""`` remains useful for direct CanonicalAxis
    querysets.
    """
    parsed_pattern = parse_axis_recurrent_pattern(pattern)

    query = Q()

    for field_name, value in parsed_pattern.items():
        lookup_base = f"{field_prefix}{field_name}"

        if value == "*":
            continue

        if value == "":
            query &= Q(**{
                lookup_base: "",
            })
            continue

        if "*" in value:
            query &= Q(**{
                f"{lookup_base}__iregex":
                    wildcard_value_to_regex(value),
            })
            continue

        query &= Q(**{
            f"{lookup_base}__iexact": value,
        })

    return query


def apply_axis_recurrent_pattern(
    queryset,
    pattern: str,
    *,
    field_prefix: str = "axis__",
):
    pattern = normalize_pattern_value(pattern)

    if not pattern:
        return queryset

    try:
        pattern_query = build_axis_recurrent_pattern_query(
            pattern,
            field_prefix=field_prefix,
        )
    except AxisRecurrentPatternError as exc:
        raise ValidationError({
            "pattern": str(exc),
        }) from exc

    return queryset.filter(pattern_query)
