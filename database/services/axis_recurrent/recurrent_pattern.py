import re

from django.db.models import Q
from rest_framework.exceptions import ValidationError


AXIS_PATTERN_FIELDS = [
    "miRNA",
    "mRNA",
    "lncRNA",
    "circRNA",
]

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
    Parse an RNA-only axis pattern:

        miRNA|mRNA|lncRNA|circRNA

    Rules:
        "*"             -> any value, including empty
        empty segment   -> field must equal ""
        plain text      -> exact match
        text with "*"   -> wildcard match

    Examples:
        *|BRD7|BAZ1A|
        hsa-mir-*|*||hsa_circ_*
        hsa-mir-99a|BRD7|BAZ1A|
    """
    pattern = normalize_pattern_value(pattern)

    if not pattern:
        return {}

    if len(pattern) > MAX_AXIS_PATTERN_LENGTH:
        raise AxisRecurrentPatternError(
            "Pattern length cannot exceed "
            f"{MAX_AXIS_PATTERN_LENGTH} characters."
        )

    # Empty segments must be preserved:
    #
    #     hsa-mir-99a|BRD7||hsa_circ_001
    #
    # means:
    #     lncRNA = ""
    parts = pattern.split("|")

    if len(parts) != AXIS_PATTERN_PART_COUNT:
        raise AxisRecurrentPatternError(
            "Invalid axis pattern. Expected exactly four parts: "
            "miRNA|mRNA|lncRNA|circRNA."
        )

    parsed_pattern = {}

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

    Examples:
        hsa-mir-* -> ^hsa\\-mir\\-.*$
        *BRD7*    -> ^.*BRD7.*$

    Other regex characters are escaped, so users cannot inject
    arbitrary regular expressions.
    """
    escaped_parts = [
        re.escape(part)
        for part in value.split("*")
    ]

    return "^" + ".*".join(escaped_parts) + "$"


def build_axis_recurrent_pattern_query(
    pattern: str,
) -> Q:
    """
    Build a Q object from:

        miRNA|mRNA|lncRNA|circRNA
    """
    parsed_pattern = parse_axis_recurrent_pattern(pattern)

    query = Q()

    for field_name, value in parsed_pattern.items():
        # "*" matches every value, including empty strings.
        if value == "*":
            continue

        # Empty segment means the RNA field must be empty.
        if value == "":
            query &= Q(**{
                field_name: "",
            })
            continue

        # Internal wildcard.
        if "*" in value:
            query &= Q(**{
                f"{field_name}__iregex":
                    wildcard_value_to_regex(value),
            })
            continue

        # Exact matching. Use iexact if RNA names should be
        # matched case-insensitively.
        query &= Q(**{
            f"{field_name}__iexact": value,
        })

    return query


def apply_axis_recurrent_pattern(
    queryset,
    pattern: str,
):
    pattern = normalize_pattern_value(pattern)

    if not pattern:
        return queryset

    try:
        pattern_query = (
            build_axis_recurrent_pattern_query(pattern)
        )
    except AxisRecurrentPatternError as exc:
        raise ValidationError({
            "pattern": str(exc),
        }) from exc

    return queryset.filter(pattern_query)
