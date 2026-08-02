from __future__ import annotations

from database.models import (
    AxisFinalRecurrentSummary,
    AxisStructureRecurrentSummary,
)


AXIS_RECURRENT_FILTER_TYPE_ITEMS = "items"
AXIS_RECURRENT_FILTER_TYPE_NUMBER = "number"


def get_distinct_axis_values(
    *,
    field_name: str,
) -> list[str]:
    lookup = f"axis__{field_name}"

    return list(
        AxisStructureRecurrentSummary.objects
        .exclude(**{lookup: ""})
        .exclude(**{f"{lookup}__isnull": True})
        .order_by(lookup)
        .values_list(lookup, flat=True)
        .distinct()
    )


def get_distinct_regulation_values() -> list[str]:
    return list(
        AxisFinalRecurrentSummary.objects
        .exclude(dominant_axis_regulation="")
        .exclude(
            dominant_axis_regulation__isnull=True
        )
        .order_by("dominant_axis_regulation")
        .values_list(
            "dominant_axis_regulation",
            flat=True,
        )
        .distinct()
    )


def build_axis_recurrent_meta() -> dict:
    axis_types = get_distinct_axis_values(
        field_name="axis_type",
    )

    regulation_patterns = (
        get_distinct_regulation_values()
    )

    return {
        "table_name":
            "axis_structure_recur_summary",

        "default_filters": {},

        "default_sort": {
            "sort_field": "dataset_count",
            "sort_order": "descend",
        },

        "pattern": {
            "format":
                "miRNA|mRNA|lncRNA|circRNA",
            "wildcard":
                "*",
            "empty_segment":
                True,
            "placeholder":
                "hsa-mir-*|ESPL1|*|*",
            "examples": [
                "*|BRD7|*|*",
                "hsa-mir-*|BRD7|BAZ1A|",
                "*|*||hsa_circ_*",
            ],
        },

        "columns": [
            {
                "field_name": "axis",
                "field_label": "Axis",
            },
            {
                "field_name": "dataset_count",
                "field_label": "Datasets",
                "sortable": True,
            },
            {
                "field_name": "context_count",
                "field_label": "Contexts",
                "sortable": True,
            },
            {
                "field_name": "source_counts",
                "field_label": "Sources",
            },
            {
                "field_name": "result_coverage",
                "field_label": "Result Coverage",
            },
            {
                "field_name": "axis_final_summary",
                "field_label": "Regulation",
            },
        ],

        "fields": [
            {
                "field_name": "axis_type",
                "field_label": "Axis Type",
                "field_type":
                    AXIS_RECURRENT_FILTER_TYPE_ITEMS,
                "options": [
                    {
                        "label": value,
                        "value": value,
                    }
                    for value in axis_types
                ],
            },
            {
                "field_name": "source",
                "field_label": "Source",
                "field_type":
                    AXIS_RECURRENT_FILTER_TYPE_ITEMS,
                "options": [
                    {
                        "label": "TCGA",
                        "value": "TCGA",
                    },
                    {
                        "label": "TIMEDB",
                        "value": "TIMEDB",
                    },
                ],
            },
            {
                "field_name": "has_axis_final",
                "field_label": "Axis Final",
                "field_type":
                    AXIS_RECURRENT_FILTER_TYPE_ITEMS,
                "options": [
                    {
                        "label": "Available",
                        "value": True,
                    },
                    {
                        "label": "Unavailable",
                        "value": False,
                    },
                ],
            },
            {
                "field_name": "has_sponge",
                "field_label": "Sponge",
                "field_type":
                    AXIS_RECURRENT_FILTER_TYPE_ITEMS,
                "options": [
                    {
                        "label": "Available",
                        "value": True,
                    },
                    {
                        "label": "Unavailable",
                        "value": False,
                    },
                ],
            },
            {
                "field_name":
                    "has_both_result_context",
                "field_label":
                    "Both Results in Same Context",
                "field_type":
                    AXIS_RECURRENT_FILTER_TYPE_ITEMS,
                "options": [
                    {
                        "label": "Yes",
                        "value": True,
                    },
                    {
                        "label": "No",
                        "value": False,
                    },
                ],
            },
            {
                "field_name":
                    "dominant_axis_regulation",
                "field_label":
                    "Dominant Regulation",
                "field_type":
                    AXIS_RECURRENT_FILTER_TYPE_ITEMS,
                "options": [
                    {
                        "label": value,
                        "value": value,
                    }
                    for value in regulation_patterns
                ],
            },
            {
                "field_name":
                    "regulation_consistent",
                "field_label":
                    "Regulation Consistency",
                "field_type":
                    AXIS_RECURRENT_FILTER_TYPE_ITEMS,
                "options": [
                    {
                        "label": "Consistent",
                        "value": True,
                    },
                    {
                        "label": "Inconsistent",
                        "value": False,
                    },
                ],
            },
            {
                "field_name":
                    "min_dataset_count",
                "field_label":
                    "Minimum Dataset Count",
                "field_type":
                    AXIS_RECURRENT_FILTER_TYPE_NUMBER,
                "minimum": 1,
            },
            {
                "field_name":
                    "min_context_count",
                "field_label":
                    "Minimum Context Count",
                "field_type":
                    AXIS_RECURRENT_FILTER_TYPE_NUMBER,
                "minimum": 1,
            },
        ],
    }
