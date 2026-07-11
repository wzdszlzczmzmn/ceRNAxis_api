# database/utils/axis_recurrent_meta_utils.py

from database.models import AxisRecurrentSummary


AXIS_RECURRENT_FILTER_TYPE_ITEMS = "items"


def get_distinct_non_empty_values(
    *,
    field_name: str,
) -> list[str]:
    return list(
        AxisRecurrentSummary.objects
        .exclude(**{field_name: ""})
        .exclude(**{f"{field_name}__isnull": True})
        .order_by(field_name)
        .values_list(field_name, flat=True)
        .distinct()
    )


def build_axis_recurrent_meta() -> dict:
    axis_types = get_distinct_non_empty_values(
        field_name="axis_type",
    )

    regulation_patterns = get_distinct_non_empty_values(
        field_name="dominant_axis_regulation",
    )

    return {
        "table_name": "axis_recurrent_summary",
        "pattern": {
            "format": "miRNA|mRNA|lncRNA|circRNA",
            "wildcard": "*",
            "empty_segment": True,
            "placeholder": "hsa-mir-*|ESPL1|*|*",
            "examples": [
                "*|BRD7|*|*",
                "hsa-mir-*|BRD7|BAZ1A|",
                "*|*||hsa_circ_*",
            ],
        },
        "fields": [
            {
                "field_name": "axis_type",
                "field_label": "Axis Type",
                "field_type": AXIS_RECURRENT_FILTER_TYPE_ITEMS,
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
                "field_type": AXIS_RECURRENT_FILTER_TYPE_ITEMS,
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
                "field_name": "dominant_axis_regulation",
                "field_label": "Dominant Regulation",
                "field_type": AXIS_RECURRENT_FILTER_TYPE_ITEMS,
                "options": [
                    {
                        "label": value,
                        "value": value,
                    }
                    for value in regulation_patterns
                ],
            },
            {
                "field_name": "regulation_consistent",
                "field_label": "Regulation Consistency",
                "field_type": AXIS_RECURRENT_FILTER_TYPE_ITEMS,
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
        ],
    }
