from rest_framework import serializers

from database.models import AxisRecurrentSummary


class AxisRecurrentSummarySerializer(
    serializers.ModelSerializer
):
    class Meta:
        model = AxisRecurrentSummary
        fields = [
            "id",
            "axis_signature",

            "axis_type",
            "miRNA",
            "mRNA",
            "lncRNA",
            "circRNA",

            "project_count",
            "dataset_count",
            "tcga_project_count",
            "timedb_project_count",

            "regulation_pattern_count",
            "dominant_axis_regulation",
            "dominant_regulation_count",
            "regulation_consistent",

            "updated_at",
        ]
