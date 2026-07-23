from rest_framework import serializers

from database.models import (
    AxisRecurrentSummary,
    AxisSignatureProjectIndex,
)


class AxisRecurrentDetailSummarySerializer(
    serializers.ModelSerializer
):
    class Meta:
        model = AxisRecurrentSummary
        fields = [
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

            "created_at",
            "updated_at",
        ]


class AxisRecurrentProjectRecordSerializer(
    serializers.ModelSerializer
):
    row_index = serializers.IntegerField(
        source="occurrence.row_index",
        read_only=True,
        allow_null=True,
    )

    miRNA = serializers.CharField(
        source="occurrence.miRNA",
        read_only=True,
        allow_null=True,
    )

    mRNA = serializers.CharField(
        source="occurrence.mRNA",
        read_only=True,
        allow_null=True,
    )

    lncRNA = serializers.CharField(
        source="occurrence.lncRNA",
        read_only=True,
        allow_null=True,
    )

    circRNA = serializers.CharField(
        source="occurrence.circRNA",
        read_only=True,
        allow_null=True,
    )

    mRNA_log2FC = serializers.FloatField(
        source="occurrence.mRNA_log2FC",
        read_only=True,
        allow_null=True,
    )

    mRNA_regulation = serializers.CharField(
        source="occurrence.mRNA_regulation",
        read_only=True,
        allow_null=True,
    )

    miRNA_log2FC = serializers.FloatField(
        source="occurrence.miRNA_log2FC",
        read_only=True,
        allow_null=True,
    )

    miRNA_regulation = serializers.CharField(
        source="occurrence.miRNA_regulation",
        read_only=True,
        allow_null=True,
    )

    lncRNA_log2FC = serializers.FloatField(
        source="occurrence.lncRNA_log2FC",
        read_only=True,
        allow_null=True,
    )

    lncRNA_regulation = serializers.CharField(
        source="occurrence.lncRNA_regulation",
        read_only=True,
        allow_null=True,
    )

    circRNA_log2FC = serializers.FloatField(
        source="occurrence.circRNA_log2FC",
        read_only=True,
        allow_null=True,
    )

    circRNA_regulation = serializers.CharField(
        source="occurrence.circRNA_regulation",
        read_only=True,
        allow_null=True,
    )

    project_is_active = serializers.BooleanField(
        source="project.is_active",
        read_only=True,
    )

    project_row_count = serializers.IntegerField(
        source="project.row_count",
        read_only=True,
    )

    annotation_dir_name = serializers.CharField(
        source="project.annotation_dir_name",
        read_only=True,
    )

    annotation_file_prefix = serializers.CharField(
        source="project.annotation_file_prefix",
        read_only=True,
    )

    axis_final_file_name = serializers.CharField(
        source="project.axis_final_file_name",
        read_only=True,
    )

    imported_at = serializers.DateTimeField(
        source="project.imported_at",
        read_only=True,
    )

    class Meta:
        model = AxisSignatureProjectIndex
        fields = [
            "id",
            "project_id",
            "occurrence_id",

            "source",
            "module",
            "dataset_name",
            "group_type",
            "group_by",

            "axis_id",
            "axis_type",
            "axis_regulation",
            "row_index",

            "miRNA",
            "mRNA",
            "lncRNA",
            "circRNA",

            "mRNA_log2FC",
            "mRNA_regulation",

            "miRNA_log2FC",
            "miRNA_regulation",

            "lncRNA_log2FC",
            "lncRNA_regulation",

            "circRNA_log2FC",
            "circRNA_regulation",

            "project_is_active",
            "project_row_count",

            "annotation_dir_name",
            "annotation_file_prefix",
            "axis_final_file_name",

            "imported_at",
            "created_at",
        ]
