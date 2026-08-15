from __future__ import annotations

from rest_framework import serializers

from database.models import (
    AxisDatasetSource,
    AxisFinalRecurrentSummary,
    AxisStructureRecurrentSummary,
    CanonicalAxisType,
)


class AxisRecurrentSearchFiltersSerializer(
    serializers.Serializer
):
    axis_type = serializers.ListField(
        child=serializers.ChoiceField(
            choices=CanonicalAxisType.values,
        ),
        required=False,
        allow_empty=True,
    )

    source = serializers.ListField(
        child=serializers.ChoiceField(
            choices=AxisDatasetSource.values,
        ),
        required=False,
        allow_empty=True,
    )

    min_dataset_count = serializers.IntegerField(
        min_value=1,
        required=False,
    )

    min_context_count = serializers.IntegerField(
        min_value=1,
        required=False,
    )

    has_axis_final = serializers.BooleanField(
        required=False,
    )

    has_sponge = serializers.BooleanField(
        required=False,
    )

    has_both_result_context = serializers.BooleanField(
        required=False,
    )

    regulation_available = serializers.BooleanField(
        required=False,
    )

    dominant_axis_regulation = serializers.ListField(
        child=serializers.CharField(
            allow_blank=False,
            trim_whitespace=True,
            max_length=64,
        ),
        required=False,
        allow_empty=True,
    )

    regulation_consistent = serializers.ListField(
        child=serializers.BooleanField(),
        required=False,
        allow_empty=True,
    )


class AxisRecurrentSearchRequestSerializer(
    serializers.Serializer
):
    page = serializers.IntegerField(
        min_value=1,
        required=False,
    )

    page_size = serializers.IntegerField(
        min_value=1,
        required=False,
    )

    pattern = serializers.CharField(
        required=False,
        allow_blank=True,
        trim_whitespace=True,
        max_length=1600,
    )

    filters = AxisRecurrentSearchFiltersSerializer(
        required=False,
    )

    sort_field = serializers.ChoiceField(
        choices=(
            "axis_signature",
            "axis_type",
            "miRNA",
            "mRNA",
            "lncRNA",
            "circRNA",
            "dataset_count",
            "context_count",
            "tcga_dataset_count",
            "timedb_dataset_count",
            "sc_dataset_count",
            "st_dataset_count",
            "axis_final_context_count",
            "sponge_context_count",
            "both_result_context_count",
            "regulation_pattern_count",
            "dominant_regulation_count",
            "updated_at",
        ),
        required=False,
        allow_blank=True,
    )

    sort_order = serializers.ChoiceField(
        choices=(
            "ascend",
            "descend",
        ),
        required=False,
        allow_blank=True,
    )


class AxisFinalRecurrentSummarySerializer(
    serializers.ModelSerializer
):
    """
    Default nested representation of the Axis Final recurrent summary.
    """

    dominant_regulation_ratio = (
        serializers.SerializerMethodField()
    )

    class Meta:
        model = AxisFinalRecurrentSummary
        fields = [
            "context_count",
            "observation_count",
            "regulation_pattern_count",
            "dominant_axis_regulation",
            "dominant_regulation_count",
            "dominant_regulation_ratio",
            "regulation_consistent",
            "updated_at",
        ]

    def get_dominant_regulation_ratio(
        self,
        obj,
    ) -> float | None:
        observation_count = int(
            obj.observation_count or 0
        )

        if observation_count <= 0:
            return None

        return (
            obj.dominant_regulation_count
            / observation_count
        )


class AxisStructureRecurrentSummarySerializer(
    serializers.ModelSerializer
):
    id = serializers.IntegerField(
        source="axis_id",
        read_only=True,
    )

    axis_key = serializers.CharField(
        source="axis.axis_key",
        read_only=True,
    )

    axis_signature = serializers.CharField(
        source="axis.axis_signature",
        read_only=True,
    )

    axis_type = serializers.CharField(
        source="axis.axis_type",
        read_only=True,
    )

    miRNA = serializers.CharField(
        source="axis.miRNA",
        read_only=True,
    )

    mRNA = serializers.CharField(
        source="axis.mRNA",
        read_only=True,
    )

    lncRNA = serializers.CharField(
        source="axis.lncRNA",
        read_only=True,
    )

    circRNA = serializers.CharField(
        source="axis.circRNA",
        read_only=True,
    )

    ceRNA = serializers.SerializerMethodField()
    ceRNA_type = serializers.SerializerMethodField()

    axis_final_only_context_count = (
        serializers.SerializerMethodField()
    )
    sponge_only_context_count = (
        serializers.SerializerMethodField()
    )

    has_axis_final = serializers.SerializerMethodField()
    has_sponge = serializers.SerializerMethodField()
    has_both_result_context = serializers.SerializerMethodField()

    axis_final_summary = serializers.SerializerMethodField()

    class Meta:
        model = AxisStructureRecurrentSummary
        fields = [
            "id",
            "axis_key",
            "axis_signature",
            "axis_type",

            "miRNA",
            "mRNA",
            "lncRNA",
            "circRNA",
            "ceRNA",
            "ceRNA_type",

            "dataset_count",
            "context_count",

            "tcga_dataset_count",
            "timedb_dataset_count",
            "sc_dataset_count",
            "st_dataset_count",

            "tcga_context_count",
            "timedb_context_count",
            "sc_context_count",
            "st_context_count",

            "module2_context_count",
            "module3_context_count",

            "axis_final_context_count",
            "sponge_context_count",
            "both_result_context_count",
            "axis_final_only_context_count",
            "sponge_only_context_count",

            "has_axis_final",
            "has_sponge",
            "has_both_result_context",

            "axis_final_summary",

            "summary_version",
            "updated_at",
        ]

    def get_ceRNA(self, obj) -> str:
        return obj.axis.lncRNA or obj.axis.circRNA or ""

    def get_ceRNA_type(self, obj) -> str:
        if obj.axis.lncRNA:
            return "lncRNA"

        if obj.axis.circRNA:
            return "circRNA"

        return ""

    def get_axis_final_only_context_count(
        self,
        obj,
    ) -> int:
        return max(
            0,
            obj.axis_final_context_count
            - obj.both_result_context_count,
        )

    def get_sponge_only_context_count(
        self,
        obj,
    ) -> int:
        return max(
            0,
            obj.sponge_context_count
            - obj.both_result_context_count,
        )

    def get_has_axis_final(self, obj) -> bool:
        return obj.axis_final_context_count > 0

    def get_has_sponge(self, obj) -> bool:
        return obj.sponge_context_count > 0

    def get_has_both_result_context(self, obj) -> bool:
        return obj.both_result_context_count > 0

    def get_axis_final_summary(self, obj):
        """
        Return Axis Final summary by default for every recurrent record.

        An Axis may exist only in Sponge results, so the reverse OneToOne
        relation is optional and must be handled explicitly.
        """
        try:
            summary = (
                obj.axis.axis_final_recurrent_summary
            )
        except AxisFinalRecurrentSummary.DoesNotExist:
            return None

        return AxisFinalRecurrentSummarySerializer(
            summary,
            context=self.context,
        ).data
