from __future__ import annotations

from django.core.exceptions import ObjectDoesNotExist
from rest_framework import serializers

from database.models import (
    AxisContextPresence,
    AxisFinalEvidence,
    AxisFinalRecurrentSummary,
    AxisObservation,
    AxisResultArtifact,
    AxisStructureRecurrentSummary,
    SpongeEvidence,
)


class AxisFinalRecurrentSummarySerializer(
    serializers.ModelSerializer
):
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
            "created_at",
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


class AxisRecurrentDetailSummarySerializer(
    serializers.ModelSerializer
):
    axis_id = serializers.IntegerField(
        source="axis.id",
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

    axis_final_summary = serializers.SerializerMethodField()

    class Meta:
        model = AxisStructureRecurrentSummary
        fields = [
            "axis_id",
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
            "tcga_context_count",
            "timedb_context_count",

            "module2_context_count",
            "module3_context_count",

            "axis_final_context_count",
            "sponge_context_count",
            "both_result_context_count",

            "axis_final_summary",

            "summary_version",
            "created_at",
            "updated_at",
        ]

    def get_ceRNA(self, obj) -> str:
        return (
            obj.axis.lncRNA
            or obj.axis.circRNA
            or ""
        )

    def get_ceRNA_type(self, obj) -> str:
        if obj.axis.lncRNA:
            return "lncRNA"

        if obj.axis.circRNA:
            return "circRNA"

        return ""

    def get_axis_final_summary(self, obj):
        try:
            summary = (
                obj.axis.axis_final_recurrent_summary
            )
        except ObjectDoesNotExist:
            return None

        return AxisFinalRecurrentSummarySerializer(
            summary,
            context=self.context,
        ).data


class AxisResultArtifactSerializer(
    serializers.ModelSerializer
):
    class Meta:
        model = AxisResultArtifact
        fields = [
            "id",
            "result_kind",
            "file_name",
            "file_sha256",
            "row_count",
            "schema_version",
            "is_active",
            "imported_at",
            "updated_at",
        ]


class AxisFinalEvidenceSerializer(
    serializers.ModelSerializer
):
    class Meta:
        model = AxisFinalEvidence
        fields = [
            "axis_regulation",

            "mRNA_log2FC",
            "mRNA_regulation",

            "miRNA_log2FC",
            "miRNA_regulation",

            "lncRNA_log2FC",
            "lncRNA_regulation",

            "circRNA_log2FC",
            "circRNA_regulation",
        ]


class SpongeEvidenceSerializer(
    serializers.ModelSerializer
):
    class Meta:
        model = SpongeEvidence
        fields = [
            "cor",
            "pcor",
            "mscor",
        ]


class AxisFinalObservationSerializer(
    serializers.ModelSerializer
):
    artifact = AxisResultArtifactSerializer(
        read_only=True,
    )

    evidence = serializers.SerializerMethodField()

    class Meta:
        model = AxisObservation
        fields = [
            "id",
            "row_index",
            "source_axis_id",
            "source_axis_type",
            "extra_data",
            "artifact",
            "evidence",
            "created_at",
        ]

    def get_evidence(self, obj):
        try:
            evidence = obj.axis_final_evidence
        except ObjectDoesNotExist:
            return None

        return AxisFinalEvidenceSerializer(
            evidence,
            context=self.context,
        ).data


class SpongeObservationSerializer(
    serializers.ModelSerializer
):
    artifact = AxisResultArtifactSerializer(
        read_only=True,
    )

    evidence = serializers.SerializerMethodField()

    class Meta:
        model = AxisObservation
        fields = [
            "id",
            "row_index",
            "source_axis_id",
            "source_axis_type",
            "extra_data",
            "artifact",
            "evidence",
            "created_at",
        ]

    def get_evidence(self, obj):
        try:
            evidence = obj.sponge_evidence
        except ObjectDoesNotExist:
            return None

        return SpongeEvidenceSerializer(
            evidence,
            context=self.context,
        ).data


class AxisRecurrentContextRecordSerializer(
    serializers.ModelSerializer
):
    context_id = serializers.IntegerField(
        source="context.id",
        read_only=True,
    )

    dataset_source = serializers.CharField(
        source="context.dataset_source",
        read_only=True,
    )

    module = serializers.CharField(
        source="context.module",
        read_only=True,
    )

    dataset_name = serializers.CharField(
        source="context.dataset_name",
        read_only=True,
    )

    group_type = serializers.CharField(
        source="context.group_type",
        read_only=True,
    )

    group_by = serializers.CharField(
        source="context.group_by",
        read_only=True,
    )

    annotation_dir_name = serializers.CharField(
        source="context.annotation_dir_name",
        read_only=True,
    )

    annotation_file_prefix = serializers.CharField(
        source="context.annotation_file_prefix",
        read_only=True,
    )

    context_is_active = serializers.BooleanField(
        source="context.is_active",
        read_only=True,
    )

    has_both_results = serializers.SerializerMethodField()

    axis_final = serializers.SerializerMethodField()
    sponge = serializers.SerializerMethodField()

    class Meta:
        model = AxisContextPresence
        fields = [
            "id",
            "context_id",

            "dataset_source",
            "module",
            "dataset_name",
            "group_type",
            "group_by",

            "annotation_dir_name",
            "annotation_file_prefix",
            "context_is_active",

            "observation_count",
            "axis_final_observation_count",
            "sponge_observation_count",

            "has_axis_final",
            "has_sponge",
            "has_both_results",

            "axis_final",
            "sponge",

            "created_at",
            "updated_at",
        ]

    def get_has_both_results(self, obj) -> bool:
        return bool(
            obj.has_axis_final
            and obj.has_sponge
        )

    def get_axis_final(self, obj):
        observations_by_context = self.context.get(
            "observations_by_context",
            {},
        )

        observation = (
            observations_by_context
            .get(obj.context_id, {})
            .get("axis_final")
        )

        if observation is None:
            return None

        return AxisFinalObservationSerializer(
            observation,
            context=self.context,
        ).data

    def get_sponge(self, obj):
        observations_by_context = self.context.get(
            "observations_by_context",
            {},
        )

        observation = (
            observations_by_context
            .get(obj.context_id, {})
            .get("sponge")
        )

        if observation is None:
            return None

        return SpongeObservationSerializer(
            observation,
            context=self.context,
        ).data
