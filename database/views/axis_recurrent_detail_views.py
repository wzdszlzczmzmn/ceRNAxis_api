from __future__ import annotations

from collections import Counter

from django.core.exceptions import ObjectDoesNotExist
from django.db.models import Q

from rest_framework.exceptions import (
    NotFound,
    ValidationError,
)
from rest_framework.response import Response
from rest_framework.views import APIView

from database.models import (
    AxisContextPresence,
    AxisDatasetSource,
    AxisFinalEvidence,
    AxisGroupType,
    AxisModule,
    AxisObservation,
    AxisResultKind,
    AxisStructureRecurrentSummary,
)
from database.serializers.axis_recurrent_detail_serializers import (
    AxisRecurrentContextRecordSerializer,
    AxisRecurrentDetailSummarySerializer,
)


class AxisRecurrentDetailView(APIView):
    """
    Return one recurrent Axis summary and its active context-level records.

    Request:
        GET /api/database/axis_recurrent_detail/
            ?signature=<axis_signature>

    Optional context filters:
        source=TCGA
        module=module2
        dataset_name=TCGA_BRCA_mRNA
        group_type=none
        group_by=
        group_value=
        has_axis_final=true
        has_sponge=true
        has_both_result_context=true
        axis_regulation=up_down_up

    The endpoint intentionally represents active summary data only.
    Historical inactive artifacts should be exposed by a separate history API.
    """

    MAX_SIGNATURE_LENGTH = 1600

    ALLOWED_SOURCES = set(
        AxisDatasetSource.values
    )
    ALLOWED_MODULES = set(
        AxisModule.values
    )
    ALLOWED_GROUP_TYPES = set(
        AxisGroupType.values
    )

    @staticmethod
    def normalize_query_value(
        request,
        field_name: str,
    ) -> str:
        return str(
            request.query_params.get(field_name) or ""
        ).strip()

    def get_signature(self, request) -> str:
        signature = self.normalize_query_value(
            request,
            "signature",
        )

        if not signature:
            raise ValidationError({
                "signature": (
                    "Missing required parameter: signature."
                ),
            })

        if len(signature) > self.MAX_SIGNATURE_LENGTH:
            raise ValidationError({
                "signature": (
                    "signature cannot exceed "
                    f"{self.MAX_SIGNATURE_LENGTH} characters."
                ),
            })

        return signature

    @staticmethod
    def parse_boolean(
        value,
        *,
        field_name: str,
    ) -> bool | None:
        if value in {None, ""}:
            return None

        if isinstance(value, bool):
            return value

        normalized = str(value).strip().lower()

        if normalized in {
            "1",
            "true",
            "yes",
        }:
            return True

        if normalized in {
            "0",
            "false",
            "no",
        }:
            return False

        raise ValidationError({
            field_name: (
                f"{field_name} must be a boolean."
            ),
        })

    @staticmethod
    def get_summary(
        *,
        signature: str,
    ) -> AxisStructureRecurrentSummary:
        try:
            return (
                AxisStructureRecurrentSummary.objects
                .select_related(
                    "axis",
                    "axis__axis_final_recurrent_summary",
                )
                .get(
                    axis__axis_signature=signature,
                )
            )
        except AxisStructureRecurrentSummary.DoesNotExist as exc:
            raise NotFound(
                "Recurrent Axis was not found."
            ) from exc
        except (
            AxisStructureRecurrentSummary
            .MultipleObjectsReturned
        ) as exc:
            raise RuntimeError(
                "Multiple recurrent summaries were found "
                "for the same axis signature."
            ) from exc

    def apply_filters(
        self,
        queryset,
        *,
        request,
        axis_id: int,
    ):
        source = self.normalize_query_value(
            request,
            "source",
        ).upper()

        module = self.normalize_query_value(
            request,
            "module",
        ).lower()

        dataset_name = self.normalize_query_value(
            request,
            "dataset_name",
        )

        group_type = self.normalize_query_value(
            request,
            "group_type",
        ).lower()

        group_by = self.normalize_query_value(
            request,
            "group_by",
        )

        group_value = self.normalize_query_value(
            request,
            "group_value",
        )

        axis_regulation = self.normalize_query_value(
            request,
            "axis_regulation",
        )

        has_axis_final = self.parse_boolean(
            request.query_params.get(
                "has_axis_final"
            ),
            field_name="has_axis_final",
        )

        has_sponge = self.parse_boolean(
            request.query_params.get(
                "has_sponge"
            ),
            field_name="has_sponge",
        )

        has_both_result_context = self.parse_boolean(
            request.query_params.get(
                "has_both_result_context"
            ),
            field_name="has_both_result_context",
        )

        if source:
            if source not in self.ALLOWED_SOURCES:
                raise ValidationError({
                    "source": (
                        "Unsupported source. Allowed values: "
                        + ", ".join(
                            sorted(self.ALLOWED_SOURCES)
                        )
                        + "."
                    ),
                })

            queryset = queryset.filter(
                context__dataset_source=source,
            )

        if module:
            if module not in self.ALLOWED_MODULES:
                raise ValidationError({
                    "module": (
                        "Unsupported module. Allowed values: "
                        + ", ".join(
                            sorted(self.ALLOWED_MODULES)
                        )
                        + "."
                    ),
                })

            queryset = queryset.filter(
                context__module=module,
            )

        if dataset_name:
            queryset = queryset.filter(
                context__dataset_metadata_id=dataset_name,
            )

        if group_type:
            if group_type not in self.ALLOWED_GROUP_TYPES:
                raise ValidationError({
                    "group_type": (
                        "Unsupported group_type. Allowed values: "
                        + ", ".join(
                            sorted(self.ALLOWED_GROUP_TYPES)
                        )
                        + "."
                    ),
                })

            queryset = queryset.filter(
                context__group_type=group_type,
            )

        if "group_by" in request.query_params:
            queryset = queryset.filter(
                context__group_by=group_by,
            )

        if "group_value" in request.query_params:
            queryset = queryset.filter(
                context__group_value=group_value,
            )

        if has_axis_final is not None:
            queryset = queryset.filter(
                has_axis_final=has_axis_final,
            )

        if has_sponge is not None:
            queryset = queryset.filter(
                has_sponge=has_sponge,
            )

        if has_both_result_context is True:
            queryset = queryset.filter(
                has_axis_final=True,
                has_sponge=True,
            )
        elif has_both_result_context is False:
            queryset = queryset.filter(
                Q(has_axis_final=False)
                | Q(has_sponge=False)
            )

        if axis_regulation:
            queryset = queryset.filter(
                context__result_artifacts__is_active=True,
                context__result_artifacts__result_kind=(
                    AxisResultKind.AXIS_FINAL
                ),
                **{
                    (
                        "context__result_artifacts__"
                        "observations__axis_id"
                    ): axis_id,
                    (
                        "context__result_artifacts__"
                        "observations__axis_final_evidence__"
                        "axis_regulation__iexact"
                    ): axis_regulation,
                },
            )

        return queryset.distinct()

    @staticmethod
    def load_observations(
        *,
        axis_id: int,
        context_ids: list[int],
    ) -> tuple[
        dict[int, dict[str, AxisObservation | None]],
        list[AxisObservation],
    ]:
        observations_by_context = {
            context_id: {
                AxisResultKind.AXIS_FINAL: None,
                AxisResultKind.SPONGE: None,
            }
            for context_id in context_ids
        }

        if not context_ids:
            return observations_by_context, []

        observations = list(
            AxisObservation.objects
            .filter(
                axis_id=axis_id,
                artifact__context_id__in=context_ids,
                artifact__is_active=True,
                artifact__context__is_active=True,
            )
            .select_related(
                "artifact",
                "axis_final_evidence",
                "sponge_evidence",
            )
            .order_by(
                "artifact__context_id",
                "artifact__result_kind",
                "id",
            )
        )

        for observation in observations:
            context_id = observation.artifact.context_id
            result_kind = observation.artifact.result_kind

            if result_kind not in {
                AxisResultKind.AXIS_FINAL,
                AxisResultKind.SPONGE,
            }:
                continue

            current = observations_by_context[
                context_id
            ][result_kind]

            if current is not None:
                raise RuntimeError(
                    "Multiple active observations were found "
                    "for the same axis, context and result kind: "
                    f"axis={axis_id}, "
                    f"context={context_id}, "
                    f"result_kind={result_kind}."
                )

            observations_by_context[
                context_id
            ][result_kind] = observation

        return observations_by_context, observations

    @staticmethod
    def build_statistics(
        *,
        records: list[AxisContextPresence],
        observations: list[AxisObservation],
    ) -> dict:
        source_counter = Counter()
        module_counter = Counter()
        group_type_counter = Counter()
        regulation_counter = Counter()

        dataset_names = set()

        for record in records:
            context = record.context

            dataset_names.add(
                context.dataset_name
            )
            source_counter[
                context.dataset_source
            ] += 1
            module_counter[
                context.module
            ] += 1
            group_type_counter[
                context.group_type
            ] += 1

        for observation in observations:
            if (
                observation.artifact.result_kind
                != AxisResultKind.AXIS_FINAL
            ):
                continue

            try:
                regulation = (
                    observation
                    .axis_final_evidence
                    .axis_regulation
                )
            except (
                AxisFinalEvidence.DoesNotExist,
                ObjectDoesNotExist,
            ):
                continue

            regulation = str(
                regulation or ""
            ).strip()

            if regulation:
                regulation_counter[regulation] += 1

        return {
            "context_count": len(records),
            "dataset_count": len(dataset_names),

            "observation_count": sum(
                record.observation_count
                for record in records
            ),
            "axis_final_observation_count": sum(
                record.axis_final_observation_count
                for record in records
            ),
            "sponge_observation_count": sum(
                record.sponge_observation_count
                for record in records
            ),

            "axis_final_context_count": sum(
                1
                for record in records
                if record.has_axis_final
            ),
            "sponge_context_count": sum(
                1
                for record in records
                if record.has_sponge
            ),
            "both_result_context_count": sum(
                1
                for record in records
                if (
                    record.has_axis_final
                    and record.has_sponge
                )
            ),

            "source_counts": dict(
                sorted(source_counter.items())
            ),
            "module_counts": dict(
                sorted(module_counter.items())
            ),
            "group_type_counts": dict(
                sorted(group_type_counter.items())
            ),
            "regulation_counts": dict(
                sorted(
                    regulation_counter.items(),
                    key=lambda item: (
                        -item[1],
                        item[0],
                    ),
                )
            ),
        }

    def get(self, request):
        signature = self.get_signature(request)

        summary = self.get_summary(
            signature=signature,
        )
        axis = summary.axis

        queryset = (
            AxisContextPresence.objects
            .filter(
                axis_id=axis.id,
                context__is_active=True,
            )
            .select_related(
                "context",
                "context__dataset_metadata",
            )
        )

        queryset = self.apply_filters(
            queryset,
            request=request,
            axis_id=axis.id,
        )

        records = list(
            queryset.order_by(
                "context__dataset_source",
                "context__module",
                "context__dataset_metadata_id",
                "context__group_type",
                "context__group_by",
                "context__group_value",
                "context_id",
            )
        )

        context_ids = [
            record.context_id
            for record in records
        ]

        (
            observations_by_context,
            observations,
        ) = self.load_observations(
            axis_id=axis.id,
            context_ids=context_ids,
        )

        serializer_context = {
            "request": request,
            "observations_by_context":
                observations_by_context,
        }

        return Response({
            "success": True,
            "summary": (
                AxisRecurrentDetailSummarySerializer(
                    summary,
                    context=serializer_context,
                ).data
            ),
            "statistics": self.build_statistics(
                records=records,
                observations=observations,
            ),
            "count": len(records),
            "results": (
                AxisRecurrentContextRecordSerializer(
                    records,
                    many=True,
                    context=serializer_context,
                ).data
            ),
        })
