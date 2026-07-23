from collections import Counter

from django.db.models import Q

from rest_framework.exceptions import (
    NotFound,
    ValidationError,
)
from rest_framework.response import Response
from rest_framework.views import APIView

from database.models import (
    AxisRecurrentSummary,
    AxisSignatureProjectIndex,
)
from database.serializers.axis_recurrent_detail_serializers import (
    AxisRecurrentDetailSummarySerializer,
    AxisRecurrentProjectRecordSerializer,
)


class AxisRecurrentDetailView(APIView):
    """
    Return one recurrent ceRNA axis summary and all project-level
    occurrence records.

    Request:
        GET /api/database/axis_recurrent_detail/
            ?signature=<axis_signature>

    Optional filters:
        source=TCGA
        module=module2
        dataset_name=TCGA_BRCA
        group_type=none
        group_by=
        axis_regulation=Up-Down-Up
        active_only=true
    """

    MAX_SIGNATURE_LENGTH = 1024

    ALLOWED_SOURCES = {
        "TCGA",
        "TIMEDB",
    }

    ALLOWED_MODULES = {
        "module2",
        "module3",
    }

    ALLOWED_GROUP_TYPES = {
        "none",
        "common",
        "grade",
        "stage",
    }

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
        default: bool | None = None,
    ) -> bool | None:
        if value in {None, ""}:
            return default

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
    ) -> AxisRecurrentSummary:
        try:
            return AxisRecurrentSummary.objects.get(
                axis_signature=signature,
            )
        except AxisRecurrentSummary.DoesNotExist as exc:
            raise NotFound(
                "Recurrent ceRNA axis was not found."
            ) from exc

    def apply_filters(
        self,
        queryset,
        *,
        request,
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

        axis_regulation = self.normalize_query_value(
            request,
            "axis_regulation",
        )

        active_only = self.parse_boolean(
            request.query_params.get("active_only"),
            field_name="active_only",
            default=True,
        )

        if source:
            if source not in self.ALLOWED_SOURCES:
                raise ValidationError({
                    "source": (
                        "Unsupported source. "
                        "Allowed values are TCGA and TIMEDB."
                    ),
                })

            queryset = queryset.filter(
                source=source,
            )

        if module:
            if module not in self.ALLOWED_MODULES:
                raise ValidationError({
                    "module": (
                        "Unsupported module. "
                        "Allowed values are module2 and module3."
                    ),
                })

            queryset = queryset.filter(
                module=module,
            )

        if dataset_name:
            queryset = queryset.filter(
                dataset_name=dataset_name,
            )

        if group_type:
            if group_type not in self.ALLOWED_GROUP_TYPES:
                raise ValidationError({
                    "group_type": (
                        "Unsupported group_type. "
                        "Allowed values are none, common, "
                        "grade and stage."
                    ),
                })

            queryset = queryset.filter(
                group_type=group_type,
            )

        if "group_by" in request.query_params:
            queryset = queryset.filter(
                group_by=group_by,
            )

        if axis_regulation:
            queryset = queryset.filter(
                axis_regulation=axis_regulation,
            )

        if active_only:
            queryset = queryset.filter(
                project__is_active=True,
            )

        return queryset

    @staticmethod
    def build_statistics(
        records: list[AxisSignatureProjectIndex],
    ) -> dict:
        source_counter = Counter()
        module_counter = Counter()
        group_type_counter = Counter()
        regulation_counter = Counter()

        dataset_names = set()
        project_ids = set()

        for record in records:
            project_ids.add(record.project_id)
            dataset_names.add(record.dataset_name)

            if record.source:
                source_counter[record.source] += 1

            if record.module:
                module_counter[record.module] += 1

            if record.group_type:
                group_type_counter[record.group_type] += 1

            if record.axis_regulation:
                regulation_counter[
                    record.axis_regulation
                ] += 1

        return {
            "record_count": len(records),
            "project_count": len(project_ids),
            "dataset_count": len(dataset_names),
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

        queryset = (
            AxisSignatureProjectIndex.objects
            .filter(
                axis_signature=signature,
            )
            .select_related(
                "project",
                "occurrence",
            )
        )

        queryset = self.apply_filters(
            queryset,
            request=request,
        )

        queryset = queryset.order_by(
            "source",
            "module",
            "dataset_name",
            "group_type",
            "group_by",
            "axis_id",
            "project_id",
        )

        records = list(queryset)

        return Response({
            "success": True,
            "summary": (
                AxisRecurrentDetailSummarySerializer(
                    summary
                ).data
            ),
            "statistics": self.build_statistics(
                records
            ),
            "count": len(records),
            "results": (
                AxisRecurrentProjectRecordSerializer(
                    records,
                    many=True,
                ).data
            ),
        })
