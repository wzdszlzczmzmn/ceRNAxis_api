from django.db.models import Q

from rest_framework.generics import ListAPIView
from rest_framework.exceptions import ValidationError
from rest_framework.response import Response
from rest_framework.views import APIView

from database.models import AxisRecurrentSummary
from database.pagination.standard_pagination import (
    StandardPageNumberPagination,
)
from database.serializers.axis_recurrent_serializers import (
    AxisRecurrentSummarySerializer,
)
from database.services.axis_recurrent import (
    apply_axis_recurrent_pattern,
)
from database.utils.axis_recurrent_meta_utils import build_axis_recurrent_meta


class AxisRecurrentSummarySearchView(ListAPIView):
    """
    Search and paginate recurrent ceRNA axis summaries.

    Pattern format:
        miRNA|mRNA|lncRNA|circRNA

    Pattern rules:
        *    -> any value
        ||   -> the corresponding RNA field must be empty

    Example:
        {
            "page": 1,
            "page_size": 20,
            "pattern": "hsa-mir-*|BRD7|*|",
            "filters": {
                "axis_type": [
                    "mRNA-miRNA-lncRNA"
                ],
                "source": [
                    "TCGA"
                ],
                "min_project_count": 2
            },
            "sort_field": "project_count",
            "sort_order": "descend"
        }
    """

    serializer_class = AxisRecurrentSummarySerializer
    pagination_class = StandardPageNumberPagination

    SORT_FIELD_MAP = {
        "axis_signature": "axis_signature",
        "axis_type": "axis_type",

        "miRNA": "miRNA",
        "mRNA": "mRNA",
        "lncRNA": "lncRNA",
        "circRNA": "circRNA",

        "project_count": "project_count",
        "dataset_count": "dataset_count",
        "tcga_project_count": "tcga_project_count",
        "timedb_project_count": "timedb_project_count",

        "regulation_pattern_count": (
            "regulation_pattern_count"
        ),
        "dominant_regulation_count": (
            "dominant_regulation_count"
        ),
        "updated_at": "updated_at",
    }

    SOURCE_VALUES = {
        "TCGA",
        "TIMEDB",
    }

    def post(self, request, *args, **kwargs):
        return self.list(request, *args, **kwargs)

    def get_request_data(self) -> dict:
        data = self.request.data

        if not isinstance(data, dict):
            raise ValidationError({
                "detail": "Request body must be a JSON object.",
            })

        return data

    @staticmethod
    def parse_non_negative_integer(
        value,
        *,
        field_name: str,
    ) -> int | None:
        if value in {None, ""}:
            return None

        try:
            value = int(value)
        except (TypeError, ValueError) as exc:
            raise ValidationError({
                field_name: (
                    f"{field_name} must be an integer."
                ),
            }) from exc

        if value < 0:
            raise ValidationError({
                field_name: (
                    f"{field_name} cannot be negative."
                ),
            })

        return value

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

        if normalized in {"1", "true", "yes"}:
            return True

        if normalized in {"0", "false", "no"}:
            return False

        raise ValidationError({
            field_name: (
                f"{field_name} must be a boolean."
            ),
        })

    def apply_source_filter(
        self,
        queryset,
        source_values,
    ):
        if not source_values:
            return queryset

        if not isinstance(source_values, list):
            raise ValidationError({
                "filters.source": (
                    "source must be an array."
                ),
            })

        normalized_sources = {
            str(value).strip().upper()
            for value in source_values
            if str(value).strip()
        }

        unsupported_sources = (
            normalized_sources - self.SOURCE_VALUES
        )

        if unsupported_sources:
            raise ValidationError({
                "filters.source": (
                    "Unsupported source value(s): "
                    f"{', '.join(sorted(unsupported_sources))}."
                ),
            })

        source_query = Q()

        if "TCGA" in normalized_sources:
            source_query |= Q(tcga_project_count__gt=0)

        if "TIMEDB" in normalized_sources:
            source_query |= Q(timedb_project_count__gt=0)

        if source_query:
            queryset = queryset.filter(source_query)

        return queryset

    def apply_filters(
            self,
            queryset,
            filters: dict,
    ):
        if not isinstance(filters, dict):
            raise ValidationError({
                "filters": "filters must be an object.",
            })

        supported_fields = {
            "axis_type",
            "source",
            "dominant_axis_regulation",
            "regulation_consistent",
        }

        unsupported_fields = set(filters) - supported_fields

        if unsupported_fields:
            raise ValidationError({
                "filters": (
                    "Unsupported filter field(s): "
                    f"{', '.join(sorted(unsupported_fields))}."
                ),
            })

        axis_types = filters.get("axis_type") or []

        if axis_types:
            queryset = queryset.filter(
                axis_type__in=axis_types,
            )

        dominant_regulations = (
                filters.get("dominant_axis_regulation") or []
        )

        if dominant_regulations:
            queryset = queryset.filter(
                dominant_axis_regulation__in=(
                    dominant_regulations
                ),
            )

        consistency_values = (
                filters.get("regulation_consistent") or []
        )

        if consistency_values:
            normalized_values = []

            for value in consistency_values:
                if isinstance(value, bool):
                    normalized_values.append(value)
                    continue

                normalized = str(value).strip().lower()

                if normalized in {"true", "1", "yes"}:
                    normalized_values.append(True)
                elif normalized in {"false", "0", "no"}:
                    normalized_values.append(False)
                else:
                    raise ValidationError({
                        "filters.regulation_consistent": (
                            f"Invalid boolean value: {value}"
                        ),
                    })

            queryset = queryset.filter(
                regulation_consistent__in=(
                    normalized_values
                ),
            )

        sources = {
            str(value).strip().upper()
            for value in (filters.get("source") or [])
            if str(value).strip()
        }

        unsupported_sources = sources - {
            "TCGA",
            "TIMEDB",
        }

        if unsupported_sources:
            raise ValidationError({
                "filters.source": (
                    "Unsupported source value(s): "
                    f"{', '.join(sorted(unsupported_sources))}."
                ),
            })

        if sources == {"TCGA"}:
            queryset = queryset.filter(
                tcga_project_count__gt=0,
            )

        elif sources == {"TIMEDB"}:
            queryset = queryset.filter(
                timedb_project_count__gt=0,
            )

        elif sources == {"TCGA", "TIMEDB"}:
            # OR semantics:
            # This normally retains all rows with either source.
            queryset = queryset.filter(
                Q(tcga_project_count__gt=0)
                | Q(timedb_project_count__gt=0)
            )

        return queryset

    def apply_sorting(
        self,
        queryset,
        *,
        sort_field,
        sort_order,
    ):
        if not sort_field:
            return queryset.order_by(
                "-project_count",
                "-dataset_count",
                "axis_signature",
            )

        sort_lookup = self.SORT_FIELD_MAP.get(
            sort_field
        )

        if sort_lookup is None:
            raise ValidationError({
                "sort_field": (
                    f"Unsupported sort field: "
                    f"{sort_field}"
                ),
            })

        if sort_order not in {
            None,
            "",
            "ascend",
            "descend",
        }:
            raise ValidationError({
                "sort_order": (
                    "sort_order must be either "
                    "'ascend' or 'descend'."
                ),
            })

        if sort_order == "descend":
            sort_lookup = f"-{sort_lookup}"

        return queryset.order_by(
            sort_lookup,
            "axis_signature",
        )

    def get_queryset(self):
        data = self.get_request_data()

        pattern = str(
            data.get("pattern") or ""
        ).strip()

        filters = data.get("filters", {}) or {}

        sort_field = data.get("sort_field")
        sort_order = data.get("sort_order")

        queryset = AxisRecurrentSummary.objects.all()

        queryset = apply_axis_recurrent_pattern(
            queryset,
            pattern,
        )

        queryset = self.apply_filters(
            queryset,
            filters,
        )

        queryset = self.apply_sorting(
            queryset,
            sort_field=sort_field,
            sort_order=sort_order,
        )

        return queryset


class AxisRecurrentMetaView(APIView):
    def get(self, request):
        return Response(
            build_axis_recurrent_meta()
        )
