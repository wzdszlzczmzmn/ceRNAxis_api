from __future__ import annotations

from django.db.models import Q

from rest_framework.generics import ListAPIView
from rest_framework.response import Response
from rest_framework.views import APIView

from database.models import AxisStructureRecurrentSummary
from database.pagination.standard_pagination import (
    StandardPageNumberPagination,
)
from database.serializers.axis_recurrent_serializers import (
    AxisRecurrentSearchRequestSerializer,
    AxisStructureRecurrentSummarySerializer,
)
from database.services.axis_recurrent import (
    apply_axis_recurrent_pattern,
)
from database.utils.axis_recurrent_meta_utils import (
    build_axis_recurrent_meta,
)


class AxisRecurrentSummarySearchView(ListAPIView):
    """
    Search the rebuilt recurrent Axis summaries.

    Return records already materialized in the recurrent summary tables.

    Axis Final summary data is included by the response serializer whenever
    the corresponding optional OneToOne summary exists.

    Request example:
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
                "min_dataset_count": 2,
                "has_axis_final": true,
                "regulation_consistent": [
                    true
                ]
            },
            "sort_field": "dataset_count",
            "sort_order": "descend"
        }
    """

    serializer_class = (
        AxisStructureRecurrentSummarySerializer
    )
    pagination_class = StandardPageNumberPagination

    SORT_FIELD_MAP = {
        "axis_signature":
            "axis__axis_signature",
        "axis_type":
            "axis__axis_type",

        "miRNA":
            "axis__miRNA",
        "mRNA":
            "axis__mRNA",
        "lncRNA":
            "axis__lncRNA",
        "circRNA":
            "axis__circRNA",

        "dataset_count":
            "dataset_count",
        "context_count":
            "context_count",

        "tcga_dataset_count":
            "tcga_dataset_count",
        "timedb_dataset_count":
            "timedb_dataset_count",
        "sc_dataset_count":
            "sc_dataset_count",
        "st_dataset_count":
            "st_dataset_count",

        "axis_final_context_count":
            "axis_final_context_count",
        "sponge_context_count":
            "sponge_context_count",
        "both_result_context_count":
            "both_result_context_count",

        "regulation_pattern_count": (
            "axis__axis_final_recurrent_summary__"
            "regulation_pattern_count"
        ),
        "dominant_regulation_count": (
            "axis__axis_final_recurrent_summary__"
            "dominant_regulation_count"
        ),

        "updated_at":
            "updated_at",
    }

    def post(self, request, *args, **kwargs):
        return self.list(
            request,
            *args,
            **kwargs,
        )

    def get_search_data(self) -> dict:
        cached = getattr(
            self,
            "_validated_search_data",
            None,
        )

        if cached is not None:
            return cached

        serializer = AxisRecurrentSearchRequestSerializer(
            data=self.request.data,
        )
        serializer.is_valid(
            raise_exception=True,
        )

        self._validated_search_data = (
            serializer.validated_data
        )

        return self._validated_search_data

    def apply_source_filter(
            self,
            queryset,
            sources,
    ):
        if not sources:
            return queryset

        source_query = Q()

        if "TCGA" in sources:
            source_query |= Q(
                tcga_dataset_count__gt=0,
            )

        if "TIMEDB" in sources:
            source_query |= Q(
                timedb_dataset_count__gt=0,
            )

        if "SC" in sources:
            source_query |= Q(
                sc_dataset_count__gt=0,
            )

        if "ST" in sources:
            source_query |= Q(
                st_dataset_count__gt=0,
            )

        return queryset.filter(
            source_query
        )

    def apply_filters(
            self,
            queryset,
            filters: dict,
    ):
        axis_types = filters.get(
            "axis_type",
            [],
        )

        if axis_types:
            queryset = queryset.filter(
                axis__axis_type__in=axis_types,
            )

        queryset = self.apply_source_filter(
            queryset,
            filters.get("source", []),
        )

        min_dataset_count = filters.get(
            "min_dataset_count"
        )

        if min_dataset_count is not None:
            queryset = queryset.filter(
                dataset_count__gte=min_dataset_count,
            )

        min_context_count = filters.get(
            "min_context_count"
        )

        if min_context_count is not None:
            queryset = queryset.filter(
                context_count__gte=min_context_count,
            )

        has_axis_final = filters.get(
            "has_axis_final"
        )

        if has_axis_final is True:
            queryset = queryset.filter(
                axis_final_context_count__gt=0,
            )
        elif has_axis_final is False:
            queryset = queryset.filter(
                axis_final_context_count=0,
            )

        has_sponge = filters.get(
            "has_sponge"
        )

        if has_sponge is True:
            queryset = queryset.filter(
                sponge_context_count__gt=0,
            )
        elif has_sponge is False:
            queryset = queryset.filter(
                sponge_context_count=0,
            )

        has_both_result_context = filters.get(
            "has_both_result_context"
        )

        if has_both_result_context is True:
            queryset = queryset.filter(
                both_result_context_count__gt=0,
            )
        elif has_both_result_context is False:
            queryset = queryset.filter(
                both_result_context_count=0,
            )

        dominant_regulations = filters.get(
            "dominant_axis_regulation",
            [],
        )

        if dominant_regulations:
            queryset = queryset.filter(**{
                (
                    "axis__axis_final_recurrent_summary__"
                    "dominant_axis_regulation__in"
                ): dominant_regulations,
            })

        consistency_values = filters.get(
            "regulation_consistent",
            [],
        )

        if consistency_values:
            queryset = queryset.filter(**{
                (
                    "axis__axis_final_recurrent_summary__"
                    "regulation_consistent__in"
                ): consistency_values,
            })

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
                "-dataset_count",
                "-context_count",
                "axis__axis_signature",
            )

        sort_lookup = self.SORT_FIELD_MAP[
            sort_field
        ]

        if sort_order == "descend":
            sort_lookup = f"-{sort_lookup}"

        return queryset.order_by(
            sort_lookup,
            "axis__axis_signature",
        )

    def get_queryset(self):
        data = self.get_search_data()

        pattern = data.get(
            "pattern",
            "",
        )

        filters = data.get(
            "filters",
            {},
        )

        queryset = (
            AxisStructureRecurrentSummary.objects
            .select_related(
                "axis",
                "axis__axis_final_recurrent_summary",
            )
            .all()
        )

        queryset = apply_axis_recurrent_pattern(
            queryset,
            pattern,
            field_prefix="axis__",
        )

        queryset = self.apply_filters(
            queryset,
            filters,
        )

        queryset = self.apply_sorting(
            queryset,
            sort_field=data.get("sort_field"),
            sort_order=data.get("sort_order"),
        )

        return queryset


class AxisRecurrentMetaView(APIView):
    def get(self, request):
        return Response(
            build_axis_recurrent_meta()
        )
