from django.db.models import Exists, OuterRef, Prefetch

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework.generics import ListAPIView
from rest_framework.exceptions import ValidationError

from database.models import FilterField, RNAInteraction, InteractionDatabase
from database.serializers.cerna_axis_serializers import FilterFieldSerializer, RNAInteractionSerializer
from database.pagination.standard_pagination import StandardPageNumberPagination


class FilterOptionsView(APIView):
    def get(self, request):
        table_name = "rna_interaction"

        queryset = (
            FilterField.objects
            .filter(
                table_name=table_name,
                is_active=True,
            )
            .prefetch_related("options")
            .order_by("sort_order", "field_name")
        )

        serializer = FilterFieldSerializer(queryset, many=True)

        return Response({
            "table_name": table_name,
            "fields": serializer.data,
        })


class RNAInteractionSearchView(ListAPIView):
    serializer_class = RNAInteractionSerializer
    pagination_class = StandardPageNumberPagination

    SEARCH_FIELD_MAP = {
        "miRNA": "source__name__icontains",
        "ceRNA": "target__name__icontains",
    }

    FILTER_FIELD_MAP = {
        "species": "species__in",
        "type": "interaction_type__in",
    }

    SORT_FIELD_MAP = {
        "miRNA": "source__name",
        "ceRNA": "target__name",
        "species": "species",
        "type": "interaction_type",
    }

    def post(self, request, *args, **kwargs):
        return self.list(request, *args, **kwargs)

    def get_queryset(self):
        data = self.request.data

        search_field = data.get("search_field")
        search_value = data.get("search_value", "").strip()

        filters = data.get("filters", {}) or {}

        sort_field = data.get("sort_field")
        sort_order = data.get("sort_order")

        queryset = (
            RNAInteraction.objects
            .select_related("source", "target")
            .prefetch_related(
                Prefetch(
                    "databases",
                    queryset=InteractionDatabase.objects.only("id", "name"),
                )
            )
        )

        # 搜索
        if search_field and search_value:
            lookup = self.SEARCH_FIELD_MAP.get(search_field)

            if lookup is None:
                raise ValidationError({
                    "search_field": f"Unsupported search field: {search_field}"
                })

            queryset = queryset.filter(**{
                lookup: search_value
            })

        # 普通字段筛选
        for field, values in filters.items():
            if not values:
                continue

            if field == "database":
                continue

            lookup = self.FILTER_FIELD_MAP.get(field)

            if lookup is None:
                raise ValidationError({
                    "filters": f"Unsupported filter field: {field}"
                })

            queryset = queryset.filter(**{
                lookup: values
            })

        # database ManyToMany 筛选：用 Exists，避免 JOIN 后 distinct
        database_values = filters.get("database")

        if database_values:
            through_model = RNAInteraction.databases.through

            queryset = queryset.filter(
                Exists(
                    through_model.objects.filter(
                        rnainteraction_id=OuterRef("pk"),
                        interactiondatabase__name__in=database_values,
                    )
                )
            )

        # 排序
        if sort_field:
            sort_lookup = self.SORT_FIELD_MAP.get(sort_field)

            if sort_lookup is None:
                raise ValidationError({
                    "sort_field": f"Unsupported sort field: {sort_field}"
                })

            if sort_order == "descend":
                sort_lookup = f"-{sort_lookup}"

            queryset = queryset.order_by(sort_lookup, "id")
        else:
            queryset = queryset.order_by("id")

        return queryset
