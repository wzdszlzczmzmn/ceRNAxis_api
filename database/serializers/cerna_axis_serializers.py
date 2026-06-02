from rest_framework import serializers
from database.models import FilterOption, FilterField, RNANode, RNAInteraction


class FilterOptionSerializer(serializers.ModelSerializer):
    class Meta:
        model = FilterOption
        fields = ["value", "sort_order"]


class FilterFieldSerializer(serializers.ModelSerializer):
    options = serializers.SerializerMethodField()

    class Meta:
        model = FilterField
        fields = [
            "id",
            "table_name",
            "field_name",
            "field_label",
            "sort_order",
            "options",
        ]

    def get_options(self, obj):
        options = obj.options.filter(is_active=True).order_by("sort_order", "value")
        return [option.value for option in options]


class RNANodeSerializer(serializers.ModelSerializer):
    class Meta:
        model = RNANode
        fields = [
            "id",
            "name",
            "rna_type",
            "species",
        ]


class RNAInteractionSerializer(serializers.ModelSerializer):
    miRNA = serializers.CharField(source="source.name", read_only=True)
    ceRNA = serializers.CharField(source="target.name", read_only=True)
    ceRNA_type = serializers.CharField(source="target.rna_type", read_only=True)
    type = serializers.CharField(source="interaction_type", read_only=True)

    databases = serializers.SlugRelatedField(
        many=True,
        read_only=True,
        slug_field="name",
    )

    class Meta:
        model = RNAInteraction
        fields = [
            "id",
            "miRNA",
            "ceRNA",
            "ceRNA_type",
            "species",
            "type",
            "databases",
        ]
