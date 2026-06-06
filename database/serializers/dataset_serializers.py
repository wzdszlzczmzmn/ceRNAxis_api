from rest_framework import serializers

from database.models import DatasetMetadata


class DatasetMetadataSerializer(serializers.ModelSerializer):
    class Meta:
        model = DatasetMetadata
        fields = [
            "id",
            "dataset",
            "programme",
            "obs_type",
            "reference",
            "cancer_type",
            "cancer_type_full_name",
            "gene_bio_type",
            "workflow",
            "sample_nums",
        ]
