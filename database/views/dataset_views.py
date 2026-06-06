import pandas as pd

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from database.models import DatasetMetadata
from database.serializers.dataset_serializers import DatasetMetadataSerializer
from database.utils.expression_file_utils import get_available_expression_types, EXPRESSION_FILE_TYPES, \
    validate_expression_file, MAX_SELECTED_GENES
from database.utils.meta_file_utils import get_dataset_meta_file


class DatasetMetadataListView(APIView):
    VALID_GENE_BIO_TYPES = {"miRNA", "mRNA", "lncRNA", "circRNA"}

    def get(self, request):
        gene_bio_type = request.query_params.get("gene_bio_type")

        if not gene_bio_type:
            return Response(
                {
                    "detail": "Missing query parameter: gene_bio_type."
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        if gene_bio_type not in self.VALID_GENE_BIO_TYPES:
            return Response(
                {
                    "detail": (
                        "Invalid gene_bio_type. "
                        "Allowed values are: miRNA, mRNA, lncRNA, circRNA."
                    )
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        queryset = (
            DatasetMetadata.objects
            .filter(gene_bio_type=gene_bio_type)
            .order_by("programme", "cancer_type")
        )

        serializer = DatasetMetadataSerializer(queryset, many=True)

        return Response(
            {
                "gene_bio_type": gene_bio_type,
                "count": queryset.count(),
                "results": serializer.data,
            }
        )


class DatasetMetadataDetailView(APIView):
    def get(self, request, dataset):
        try:
            metadata = DatasetMetadata.objects.get(dataset=dataset)
        except DatasetMetadata.DoesNotExist:
            return Response(
                {"detail": "Dataset metadata not found."},
                status=status.HTTP_404_NOT_FOUND,
            )

        serializer = DatasetMetadataSerializer(metadata)

        return Response(
            {
                "metadata": serializer.data,
                "available_expression_types": get_available_expression_types(
                    dataset=metadata.dataset
                ),
            },
            status=status.HTTP_200_OK,
        )


class DatasetSampleMetaView(APIView):
    def get(self, request, dataset):
        file_path = get_dataset_meta_file(dataset)

        if not file_path.exists() or not file_path.is_file():
            return Response(
                {"detail": f"Sample metadata file not found for dataset '{dataset}'."},
                status=status.HTTP_404_NOT_FOUND,
            )

        try:
            df = pd.read_csv(file_path, low_memory=False)
        except Exception as e:
            return Response(
                {"detail": f"Failed to read sample metadata file: {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        df = df.fillna("")

        return Response(
            {
                "count": len(df),
                "columns": df.columns.tolist(),
                "results": df.to_dict("records"),
            },
            status=status.HTTP_200_OK,
        )


class DatasetExpressionGeneListView(APIView):
    def get(self, request):
        dataset = request.query_params.get("dataset")
        expression_type = request.query_params.get("expression_type")

        if not dataset:
            return Response(
                {"detail": "Missing query parameter: dataset."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        if not expression_type:
            return Response(
                {"detail": "Missing query parameter: expression_type."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        if expression_type not in EXPRESSION_FILE_TYPES:
            return Response(
                {
                    "detail": (
                        "Invalid expression_type. "
                        "Allowed values are: log2count, log2fpkm, log2fpkmuq, log2tpm."
                    )
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        if not DatasetMetadata.objects.filter(dataset=dataset).exists():
            return Response(
                {"detail": f"Dataset metadata not found: {dataset}."},
                status=status.HTTP_404_NOT_FOUND,
            )

        try:
            file_path = validate_expression_file(dataset, expression_type)
        except FileNotFoundError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_404_NOT_FOUND,
            )

        try:
            df_header = pd.read_csv(file_path, nrows=0)
        except Exception as e:
            return Response(
                {"detail": f"Failed to read expression file header: {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        columns = df_header.columns.tolist()

        if not columns:
            return Response(
                {"detail": "Expression file has no columns."},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        # 默认第一列是 sample_id，其余列视为 gene
        genes = columns[1:]

        return Response(
            {
                "dataset": dataset,
                "expression_type": expression_type,
                "count": len(genes),
                "genes": genes,
            },
            status=status.HTTP_200_OK,
        )


class DatasetExpressionDataView(APIView):
    def post(self, request):
        dataset = request.data.get("dataset")
        expression_type = request.data.get("expression_type")
        genes = request.data.get("genes", [])

        if not dataset:
            return Response(
                {"detail": "Missing field: dataset."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        if not expression_type:
            return Response(
                {"detail": "Missing field: expression_type."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        if expression_type not in EXPRESSION_FILE_TYPES:
            return Response(
                {
                    "detail": (
                        "Invalid expression_type. "
                        "Allowed values are: log2count, log2fpkm, log2fpkmuq, log2tpm."
                    )
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        if not isinstance(genes, list):
            return Response(
                {"detail": "Field 'genes' must be a list."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        genes = [str(gene).strip() for gene in genes if str(gene).strip()]

        if not genes:
            return Response(
                {"detail": "At least one gene must be selected."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        if len(genes) > MAX_SELECTED_GENES:
            return Response(
                {"detail": f"At most {MAX_SELECTED_GENES} genes can be selected."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        if not DatasetMetadata.objects.filter(dataset=dataset).exists():
            return Response(
                {"detail": f"Dataset metadata not found: {dataset}."},
                status=status.HTTP_404_NOT_FOUND,
            )

        try:
            file_path = validate_expression_file(dataset, expression_type)
        except FileNotFoundError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_404_NOT_FOUND,
            )

        try:
            header_df = pd.read_csv(file_path, nrows=0)
            all_columns = header_df.columns.tolist()
        except Exception as e:
            return Response(
                {"detail": f"Failed to read expression file header: {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        if not all_columns:
            return Response(
                {"detail": "Expression file has no columns."},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        sample_column = all_columns[0]
        available_genes = set(all_columns[1:])

        missing_genes = [gene for gene in genes if gene not in available_genes]

        if missing_genes:
            return Response(
                {
                    "detail": "Some genes are not found in the expression file.",
                    "missing_genes": missing_genes,
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        usecols = [sample_column] + genes

        try:
            df = pd.read_csv(
                file_path,
                usecols=usecols,
                low_memory=False,
            )
        except Exception as e:
            return Response(
                {"detail": f"Failed to read expression data: {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        df = df.fillna("")

        return Response(
            {
                "dataset": dataset,
                "expression_type": expression_type,
                "count": len(df),
                "columns": df.columns.tolist(),
                "results": df.to_dict("records"),
            },
            status=status.HTTP_200_OK,
        )
