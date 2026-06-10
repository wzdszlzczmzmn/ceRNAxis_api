import pandas as pd
import numpy as np
import pyarrow.parquet as pq

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from database.models import DatasetMetadata
from database.serializers.dataset_serializers import DatasetMetadataSerializer
from database.utils.expression_file_utils import (
    MAX_SELECTED_GENES,
    DEFAULT_EXPRESSION_FILE_FORMAT,
    ExpressionPathError,
    validate_expression_type,
    validate_expression_file,
    get_available_expression_types,
)
from database.utils.meta_file_utils import get_dataset_meta_file
from database.utils.viz_file_utils import get_available_deg_expression_types, DEGPathError, validate_deg_file


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

        rna_type = metadata.gene_bio_type

        try:
            available_expression_types = get_available_expression_types(
                dataset=metadata.dataset,
                rna_type=rna_type,
                file_format=DEFAULT_EXPRESSION_FILE_FORMAT,
            )

            available_deg_expression_types = (
                get_available_deg_expression_types(
                    dataset=metadata.dataset,
                    rna_type=rna_type,
                )
            )

        except (ExpressionPathError, DEGPathError) as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )

        return Response(
            {
                "metadata": serializer.data,
                "expression_file_format": DEFAULT_EXPRESSION_FILE_FORMAT,
                "available_expression_types": available_expression_types,
                "available_deg_expression_types": (
                    available_deg_expression_types
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

        try:
            metadata = DatasetMetadata.objects.get(dataset=dataset)
        except DatasetMetadata.DoesNotExist:
            return Response(
                {"detail": f"Dataset metadata not found: {dataset}."},
                status=status.HTTP_404_NOT_FOUND,
            )

        rna_type = metadata.gene_bio_type

        try:
            validate_expression_type(
                rna_type=rna_type,
                expression_type=expression_type,
            )

            file_path = validate_expression_file(
                dataset=dataset,
                rna_type=rna_type,
                expression_type=expression_type,
            )
        except ExpressionPathError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )
        except FileNotFoundError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_404_NOT_FOUND,
            )

        try:
            schema = pq.read_schema(file_path)
            columns = schema.names
        except Exception as e:
            return Response(
                {"detail": f"Failed to read expression parquet schema: {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        if not columns:
            return Response(
                {"detail": "Expression file has no columns."},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        sample_column = columns[0]
        genes = columns[1:]

        return Response(
            {
                "dataset": dataset,
                "rna_type": rna_type,
                "expression_type": expression_type,
                "file_format": DEFAULT_EXPRESSION_FILE_FORMAT,
                "sample_column": sample_column,
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

        try:
            metadata = DatasetMetadata.objects.get(dataset=dataset)
        except DatasetMetadata.DoesNotExist:
            return Response(
                {"detail": f"Dataset metadata not found: {dataset}."},
                status=status.HTTP_404_NOT_FOUND,
            )

        rna_type = metadata.gene_bio_type

        try:
            validate_expression_type(
                rna_type=rna_type,
                expression_type=expression_type,
            )

            file_path = validate_expression_file(
                dataset=dataset,
                rna_type=rna_type,
                expression_type=expression_type,
            )
        except ExpressionPathError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )
        except FileNotFoundError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_404_NOT_FOUND,
            )

        try:
            schema = pq.read_schema(file_path)
            all_columns = schema.names
        except Exception as e:
            return Response(
                {"detail": f"Failed to read expression parquet schema: {e}"},
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
            df = pd.read_parquet(
                file_path,
                columns=usecols,
                engine="pyarrow",
            )
        except Exception as e:
            return Response(
                {"detail": f"Failed to read expression parquet data: {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        df = df.fillna("")

        return Response(
            {
                "dataset": dataset,
                "rna_type": rna_type,
                "expression_type": expression_type,
                "file_format": DEFAULT_EXPRESSION_FILE_FORMAT,
                "count": len(df),
                "columns": df.columns.tolist(),
                "results": df.to_dict("records"),
            },
            status=status.HTTP_200_OK,
        )


class DatasetDEGVolcanoView(APIView):
    REQUIRED_COLUMNS = {
        "gene_name",
        "log2FC",
        "padj",
        "regulation",
    }

    VALID_REGULATION_GROUPS = ["NotSig", "Down", "Up"]

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

        try:
            metadata = DatasetMetadata.objects.get(dataset=dataset)
        except DatasetMetadata.DoesNotExist:
            return Response(
                {"detail": f"Dataset metadata not found: {dataset}."},
                status=status.HTTP_404_NOT_FOUND,
            )

        rna_type = metadata.gene_bio_type

        try:
            deg_file = validate_deg_file(
                dataset=dataset,
                rna_type=rna_type,
                expression_type=expression_type,
            )
        except DEGPathError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )
        except FileNotFoundError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_404_NOT_FOUND,
            )

        try:
            df = pd.read_csv(deg_file)
        except Exception as e:
            return Response(
                {"detail": f"Failed to read DEG file: {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        missing_columns = self.REQUIRED_COLUMNS - set(df.columns)

        if missing_columns:
            return Response(
                {
                    "detail": f"Missing required columns: {sorted(missing_columns)}"
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        df = df[
            [
                "gene_name",
                "log2FC",
                "padj",
                "regulation",
            ]
        ].copy()

        df = df.replace([np.inf, -np.inf], np.nan)

        raw_count = int(df.shape[0])

        df = df.dropna(
            subset=[
                "gene_name",
                "log2FC",
                "padj",
                "regulation",
            ]
        )

        df = df[df["padj"] > 0]

        df = df[
            df["regulation"].isin(self.VALID_REGULATION_GROUPS)
        ]

        cleaned_count = int(df.shape[0])
        dropped_count = raw_count - cleaned_count

        df["neg_log10_padj"] = -np.log10(df["padj"])

        groups = {}

        for group in self.VALID_REGULATION_GROUPS:
            sub_df = df[df["regulation"] == group]

            groups[group] = [
                {
                    "gene_name": row["gene_name"],
                    "log2FC": float(row["log2FC"]),
                    "padj": float(row["padj"]),
                    "neg_log10_padj": float(row["neg_log10_padj"]),
                }
                for _, row in sub_df.iterrows()
            ]

        return Response(
            {
                "dataset": dataset,
                "rna_type": rna_type,
                "expression_type": expression_type,
                "summary": {
                    "raw_count": raw_count,
                    "cleaned_count": cleaned_count,
                    "dropped_count": dropped_count,
                    "not_sig": len(groups["NotSig"]),
                    "down": len(groups["Down"]),
                    "up": len(groups["Up"]),
                },
                "groups": groups,
            },
            status=status.HTTP_200_OK,
        )
