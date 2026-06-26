import pandas as pd
import numpy as np
import pyarrow.parquet as pq
from django.http import FileResponse

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from database.models import DatasetMetadata
from database.serializers.dataset_serializers import DatasetMetadataSerializer
from database.utils.expression_file_utils import (
    MAX_SELECTED_GENES,
    DEFAULT_EXPRESSION_FILE_FORMAT,
    ExpressionPathError,
    get_available_expression_types, get_available_aliquot_expression_files, ALIQUOT_EXPRESSION_FILE_FORMAT,
    get_forced_isoform_dataset_from_mirna_dataset, get_aliquot_expression_file_path, get_aliquot_isoform_info_file_path,
    get_default_timedb_expression_file_format, get_available_timedb_expression_types, get_expression_mode_from_metadata,
    resolve_dataset_expression_file, read_expression_columns, read_expression_data,
)
from database.utils.meta_file_utils import get_dataset_meta_file, DatasetMetaPathError
from database.utils.viz_file_utils import get_available_deg_expression_types, DEGPathError, validate_deg_file, \
    get_available_timedb_deg_expression_types, validate_dataset_deg_file


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
            .order_by("id")
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
        expression_mode = "tcga" if metadata.programme == "TCGA" else "timedb"

        try:
            if expression_mode == "tcga":
                expression_file_format = DEFAULT_EXPRESSION_FILE_FORMAT

                available_expression_types = get_available_expression_types(
                    dataset=metadata.dataset,
                    rna_type=rna_type,
                    file_format=expression_file_format,
                )

                available_deg_expression_types = (
                    get_available_deg_expression_types(
                        dataset=metadata.dataset,
                        rna_type=rna_type,
                    )
                )

            else:
                expression_file_format = (
                    get_default_timedb_expression_file_format(
                        dataset=metadata.dataset,
                        rna_type=rna_type,
                    )
                )

                available_expression_types = (
                    get_available_timedb_expression_types(
                        dataset=metadata.dataset,
                        rna_type=rna_type,
                        file_format=expression_file_format,
                    )
                )

                available_deg_expression_types = (
                    get_available_timedb_deg_expression_types(
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
                "expression_mode": expression_mode,
                "expression_file_format": expression_file_format,
                "available_expression_types": available_expression_types,
                "available_deg_expression_types": (
                    available_deg_expression_types
                ),
            },
            status=status.HTTP_200_OK,
        )


class DatasetSampleMetaView(APIView):
    def get(self, request, dataset):
        try:
            metadata = DatasetMetadata.objects.get(dataset=dataset)
        except DatasetMetadata.DoesNotExist:
            return Response(
                {"detail": "Dataset metadata not found."},
                status=status.HTTP_404_NOT_FOUND,
            )

        expression_mode = "tcga" if metadata.programme == "TCGA" else "timedb"

        try:
            file_path = get_dataset_meta_file(
                dataset=metadata.dataset,
                expression_mode=expression_mode,
            )
        except DatasetMetaPathError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )

        if not file_path.exists() or not file_path.is_file():
            return Response(
                {
                    "detail": (
                        f"Sample metadata file not found for dataset "
                        f"'{metadata.dataset}'."
                    )
                },
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
        expression_mode = get_expression_mode_from_metadata(metadata)

        try:
            file_path, file_format = resolve_dataset_expression_file(
                dataset=metadata.dataset,
                rna_type=rna_type,
                expression_type=expression_type,
                expression_mode=expression_mode,
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
            columns = read_expression_columns(
                file_path=file_path,
                file_format=file_format,
            )
        except ExpressionPathError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )
        except Exception as e:
            return Response(
                {
                    "detail": (
                        f"Failed to read expression {file_format} schema: {e}"
                    )
                },
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
                "dataset": metadata.dataset,
                "rna_type": rna_type,
                "expression_mode": expression_mode,
                "expression_type": expression_type,
                "file_format": file_format,
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

        genes = [
            str(gene).strip()
            for gene in genes
            if str(gene).strip()
        ]

        if not genes:
            return Response(
                {"detail": "At least one gene must be selected."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        if len(genes) > MAX_SELECTED_GENES:
            return Response(
                {
                    "detail": (
                        f"At most {MAX_SELECTED_GENES} genes can be selected."
                    )
                },
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
        expression_mode = get_expression_mode_from_metadata(metadata)

        try:
            file_path, file_format = resolve_dataset_expression_file(
                dataset=metadata.dataset,
                rna_type=rna_type,
                expression_type=expression_type,
                expression_mode=expression_mode,
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
            all_columns = read_expression_columns(
                file_path=file_path,
                file_format=file_format,
            )
        except ExpressionPathError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )
        except Exception as e:
            return Response(
                {
                    "detail": (
                        f"Failed to read expression {file_format} schema: {e}"
                    )
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        if not all_columns:
            return Response(
                {"detail": "Expression file has no columns."},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        sample_column = all_columns[0]
        available_genes = set(all_columns[1:])

        missing_genes = [
            gene for gene in genes
            if gene not in available_genes
        ]

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
            df = read_expression_data(
                file_path=file_path,
                file_format=file_format,
                usecols=usecols,
            )
        except ExpressionPathError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )
        except Exception as e:
            return Response(
                {
                    "detail": (
                        f"Failed to read expression {file_format} data: {e}"
                    )
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        df = df.fillna("")

        return Response(
            {
                "dataset": metadata.dataset,
                "rna_type": rna_type,
                "expression_mode": expression_mode,
                "expression_type": expression_type,
                "file_format": file_format,
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
        "pvalue",
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
            deg_file = validate_dataset_deg_file(
                metadata=metadata,
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
                    "detail": (
                        "Missing required columns: "
                        f"{sorted(missing_columns)}"
                    )
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        df = df[
            [
                "gene_name",
                "log2FC",
                "pvalue",
                "regulation",
            ]
        ].copy()

        df = df.replace([np.inf, -np.inf], np.nan)

        raw_count = int(df.shape[0])

        df = df.dropna(
            subset=[
                "gene_name",
                "log2FC",
                "pvalue",
                "regulation",
            ]
        )

        df["gene_name"] = df["gene_name"].astype(str).str.strip()
        df["regulation"] = df["regulation"].astype(str).str.strip()

        df["log2FC"] = pd.to_numeric(
            df["log2FC"],
            errors="coerce",
        )

        df["pvalue"] = pd.to_numeric(
            df["pvalue"],
            errors="coerce",
        )

        df = df.dropna(
            subset=[
                "gene_name",
                "log2FC",
                "pvalue",
                "regulation",
            ]
        )

        df = df[df["pvalue"] > 0]

        df = df[
            df["regulation"].isin(self.VALID_REGULATION_GROUPS)
        ]

        cleaned_count = int(df.shape[0])
        dropped_count = raw_count - cleaned_count

        df["neg_log10_pvalue"] = -np.log10(df["pvalue"])

        groups = {}

        for group in self.VALID_REGULATION_GROUPS:
            sub_df = df[df["regulation"] == group]

            groups[group] = [
                {
                    "gene_name": row["gene_name"],
                    "log2FC": float(row["log2FC"]),
                    "pvalue": float(row["pvalue"]),
                    "neg_log10_pvalue": float(row["neg_log10_pvalue"]),
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


class DatasetAliquotExpressionFileListView(APIView):
    def get(self, request):
        dataset = request.query_params.get("dataset")

        if not dataset:
            return Response(
                {"detail": "Missing required query parameter: dataset."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        try:
            metadata = DatasetMetadata.objects.get(dataset=dataset)
        except DatasetMetadata.DoesNotExist:
            return Response(
                {"detail": "Dataset metadata not found."},
                status=status.HTTP_404_NOT_FOUND,
            )

        rna_type = metadata.gene_bio_type

        try:
            available_aliquot_expression_files = (
                get_available_aliquot_expression_files(
                    dataset=metadata.dataset,
                    rna_type=rna_type,
                )
            )

            available_isoform_files = []

            if rna_type == "miRNA":
                isoform_dataset = get_forced_isoform_dataset_from_mirna_dataset(
                    dataset=metadata.dataset,
                )

                available_isoform_files = get_available_aliquot_expression_files(
                    dataset=isoform_dataset,
                    rna_type="isoform",
                )

        except ExpressionPathError as exc:
            return Response(
                {"detail": str(exc)},
                status=status.HTTP_400_BAD_REQUEST,
            )

        return Response(
            {
                "aliquot_expression_file_format": ALIQUOT_EXPRESSION_FILE_FORMAT,
                "available_aliquot_expression_files": (
                    available_aliquot_expression_files
                ),
                "available_isoform_files": available_isoform_files,
            },
            status=status.HTTP_200_OK,
        )


class DatasetAliquotExpressionFileDownloadView(APIView):
    def get(self, request):
        dataset = request.query_params.get("dataset")
        value_type = request.query_params.get("value_type")
        file_type = request.query_params.get("file_type", "expression")

        if not dataset:
            return Response(
                {"detail": "Missing required query parameter: dataset."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        if file_type not in {"expression", "annotation"}:
            return Response(
                {
                    "detail": (
                        "Invalid file_type. Allowed values: "
                        "['expression', 'annotation']."
                    )
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        try:
            if dataset.endswith("_isoform"):
                rna_type = "isoform"
                resolved_dataset = dataset
            else:
                try:
                    metadata = DatasetMetadata.objects.get(dataset=dataset)
                except DatasetMetadata.DoesNotExist:
                    return Response(
                        {"detail": "Dataset metadata not found."},
                        status=status.HTTP_404_NOT_FOUND,
                    )

                rna_type = metadata.gene_bio_type
                resolved_dataset = metadata.dataset

            if file_type == "expression":
                if not value_type:
                    return Response(
                        {
                            "detail": (
                                "Missing required query parameter: value_type "
                                "for expression file download."
                            )
                        },
                        status=status.HTTP_400_BAD_REQUEST,
                    )

                file_path = get_aliquot_expression_file_path(
                    dataset=resolved_dataset,
                    rna_type=rna_type,
                    value_type=value_type,
                )

            else:
                if rna_type != "isoform":
                    return Response(
                        {
                            "detail": (
                                "Annotation file is only available for "
                                "isoform files."
                            )
                        },
                        status=status.HTTP_400_BAD_REQUEST,
                    )

                file_path = get_aliquot_isoform_info_file_path(
                    dataset=resolved_dataset,
                )

            if not file_path.exists() or not file_path.is_file():
                return Response(
                    {"detail": f"Aliquot file not found: {file_path.name}"},
                    status=status.HTTP_404_NOT_FOUND,
                )

        except ExpressionPathError as exc:
            return Response(
                {"detail": str(exc)},
                status=status.HTTP_400_BAD_REQUEST,
            )

        return FileResponse(
            open(file_path, "rb"),
            as_attachment=True,
            filename=file_path.name,
            content_type="text/csv",
        )
