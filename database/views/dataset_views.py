import pandas as pd
import numpy as np
import pyarrow.parquet as pq
import pyarrow as pa
from django.http import FileResponse

from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status

from database.models import DatasetMetadata
from database.serializers.dataset_serializers import DatasetMetadataSerializer
from database.services.dataset_download.service import prepare_dataset_download, prepare_tcga_annotation_download, \
    prepare_timedb_annotation_download
from database.utils.expression_file_utils import (
    MAX_SELECTED_GENES,
    DEFAULT_EXPRESSION_FILE_FORMAT,
    ExpressionPathError,
    get_available_expression_types, get_available_aliquot_expression_files, ALIQUOT_EXPRESSION_FILE_FORMAT,
    get_forced_isoform_dataset_from_mirna_dataset, get_aliquot_expression_file_path, get_aliquot_isoform_info_file_path,
    get_default_timedb_expression_file_format, get_available_timedb_expression_types, get_expression_mode_from_metadata,
    resolve_dataset_expression_file, read_expression_columns, read_expression_data,
    get_available_tisch2_expression_types, TISCH2_EXPRESSION_FILE_FORMAT, get_available_sctml_expression_types,
    SCTML_EXPRESSION_FILE_FORMAT,
)
from database.utils.meta_file_utils import get_dataset_meta_file, DatasetMetaPathError, get_large_meta_file
from database.utils.viz_file_utils import get_available_deg_expression_types, DEGPathError, \
    get_available_timedb_deg_expression_types, validate_dataset_deg_file, get_available_tisch2_deg_expression_types, \
    validate_tisch2_deg_file


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
        expression_mode = get_expression_mode_from_metadata(metadata)

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

            elif expression_mode == "timedb":
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

            elif expression_mode == "scTML":
                expression_file_format = SCTML_EXPRESSION_FILE_FORMAT

                available_expression_types = (
                    get_available_sctml_expression_types(
                        dataset=metadata.dataset,
                        rna_type=rna_type,
                    )
                )

                available_deg_expression_types = []

            elif expression_mode == "tisch2":
                expression_file_format = TISCH2_EXPRESSION_FILE_FORMAT

                available_expression_types = (
                    get_available_tisch2_expression_types(
                        dataset=metadata.dataset,
                        rna_type=rna_type,
                    )
                )

                available_deg_expression_types = (
                    get_available_tisch2_deg_expression_types(
                        dataset=metadata.dataset,
                        rna_type=rna_type,
                    )
                )

            else:
                return Response(
                    {
                        "detail": (
                            f"Unsupported expression_mode: {expression_mode}"
                        )
                    },
                    status=status.HTTP_400_BAD_REQUEST,
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

        expression_mode = get_expression_mode_from_metadata(metadata)

        if expression_mode not in {"tcga", "timedb"}:
            return Response(
                {
                    "detail": (
                        f"Sample metadata table is not available for "
                        f"expression_mode '{expression_mode}'."
                    )
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

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


class DatasetLargeMetaView(APIView):
    DEFAULT_PAGE = 1
    DEFAULT_PAGE_SIZE = 50
    MAX_PAGE_SIZE = 200

    ALLOWED_EXPRESSION_MODES = {"tisch2", "scTML"}

    def parse_positive_int(
        self,
        value,
        default: int,
        field_name: str,
        max_value: int | None = None,
    ) -> int:
        if value in [None, ""]:
            return default

        try:
            parsed = int(value)
        except (TypeError, ValueError):
            raise ValueError(f"Invalid {field_name}. It must be an integer.")

        if parsed < 1:
            raise ValueError(f"Invalid {field_name}. It must be >= 1.")

        if max_value is not None and parsed > max_value:
            raise ValueError(
                f"Invalid {field_name}. It must be <= {max_value}."
            )

        return parsed

    def infer_delimiter(self, file_path):
        with open(file_path, "r", encoding="utf-8-sig") as f:
            first_line = f.readline()

        if "\t" in first_line:
            return "\t"

        return ","

    def count_csv_rows(self, file_path) -> int:
        with open(file_path, "r", encoding="utf-8-sig") as f:
            return max(sum(1 for _ in f) - 1, 0)

    def get(self, request, dataset):
        if not dataset:
            return Response(
                {"detail": "Missing query parameter: dataset."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        try:
            page = self.parse_positive_int(
                request.query_params.get("page"),
                default=self.DEFAULT_PAGE,
                field_name="page",
            )
            page_size = self.parse_positive_int(
                request.query_params.get("page_size"),
                default=self.DEFAULT_PAGE_SIZE,
                field_name="page_size",
                max_value=self.MAX_PAGE_SIZE,
            )
        except ValueError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )

        try:
            metadata = DatasetMetadata.objects.get(dataset=dataset)
        except DatasetMetadata.DoesNotExist:
            return Response(
                {"detail": "Dataset metadata not found."},
                status=status.HTTP_404_NOT_FOUND,
            )

        expression_mode = get_expression_mode_from_metadata(metadata)

        if expression_mode not in self.ALLOWED_EXPRESSION_MODES:
            return Response(
                {
                    "detail": (
                        "Large metadata table is only available for "
                        "TISCH2 and scTML datasets."
                    )
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        try:
            file_path = get_large_meta_file(
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
                        f"Large metadata file not found for dataset "
                        f"'{metadata.dataset}'."
                    )
                },
                status=status.HTTP_404_NOT_FOUND,
            )

        try:
            delimiter = self.infer_delimiter(file_path)

            header_df = pd.read_csv(
                file_path,
                sep=delimiter,
                nrows=0,
                encoding="utf-8-sig",
            )
            columns = header_df.columns.tolist()

            row_count = self.count_csv_rows(file_path)
            row_start = (page - 1) * page_size

            if row_start >= row_count and row_count > 0:
                return Response(
                    {
                        "detail": "Page is out of range.",
                        "count": row_count,
                    },
                    status=status.HTTP_400_BAD_REQUEST,
                )

            df = pd.read_csv(
                file_path,
                sep=delimiter,
                skiprows=range(1, row_start + 1),
                nrows=page_size,
                encoding="utf-8-sig",
                low_memory=False,
            )

            df = df.fillna("")

        except Exception as e:
            return Response(
                {"detail": f"Failed to read large metadata file: {e}"},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        return Response(
            {
                "dataset": metadata.dataset,
                "expression_mode": expression_mode,
                "file_format": "csv",
                "count": row_count,
                "page": page,
                "page_size": page_size,
                "columns": columns,
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

        all_columns = columns

        if expression_mode in {"tisch2", "scTML"}:
            sample_column = None
            genes = all_columns
        else:
            sample_column = all_columns[0]
            genes = all_columns[1:]

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
        raw_genes = request.data.get("genes", [])

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

        try:
            genes = normalize_gene_list(raw_genes)
        except ValueError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )

        gene_error = validate_selected_genes(genes)
        if gene_error is not None:
            return gene_error

        try:
            metadata = DatasetMetadata.objects.get(dataset=dataset)
        except DatasetMetadata.DoesNotExist:
            return Response(
                {"detail": f"Dataset metadata not found: {dataset}."},
                status=status.HTTP_404_NOT_FOUND,
            )

        rna_type = metadata.gene_bio_type
        expression_mode = get_expression_mode_from_metadata(metadata)

        if expression_mode in LARGE_EXPRESSION_MODES:
            return Response(
                {
                    "detail": (
                        "DatasetExpressionDataView does not support "
                        f"expression_mode '{expression_mode}' because the "
                        "expression matrix is too large. Use the paginated "
                        "large expression data endpoint instead."
                    )
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

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

        gene_exist_error = validate_genes_exist(
            genes=genes,
            available_genes=available_genes,
        )
        if gene_exist_error is not None:
            return gene_exist_error

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


LARGE_EXPRESSION_MODES = {"tisch2", "sctml"}

DEFAULT_LARGE_EXPRESSION_PAGE = 1
DEFAULT_LARGE_EXPRESSION_PAGE_SIZE = 50
MAX_LARGE_EXPRESSION_PAGE_SIZE = 500


def parse_positive_int(
    value,
    default: int,
    field_name: str,
    max_value: int | None = None,
) -> int:
    if value in [None, ""]:
        return default

    try:
        parsed = int(value)
    except (TypeError, ValueError):
        raise ValueError(f"Invalid {field_name}. It must be an integer.")

    if parsed < 1:
        raise ValueError(f"Invalid {field_name}. It must be >= 1.")

    if max_value is not None and parsed > max_value:
        raise ValueError(
            f"Invalid {field_name}. It must be <= {max_value}."
        )

    return parsed


def normalize_gene_list(raw_genes) -> list[str]:
    if not isinstance(raw_genes, list):
        raise ValueError("Field 'genes' must be a list.")

    return [
        str(gene).strip()
        for gene in raw_genes
        if str(gene).strip()
    ]


def validate_selected_genes(genes: list[str]) -> Response | None:
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

    return None


def validate_genes_exist(
    genes: list[str],
    available_genes: set[str],
) -> Response | None:
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

    return None


def read_parquet_page(
    file_path,
    columns: list[str],
    page: int,
    page_size: int,
) -> tuple[pd.DataFrame, int]:
    parquet_file = pq.ParquetFile(file_path)

    total_count = parquet_file.metadata.num_rows
    row_start = (page - 1) * page_size

    if row_start >= total_count and total_count > 0:
        raise ExpressionPathError("Page is out of range.")

    remaining_skip = row_start
    remaining_take = page_size
    batches = []

    for batch in parquet_file.iter_batches(
            batch_size=page_size,
            columns=columns,
            use_pandas_metadata=True,
    ):
        batch_rows = batch.num_rows

        if remaining_skip >= batch_rows:
            remaining_skip -= batch_rows
            continue

        if remaining_skip > 0:
            batch = batch.slice(remaining_skip)
            batch_rows = batch.num_rows
            remaining_skip = 0

        take_rows = min(batch_rows, remaining_take)
        batches.append(batch.slice(0, take_rows))
        remaining_take -= take_rows

        if remaining_take <= 0:
            break

    if not batches:
        return pd.DataFrame(columns=columns), total_count

    table = pa.Table.from_batches(batches)

    # 如果 parquet 是 pandas 保存的，to_pandas 通常会恢复 index
    df = table.to_pandas()

    return df, total_count


def get_large_expression_id_column(expression_mode: str) -> str:
    normalized_mode = str(expression_mode or "").strip().lower()

    if normalized_mode == "tisch2":
        return "cell_id"

    if normalized_mode == "sctml":
        return "spot_id"

    return "observation_id"


class DatasetLargeExpressionDataView(APIView):
    def post(self, request):
        dataset = request.data.get("dataset")
        expression_type = request.data.get("expression_type")
        raw_genes = request.data.get("genes", [])

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

        try:
            page = parse_positive_int(
                request.data.get("page"),
                default=DEFAULT_LARGE_EXPRESSION_PAGE,
                field_name="page",
            )
            page_size = parse_positive_int(
                request.data.get("page_size"),
                default=DEFAULT_LARGE_EXPRESSION_PAGE_SIZE,
                field_name="page_size",
                max_value=MAX_LARGE_EXPRESSION_PAGE_SIZE,
            )
        except ValueError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )

        try:
            genes = normalize_gene_list(raw_genes)
        except ValueError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )

        gene_error = validate_selected_genes(genes)
        if gene_error is not None:
            return gene_error

        try:
            metadata = DatasetMetadata.objects.get(dataset=dataset)
        except DatasetMetadata.DoesNotExist:
            return Response(
                {"detail": f"Dataset metadata not found: {dataset}."},
                status=status.HTTP_404_NOT_FOUND,
            )

        rna_type = metadata.gene_bio_type
        expression_mode = get_expression_mode_from_metadata(metadata)

        if str(expression_mode or "").strip().lower() not in LARGE_EXPRESSION_MODES:
            return Response(
                {
                    "detail": (
                        "DatasetLargeExpressionDataView is only available "
                        "for TISCH2 and scTML datasets."
                    )
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

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

        if file_format != "parquet":
            return Response(
                {
                    "detail": (
                        "Large expression data endpoint only supports "
                        "parquet files."
                    )
                },
                status=status.HTTP_400_BAD_REQUEST,
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

        available_genes = set(all_columns)

        gene_exist_error = validate_genes_exist(
            genes=genes,
            available_genes=available_genes,
        )
        if gene_exist_error is not None:
            return gene_exist_error

        try:
            df, total_count = read_parquet_page(
                file_path=file_path,
                columns=genes,
                page=page,
                page_size=page_size,
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
                        f"Failed to read large expression parquet data: {e}"
                    )
                },
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

        id_column = get_large_expression_id_column(expression_mode)

        df = df.reset_index()

        first_column = df.columns[0]

        if first_column != id_column:
            df = df.rename(columns={first_column: id_column})

        df = df.fillna("")

        return Response(
            {
                "dataset": metadata.dataset,
                "rna_type": rna_type,
                "expression_mode": expression_mode,
                "expression_type": expression_type,
                "file_format": file_format,
                "count": total_count,
                "page": page,
                "page_size": page_size,
                "id_column": id_column,
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


class TISCH2DEGClusterPlotView(APIView):
    REQUIRED_COLUMNS = {
        "Cluster",
        "Celltype (malignancy)",
        "Celltype (major-lineage)",
        "Celltype (minor-lineage)",
        "Gene",
        "log2FC",
        "Percentage (%)",
        "Adjusted p-value",
    }

    DEFAULT_EXPRESSION_TYPE = "exp"
    DEFAULT_PADJ_CUTOFF = 0.05
    DEFAULT_PANEL_GAP = 0
    MIN_X_ABS_LIMIT = 1.0

    def parse_float_param(
        self,
        request,
        name: str,
        default: float,
    ) -> float:
        value = request.query_params.get(name)

        if value in [None, ""]:
            return default

        try:
            parsed = float(value)
        except (TypeError, ValueError):
            raise ValueError(f"Invalid {name}. It must be a number.")

        return parsed

    def cluster_sort_key(self, value):
        text = str(value)

        try:
            return 0, int(float(text))
        except (TypeError, ValueError):
            return 1, text

    def get(self, request, dataset):
        expression_type = request.query_params.get(
            "expression_type",
            self.DEFAULT_EXPRESSION_TYPE,
        )

        padj_cutoff = self.DEFAULT_PADJ_CUTOFF

        try:
            panel_gap = self.parse_float_param(
                request=request,
                name="panel_gap",
                default=self.DEFAULT_PANEL_GAP,
            )
        except ValueError as e:
            return Response(
                {"detail": str(e)},
                status=status.HTTP_400_BAD_REQUEST,
            )

        if panel_gap < 0:
            return Response(
                {"detail": "Invalid panel_gap. It must be >= 0."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        try:
            metadata = DatasetMetadata.objects.get(dataset=dataset)
        except DatasetMetadata.DoesNotExist:
            return Response(
                {"detail": f"Dataset metadata not found: {dataset}."},
                status=status.HTTP_404_NOT_FOUND,
            )

        expression_mode = get_expression_mode_from_metadata(metadata)

        if str(expression_mode or "").strip().lower() != "tisch2":
            return Response(
                {
                    "detail": (
                        "TISCH2 DEG cluster plot is only available for "
                        "TISCH2 cell datasets."
                    )
                },
                status=status.HTTP_400_BAD_REQUEST,
            )

        rna_type = metadata.gene_bio_type

        try:
            deg_file = validate_tisch2_deg_file(
                dataset=metadata.dataset,
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
            df = pd.read_csv(deg_file, low_memory=False)
        except Exception as e:
            return Response(
                {"detail": f"Failed to read TISCH2 DEG file: {e}"},
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

        raw_count = int(df.shape[0])

        df = df[
            [
                "Cluster",
                "Celltype (malignancy)",
                "Celltype (major-lineage)",
                "Celltype (minor-lineage)",
                "Gene",
                "log2FC",
                "Percentage (%)",
                "Adjusted p-value",
            ]
        ].copy()

        df = df.rename(
            columns={
                "Cluster": "cluster",
                "Celltype (malignancy)": "celltype_malignancy",
                "Celltype (major-lineage)": "celltype_major_lineage",
                "Celltype (minor-lineage)": "celltype_minor_lineage",
                "Gene": "gene",
                "Percentage (%)": "percentage",
                "Adjusted p-value": "adjusted_p",
            }
        )

        df = df.replace([np.inf, -np.inf], np.nan)

        for col in [
            "cluster",
            "celltype_malignancy",
            "celltype_major_lineage",
            "celltype_minor_lineage",
            "gene",
        ]:
            df[col] = df[col].astype(str).str.strip()

        df["log2FC"] = pd.to_numeric(
            df["log2FC"],
            errors="coerce",
        )
        df["percentage"] = pd.to_numeric(
            df["percentage"],
            errors="coerce",
        )
        df["adjusted_p"] = pd.to_numeric(
            df["adjusted_p"],
            errors="coerce",
        )

        before_dropna_count = int(df.shape[0])

        df = df.dropna(
            subset=[
                "cluster",
                "gene",
                "log2FC",
                "percentage",
                "adjusted_p",
            ]
        )

        dropna_count = before_dropna_count - int(df.shape[0])

        zero_p_count = int((df["adjusted_p"] == 0).sum())
        invalid_p_count = int(
            ((df["adjusted_p"] < 0) | (df["adjusted_p"] > 1)).sum()
        )

        df = df[
            (df["adjusted_p"] > 0)
            & (df["adjusted_p"] <= 1)
        ].copy()

        if df.empty:
            return Response(
                {
                    "dataset": metadata.dataset,
                    "rna_type": rna_type,
                    "expression_mode": expression_mode,
                    "expression_type": expression_type,
                    "summary": {
                        "raw_count": raw_count,
                        "cleaned_count": 0,
                        "dropna_count": dropna_count,
                        "zero_p_dropped_count": zero_p_count,
                        "invalid_p_dropped_count": invalid_p_count,
                        "up_count": 0,
                        "down_count": 0,
                        "not_count": 0,
                    },
                    "thresholds": {
                        "padj_cutoff": padj_cutoff,
                        "neg_log10_padj_cutoff": float(-np.log10(padj_cutoff)),
                    },
                    "clusters": [],
                    "points": [],
                },
                status=status.HTTP_200_OK,
            )

        df["neg_log10_adjusted_p"] = -np.log10(df["adjusted_p"])

        df["regulation"] = "Not"
        df.loc[
            (df["adjusted_p"] <= padj_cutoff) & (df["log2FC"] > 0),
            "regulation",
        ] = "Up"
        df.loc[
            (df["adjusted_p"] <= padj_cutoff) & (df["log2FC"] < 0),
            "regulation",
        ] = "Down"

        df["abs_log2FC"] = df["log2FC"].abs()

        clusters = sorted(
            df["cluster"].unique().tolist(),
            key=self.cluster_sort_key,
        )

        cluster_to_index = {
            cluster: index
            for index, cluster in enumerate(clusters)
        }

        # Cluster 内查看优先：
        # 每个 cluster 内部使用自己的 log2FC 范围做归一化。
        # 前端显示位置用 scaled_log2FC，tooltip 仍然显示原始 log2FC。
        panel_half_width = 1.0
        panel_width = 2 * panel_half_width
        panel_span = panel_width + panel_gap

        cluster_panels = []
        zero_x_map = {}
        cluster_x_abs_max_map = {}

        for cluster in clusters:
            panel_index = cluster_to_index[cluster]
            sub_df = df[df["cluster"] == cluster]

            local_x_abs_max = float(
                max(
                    sub_df["log2FC"].abs().max(),
                    self.MIN_X_ABS_LIMIT,
                )
            )

            panel_start = panel_index * panel_span
            zero_x = panel_start + panel_half_width
            panel_end = panel_start + panel_width

            zero_x_map[cluster] = zero_x
            cluster_x_abs_max_map[cluster] = local_x_abs_max

            up_count = int((sub_df["regulation"] == "Up").sum())
            down_count = int((sub_df["regulation"] == "Down").sum())
            not_count = int((sub_df["regulation"] == "Not").sum())

            first_row = sub_df.iloc[0]

            cluster_panels.append(
                {
                    "cluster": cluster,
                    "panel_index": int(panel_index),
                    "panel_start": float(panel_start),
                    "panel_end": float(panel_end),
                    "panel_center": float(zero_x),
                    "zero_x": float(zero_x),
                    "count": int(sub_df.shape[0]),
                    "up_count": up_count,
                    "down_count": down_count,
                    "not_count": not_count,
                    "celltype_malignancy": first_row["celltype_malignancy"],
                    "celltype_major_lineage": (
                        first_row["celltype_major_lineage"]
                    ),
                    "celltype_minor_lineage": (
                        first_row["celltype_minor_lineage"]
                    ),
                    "x_min": float(sub_df["log2FC"].min()),
                    "x_max": float(sub_df["log2FC"].max()),
                    "x_abs_max": float(local_x_abs_max),
                    "y_max": float(sub_df["neg_log10_adjusted_p"].max()),
                }
            )

        df["panel_index"] = df["cluster"].map(cluster_to_index)
        df["zero_x"] = df["cluster"].map(zero_x_map)
        df["cluster_x_abs_max"] = df["cluster"].map(cluster_x_abs_max_map)

        df["scaled_log2FC"] = df["log2FC"] / df["cluster_x_abs_max"]
        df["plot_x"] = df["zero_x"] + df["scaled_log2FC"] * panel_half_width
        df["plot_y"] = df["neg_log10_adjusted_p"]

        cleaned_count = int(df.shape[0])
        up_count = int((df["regulation"] == "Up").sum())
        down_count = int((df["regulation"] == "Down").sum())
        not_count = int((df["regulation"] == "Not").sum())

        points = [
            {
                "cluster": row["cluster"],
                "panel_index": int(row["panel_index"]),
                "gene": row["gene"],
                "log2FC": float(row["log2FC"]),
                "scaled_log2FC": float(row["scaled_log2FC"]),
                "cluster_x_abs_max": float(row["cluster_x_abs_max"]),
                "percentage": float(row["percentage"]),
                "adjusted_p": float(row["adjusted_p"]),
                "neg_log10_adjusted_p": float(
                    row["neg_log10_adjusted_p"]
                ),
                "plot_x": float(row["plot_x"]),
                "plot_y": float(row["plot_y"]),
                "regulation": row["regulation"],
                "celltype_malignancy": row["celltype_malignancy"],
                "celltype_major_lineage": row[
                    "celltype_major_lineage"
                ],
                "celltype_minor_lineage": row[
                    "celltype_minor_lineage"
                ],
            }
            for _, row in df.iterrows()
        ]

        return Response(
            {
                "dataset": metadata.dataset,
                "rna_type": rna_type,
                "expression_mode": expression_mode,
                "expression_type": expression_type,
                "summary": {
                    "raw_count": raw_count,
                    "cleaned_count": cleaned_count,
                    "dropna_count": dropna_count,
                    "zero_p_dropped_count": zero_p_count,
                    "invalid_p_dropped_count": invalid_p_count,
                    "up_count": up_count,
                    "down_count": down_count,
                    "not_count": not_count,
                    "cluster_count": len(clusters),
                },
                "thresholds": {
                    "padj_cutoff": padj_cutoff,
                    "neg_log10_padj_cutoff": float(-np.log10(padj_cutoff)),
                    "log2fc_center": 0.0,
                    "x_scale_mode": "per_cluster_normalized",
                    "panel_half_width": float(panel_half_width),
                    "panel_width": float(panel_width),
                    "panel_gap": float(panel_gap),
                    "panel_span": float(panel_span),
                },
                "clusters": cluster_panels,
                "points": points,
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


class DatasetDownloadView(APIView):
    def get(self, request):
        dataset = request.query_params.get("dataset")

        result = prepare_dataset_download(dataset)

        return FileResponse(
            open(result.archive_path, "rb"),
            as_attachment=True,
            filename=result.archive_name,
        )


class DatasetAnnotationDownloadView(APIView):
    def get(self, request):
        dataset = request.query_params.get("dataset")

        result = prepare_tcga_annotation_download(dataset)

        return FileResponse(
            open(result.archive_path, "rb"),
            as_attachment=True,
            filename=result.archive_name,
        )


class TIMEDBAnnotationDownloadView(APIView):
    def get(self, request):
        dataset = request.query_params.get("dataset")

        result = prepare_timedb_annotation_download(dataset)

        return FileResponse(
            open(result.archive_path, "rb"),
            as_attachment=True,
            filename=result.archive_name,
        )
