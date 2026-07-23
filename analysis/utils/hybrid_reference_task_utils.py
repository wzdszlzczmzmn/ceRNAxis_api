from collections import Counter
from pathlib import Path
import csv
import json

import pandas as pd
import pyarrow.parquet as pq
from pyarrow.lib import ArrowInvalid

from analysis.utils.paired_cohort_task_utils import (
    validate_safe_name,
)

from django.conf import settings

from analysis.utils.workflow_detail_utils.deg_pathway_utils import DEG_PATHWAY_REQUIRED_COLUMNS, \
    validate_deg_pathway_dataframe_columns, DEGPathwayInputError, build_deg_pathway_data_from_dataframe
from analysis.utils.workflow_detail_utils.survival_km_utils import SURVIVAL_REQUIRED_COLUMNS, DEFAULT_SURVIVAL_GROUPS, \
    validate_survival_dataframe_columns, SurvivalKMInputError, build_survival_km_data_from_dataframe

HYBRID_REFERENCE_ALLOWED_FILE_FIELDS = [
    "mrna_file",
    "meta_file",
]

HYBRID_REFERENCE_INPUT_FILENAME_MAP = {
    "mrna_file": "mrna.csv",
    "meta_file": "meta.csv",
}

HYBRID_REFERENCE_VALID_DEG_METHODS = [
    "limma",
    "deseq2",
]

HYBRID_REFERENCE_VALID_LNCRNA_TYPES = [
    "log2count",
    "log2fpkm",
    "log2fpkmuq",
    "log2tpm",
]

HYBRID_REFERENCE_VALID_TCGA_TYPES = [
    "TCGA_ACC",
    "TCGA_BLCA",
    "TCGA_BRCA",
    "TCGA_CESC",
    "TCGA_CHOL",
    "TCGA_COAD",
    "TCGA_DLBC",
    "TCGA_ESCA",
    "TCGA_GBM",
    "TCGA_HNSC",
    "TCGA_KICH",
    "TCGA_KIRC",
    "TCGA_KIRP",
    "TCGA_LAML",
    "TCGA_LGG",
    "TCGA_LIHC",
    "TCGA_LUAD",
    "TCGA_LUSC",
    "TCGA_MESO",
    "TCGA_OV",
    "TCGA_PAAD",
    "TCGA_PCPG",
    "TCGA_PRAD",
    "TCGA_READ",
    "TCGA_SARC",
    "TCGA_SKCM",
    "TCGA_STAD",
    "TCGA_TGCT",
    "TCGA_THCA",
    "TCGA_THYM",
    "TCGA_UCEC",
    "TCGA_UCS",
    "TCGA_UVM",
]

HYBRID_REFERENCE_EXPR_SAMPLE_COL = "sample_id"
HYBRID_REFERENCE_META_SAMPLE_COL = "sample_id"
HYBRID_REFERENCE_GROUP_COL = "c_group"
HYBRID_REFERENCE_CASE_LABEL = "case"
HYBRID_REFERENCE_CONTROL_LABEL = "control"

EXPRESSION_REQUIRED_COLUMNS = ["sample_id"]
META_REQUIRED_COLUMNS = ["sample_id", "c_group"]
META_REQUIRED_GROUPS = {"case", "control"}


HYBRID_REFERENCE_SURVIVAL_FILENAME_SUFFIX = "_survival_analysis.csv"

HYBRID_REFERENCE_SURVIVAL_REQUIRED_COLUMNS = SURVIVAL_REQUIRED_COLUMNS

HYBRID_REFERENCE_SURVIVAL_GROUPS = DEFAULT_SURVIVAL_GROUPS


HYBRID_REFERENCE_MRNA_GSEA_FILENAME_SUFFIX = "_mRNA_gsea.csv"

HYBRID_REFERENCE_MRNA_GSEA_REQUIRED_COLUMNS = DEG_PATHWAY_REQUIRED_COLUMNS


SCST_HYBRID_REFERENCE_ALLOWED_FILE_FIELDS = [
    "exp_file",
    "meta_file",
]

SCST_HYBRID_REFERENCE_INPUT_FILENAME_MAP = {
    "exp_file": "expression.parquet",
    "meta_file": "meta.csv",
}

SCST_HYBRID_REFERENCE_ALLOWED_FILE_SUFFIXES = {
    "exp_file": {".parquet"},
    "meta_file": {".csv"},
}

SCST_HYBRID_REFERENCE_VALID_DATA_TYPES = [
    "sc",
    "st",
]

SCST_HYBRID_REFERENCE_ID_COLUMN_MAP = {
    "sc": "cell_id",
    "st": "spot_id",
}


class HybridReferenceTaskInputError(ValueError):
    pass


class HybridReferenceTaskPathError(ValueError):
    pass


class SCSTHybridReferenceTaskInputError(HybridReferenceTaskInputError):
    pass


class SCSTHybridReferenceTaskPathError(HybridReferenceTaskPathError):
    pass


def get_hybrid_reference_task_input_dir(task) -> Path:
    return Path(task.get_input_dir_absolute_path()).resolve()


def get_hybrid_reference_task_output_dir(task) -> Path:
    return Path(task.get_output_dir_absolute_path()).resolve()


def get_hybrid_reference_input_file_path(task, field_name: str) -> Path:
    if field_name not in HYBRID_REFERENCE_INPUT_FILENAME_MAP:
        raise HybridReferenceTaskPathError(
            f"Invalid input file field: {field_name}."
        )

    input_dir = get_hybrid_reference_task_input_dir(task)

    file_path = (
        input_dir / HYBRID_REFERENCE_INPUT_FILENAME_MAP[field_name]
    ).resolve()

    if not str(file_path).startswith(str(input_dir)):
        raise HybridReferenceTaskPathError(
            "Invalid hybrid reference input file path."
        )

    return file_path


def prepare_hybrid_reference_workspace(task) -> dict:
    input_dir = get_hybrid_reference_task_input_dir(task)
    output_dir = get_hybrid_reference_task_output_dir(task)

    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    return {
        "input_dir": input_dir,
        "output_dir": output_dir,
    }


def save_hybrid_reference_uploaded_input_files(task, files) -> dict:
    input_dir = get_hybrid_reference_task_input_dir(task)
    input_dir.mkdir(parents=True, exist_ok=True)

    saved_files = {}

    for field_name in HYBRID_REFERENCE_ALLOWED_FILE_FIELDS:
        if field_name not in files:
            raise HybridReferenceTaskInputError(
                f"Missing uploaded file: {field_name}."
            )

        uploaded_file = files[field_name]
        file_path = get_hybrid_reference_input_file_path(task, field_name)

        with file_path.open("wb+") as destination:
            for chunk in uploaded_file.chunks():
                destination.write(chunk)

        saved_files[field_name] = file_path.name

    return saved_files


def validate_hybrid_reference_input_files(task) -> dict:
    validated_files = {}

    for field_name in HYBRID_REFERENCE_ALLOWED_FILE_FIELDS:
        file_path = get_hybrid_reference_input_file_path(task, field_name)

        if not file_path.exists() or not file_path.is_file():
            raise FileNotFoundError(
                f"Hybrid reference input file not found: {file_path}"
            )

        validated_files[field_name] = file_path

    return validated_files


def read_csv_header(file_path: Path) -> list[str]:
    try:
        with file_path.open("r", newline="", encoding="utf-8-sig") as f:
            reader = csv.reader(f)
            header = next(reader, None)
    except UnicodeDecodeError:
        raise HybridReferenceTaskInputError(
            f"File must be UTF-8 encoded: {file_path.name}."
        )
    except csv.Error as e:
        raise HybridReferenceTaskInputError(
            f"Invalid CSV file: {file_path.name}. {str(e)}"
        )

    if not header:
        raise HybridReferenceTaskInputError(
            f"CSV file is empty or missing header: {file_path.name}."
        )

    return [str(col).strip() for col in header]


def validate_required_columns(
    file_path: Path,
    required_columns: list[str],
    file_label: str,
) -> None:
    header = read_csv_header(file_path)

    missing_columns = [
        col for col in required_columns
        if col not in header
    ]

    if missing_columns:
        raise HybridReferenceTaskInputError(
            f"{file_label} is missing required column(s): "
            f"{', '.join(missing_columns)}."
        )


def validate_expression_file_columns(file_path: Path, file_label: str) -> None:
    validate_required_columns(
        file_path=file_path,
        required_columns=EXPRESSION_REQUIRED_COLUMNS,
        file_label=file_label,
    )


def validate_meta_file_columns_and_groups(file_path: Path) -> None:
    try:
        with file_path.open("r", newline="", encoding="utf-8-sig") as f:
            reader = csv.reader(f)
            header = next(reader, None)

            if not header:
                raise HybridReferenceTaskInputError(
                    "Meta file is empty or missing header."
                )

            normalized_header = [str(col).strip() for col in header]

            missing_columns = [
                col for col in META_REQUIRED_COLUMNS
                if col not in normalized_header
            ]

            if missing_columns:
                raise HybridReferenceTaskInputError(
                    "Meta file is missing required column(s): "
                    f"{', '.join(missing_columns)}."
                )

            dict_reader = csv.DictReader(
                f,
                fieldnames=normalized_header,
            )

            observed_required_groups = set()
            row_count = 0

            for row in dict_reader:
                row_count += 1

                group = str(row.get(HYBRID_REFERENCE_GROUP_COL, "")).strip()

                if group in META_REQUIRED_GROUPS:
                    observed_required_groups.add(group)

    except UnicodeDecodeError:
        raise HybridReferenceTaskInputError(
            f"File must be UTF-8 encoded: {file_path.name}."
        )
    except csv.Error as e:
        raise HybridReferenceTaskInputError(
            f"Invalid CSV file: {file_path.name}. {str(e)}"
        )

    if row_count == 0:
        raise HybridReferenceTaskInputError(
            "Meta file has no data rows."
        )

    missing_groups = META_REQUIRED_GROUPS - observed_required_groups

    if missing_groups:
        raise HybridReferenceTaskInputError(
            "Meta file column 'c_group' must contain at least one case "
            "and at least one control sample. "
            f"Missing group(s): {', '.join(sorted(missing_groups))}."
        )


def validate_hybrid_reference_file_contents(task) -> dict:
    input_files = validate_hybrid_reference_input_files(task)

    validate_expression_file_columns(
        file_path=input_files["mrna_file"],
        file_label="mRNA expression file",
    )

    validate_meta_file_columns_and_groups(
        file_path=input_files["meta_file"],
    )

    return input_files


def validate_hybrid_reference_task_params(
    *,
    task_name: str,
    tcga_type: str,
    lncrna_type: str,
    deg_method: str,
    use_padj: bool = True,
) -> None:
    validate_safe_name(task_name, "task_name")

    if tcga_type not in HYBRID_REFERENCE_VALID_TCGA_TYPES:
        raise HybridReferenceTaskInputError(
            "Invalid field: tcga_type. Allowed values are: "
            f"{', '.join(HYBRID_REFERENCE_VALID_TCGA_TYPES)}."
        )

    if lncrna_type not in HYBRID_REFERENCE_VALID_LNCRNA_TYPES:
        raise HybridReferenceTaskInputError(
            "Invalid field: lncrna_type. Allowed values are: "
            f"{', '.join(HYBRID_REFERENCE_VALID_LNCRNA_TYPES)}."
        )

    if deg_method not in HYBRID_REFERENCE_VALID_DEG_METHODS:
        raise HybridReferenceTaskInputError(
            "Invalid field: deg_method. Allowed values are: "
            f"{', '.join(HYBRID_REFERENCE_VALID_DEG_METHODS)}."
        )

    if not isinstance(use_padj, bool):
        raise HybridReferenceTaskInputError(
            "Invalid field: use_padj. Allowed values are TRUE or FALSE."
        )


def get_hybrid_reference_tcga_expression_filename(
    *,
    tcga_type: str,
    lncrna_type: str,
    rna_type: str,
) -> str:
    tcga_type = str(tcga_type or "").strip()
    lncrna_type = str(lncrna_type or "").strip()
    rna_type = str(rna_type or "").strip()

    validate_safe_name(tcga_type, "tcga_type")
    validate_safe_name(lncrna_type, "lncrna_type")
    validate_safe_name(rna_type, "rna_type")

    if rna_type == "mRNA":
        return f"{tcga_type}_mRNA_{lncrna_type}_exp.csv"

    if rna_type == "miRNA":
        return f"{tcga_type}_miRNA_log2rpm_exp.csv"

    if rna_type == "lncRNA":
        return f"{tcga_type}_lncRNA_{lncrna_type}_exp.csv"

    if rna_type == "circRNA":
        return f"{tcga_type}_circRNA_count_exp.csv"

    raise HybridReferenceTaskPathError(
        f"Unsupported Hybrid Reference RNA type: {rna_type}."
    )


def get_hybrid_reference_tcga_expression_file_path(
    *,
    task,
    rna_type: str,
) -> Path:
    base_dir = Path(settings.TCGA_DATASET_BASE_DIR).resolve()

    filename = get_hybrid_reference_tcga_expression_filename(
        tcga_type=task.tcga_type,
        lncrna_type=task.lncrna_type,
        rna_type=rna_type,
    )

    file_path = (base_dir / filename).resolve()

    if not str(file_path).startswith(str(base_dir)):
        raise HybridReferenceTaskPathError(
            "Invalid Hybrid Reference TCGA expression file path."
        )

    return file_path


def get_hybrid_reference_survival_file_path(task) -> Path:
    task_name = str(task.task_name).strip()

    validate_safe_name(task_name, "task_name")

    output_dir = get_hybrid_reference_task_output_dir(task)

    file_path = (
        output_dir / f"{task_name}{HYBRID_REFERENCE_SURVIVAL_FILENAME_SUFFIX}"
    ).resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise HybridReferenceTaskPathError(
            "Invalid hybrid reference survival analysis file path."
        )

    return file_path


def validate_hybrid_reference_survival_file(task) -> Path:
    file_path = get_hybrid_reference_survival_file_path(task)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Hybrid reference survival analysis file not found: {file_path.name}"
        )

    return file_path


def read_hybrid_reference_survival_file(task) -> tuple[Path, pd.DataFrame]:
    file_path = validate_hybrid_reference_survival_file(task)

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise HybridReferenceTaskInputError(
            f"Failed to read hybrid reference survival analysis file: {str(e)}"
        )

    try:
        validate_survival_dataframe_columns(
            df=df,
            required_columns=HYBRID_REFERENCE_SURVIVAL_REQUIRED_COLUMNS,
        )
    except SurvivalKMInputError as e:
        raise HybridReferenceTaskInputError(
            "Hybrid reference survival analysis file is invalid. "
            f"{str(e)}"
        )

    return file_path, df


def build_hybrid_reference_survival_km_data(
    task,
    title: str = "TCGA-based ceRNA axis survival analysis",
) -> dict:
    survival_file, df = read_hybrid_reference_survival_file(task)

    try:
        return build_survival_km_data_from_dataframe(
            task=task,
            survival_file_name=survival_file.name,
            df=df,
            title=title,
            valid_groups=HYBRID_REFERENCE_SURVIVAL_GROUPS,
        )
    except SurvivalKMInputError as e:
        raise HybridReferenceTaskInputError(str(e))


def get_hybrid_reference_mrna_gsea_file_path(task) -> Path:
    task_name = str(task.task_name).strip()

    validate_safe_name(task_name, "task_name")

    output_dir = get_hybrid_reference_task_output_dir(task)

    file_path = (
        output_dir / f"{task_name}{HYBRID_REFERENCE_MRNA_GSEA_FILENAME_SUFFIX}"
    ).resolve()

    if not str(file_path).startswith(str(output_dir)):
        raise HybridReferenceTaskPathError(
            "Invalid hybrid reference mRNA GSEA result file path."
        )

    return file_path


def validate_hybrid_reference_mrna_gsea_file(task) -> Path:
    file_path = get_hybrid_reference_mrna_gsea_file_path(task)

    if not file_path.exists() or not file_path.is_file():
        raise FileNotFoundError(
            f"Hybrid reference mRNA GSEA result file not found: {file_path.name}"
        )

    return file_path


def read_hybrid_reference_mrna_gsea_file(task) -> tuple[Path, pd.DataFrame]:
    file_path = validate_hybrid_reference_mrna_gsea_file(task)

    try:
        df = pd.read_csv(file_path)
    except Exception as e:
        raise HybridReferenceTaskInputError(
            f"Failed to read hybrid reference mRNA GSEA result file: {str(e)}"
        )

    try:
        validate_deg_pathway_dataframe_columns(
            df=df,
            required_columns=HYBRID_REFERENCE_MRNA_GSEA_REQUIRED_COLUMNS,
        )
    except DEGPathwayInputError as e:
        raise HybridReferenceTaskInputError(
            "Hybrid reference mRNA GSEA result file is invalid. "
            f"{str(e)}"
        )

    return file_path, df


def build_hybrid_reference_deg_pathway_data(
    task,
    title: str = "TCGA-based DEG Pathway Enrichment",
) -> dict:
    gsea_file, df = read_hybrid_reference_mrna_gsea_file(task)

    try:
        return build_deg_pathway_data_from_dataframe(
            task=task,
            gsea_file_name=gsea_file.name,
            df=df,
            title=title,
        )
    except DEGPathwayInputError as e:
        raise HybridReferenceTaskInputError(str(e))


def get_scst_hybrid_reference_task_input_dir(task) -> Path:
    return Path(task.get_input_dir_absolute_path()).resolve()


def get_scst_hybrid_reference_task_output_dir(task) -> Path:
    return Path(task.get_output_dir_absolute_path()).resolve()


def get_scst_hybrid_reference_input_file_path(
    task,
    field_name: str,
) -> Path:
    if field_name not in SCST_HYBRID_REFERENCE_INPUT_FILENAME_MAP:
        raise SCSTHybridReferenceTaskPathError(
            f"Invalid SC/ST input file field: {field_name}."
        )

    input_dir = get_scst_hybrid_reference_task_input_dir(task)

    file_path = (
        input_dir / SCST_HYBRID_REFERENCE_INPUT_FILENAME_MAP[field_name]
    ).resolve()

    if file_path.parent != input_dir:
        raise SCSTHybridReferenceTaskPathError(
            "Invalid SC/ST Hybrid Reference input file path."
        )

    return file_path


def prepare_scst_hybrid_reference_workspace(task) -> dict:
    input_dir = get_scst_hybrid_reference_task_input_dir(task)
    output_dir = get_scst_hybrid_reference_task_output_dir(task)

    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    return {
        "input_dir": input_dir,
        "output_dir": output_dir,
    }


def save_scst_hybrid_reference_uploaded_input_files(
    task,
    files,
) -> dict:
    input_dir = get_scst_hybrid_reference_task_input_dir(task)
    input_dir.mkdir(parents=True, exist_ok=True)

    saved_files = {}

    for field_name in SCST_HYBRID_REFERENCE_ALLOWED_FILE_FIELDS:
        uploaded_file = files.get(field_name)

        if uploaded_file is None:
            raise SCSTHybridReferenceTaskInputError(
                f"Missing uploaded file: {field_name}."
            )

        file_path = get_scst_hybrid_reference_input_file_path(
            task=task,
            field_name=field_name,
        )

        with file_path.open("wb+") as destination:
            for chunk in uploaded_file.chunks():
                destination.write(chunk)

        saved_files[field_name] = file_path.name

    return saved_files


def validate_scst_hybrid_reference_input_files(task) -> dict:
    validated_files = {}

    for field_name in SCST_HYBRID_REFERENCE_ALLOWED_FILE_FIELDS:
        file_path = get_scst_hybrid_reference_input_file_path(
            task=task,
            field_name=field_name,
        )

        if not file_path.exists() or not file_path.is_file():
            raise FileNotFoundError(
                f"SC/ST Hybrid Reference input file not found: "
                f"{file_path.name}"
            )

        validated_files[field_name] = file_path

    return validated_files


def read_parquet_columns(file_path: Path) -> list[str]:
    try:
        parquet_file = pq.ParquetFile(file_path)
        columns = parquet_file.schema_arrow.names

    except (
        ArrowInvalid,
        OSError,
        ValueError,
    ) as exc:
        raise SCSTHybridReferenceTaskInputError(
            f"Invalid Parquet file: {file_path.name}. {exc}"
        )

    if not columns:
        raise SCSTHybridReferenceTaskInputError(
            f"Parquet file has no columns: {file_path.name}."
        )

    return [
        str(column).strip()
        for column in columns
    ]


def get_parquet_pandas_index_columns(
    parquet_file: pq.ParquetFile,
) -> list[str]:
    """
    从 Parquet 的 pandas metadata 中提取命名索引字段。

    例如：
        dataframe.index.name = "cell_id"
        dataframe.to_parquet(..., index=True)

    此时 pandas metadata 的 index_columns 通常为：
        ["cell_id"]

    RangeIndex 可能以 dict 形式记录，此处不会将其视为命名 ID 索引。
    """
    schema_metadata = parquet_file.schema_arrow.metadata or {}
    raw_pandas_metadata = schema_metadata.get(b"pandas")

    if not raw_pandas_metadata:
        return []

    try:
        pandas_metadata = json.loads(
            raw_pandas_metadata.decode("utf-8")
        )
    except (
        UnicodeDecodeError,
        json.JSONDecodeError,
        TypeError,
    ):
        return []

    raw_index_columns = pandas_metadata.get(
        "index_columns",
        [],
    )

    return [
        str(index_column).strip()
        for index_column in raw_index_columns
        if isinstance(index_column, str)
        and str(index_column).strip()
    ]


def validate_scst_expression_parquet_schema(
    file_path: Path,
    expected_id_column: str,
) -> dict:
    try:
        parquet_file = pq.ParquetFile(file_path)
        columns = [
            str(column).strip()
            for column in parquet_file.schema_arrow.names
        ]

    except (
        ArrowInvalid,
        OSError,
        ValueError,
    ) as exc:
        raise SCSTHybridReferenceTaskInputError(
            f"Invalid Parquet file: {file_path.name}. {exc}"
        )

    if not columns:
        raise SCSTHybridReferenceTaskInputError(
            f"Parquet file has no columns: {file_path.name}."
        )

    column_counts = Counter(columns)

    duplicate_columns = sorted(
        column
        for column, count in column_counts.items()
        if count > 1
    )

    if duplicate_columns:
        raise SCSTHybridReferenceTaskInputError(
            "exp_file contains duplicate column name(s): "
            f"{', '.join(duplicate_columns[:10])}."
        )

    pandas_index_columns = get_parquet_pandas_index_columns(
        parquet_file
    )

    id_is_first_regular_column = (
        columns[0] == expected_id_column
    )

    id_is_pandas_index = (
        expected_id_column in pandas_index_columns
    )

    id_exists_in_physical_schema = (
        expected_id_column in columns
    )

    if id_is_first_regular_column:
        id_storage = "column"

    elif id_is_pandas_index:
        if not id_exists_in_physical_schema:
            raise SCSTHybridReferenceTaskInputError(
                f"exp_file pandas metadata declares "
                f"'{expected_id_column}' as an index, but the "
                "corresponding physical Parquet field is missing."
            )

        id_storage = "pandas_index"

    else:
        actual_first_column = columns[0]

        raise SCSTHybridReferenceTaskInputError(
            f"exp_file must store '{expected_id_column}' either "
            "as the first regular column or as a named Pandas index. "
            f"The first regular column is '{actual_first_column}', "
            f"and the declared Pandas index column(s) are "
            f"{pandas_index_columns or 'none'}."
        )

    expression_columns = [
        column
        for column in columns
        if column != expected_id_column
    ]

    if not expression_columns:
        raise SCSTHybridReferenceTaskInputError(
            "exp_file must contain at least one gene "
            "expression column."
        )

    return {
        "columns": columns,
        "expression_columns": expression_columns,
        "id_column": expected_id_column,
        "id_storage": id_storage,
        "pandas_index_columns": pandas_index_columns,
    }


def get_scst_expected_id_column(data_type: str) -> str:
    try:
        return SCST_HYBRID_REFERENCE_ID_COLUMN_MAP[data_type]
    except KeyError:
        raise SCSTHybridReferenceTaskInputError(
            "Invalid field: data_type. Allowed values are: sc, st."
        )


def read_scst_parquet_identifier_set(
    file_path: Path,
    id_column: str,
) -> set[str]:
    try:
        parquet_file = pq.ParquetFile(file_path)
    except (
        ArrowInvalid,
        OSError,
        ValueError,
    ) as exc:
        raise SCSTHybridReferenceTaskInputError(
            f"Invalid Parquet file: {file_path.name}. {exc}"
        )

    physical_columns = set(
        parquet_file.schema_arrow.names
    )

    if id_column not in physical_columns:
        raise SCSTHybridReferenceTaskInputError(
            f"Unable to read identifier field '{id_column}' "
            "from exp_file because it is not present in the "
            "physical Parquet schema."
        )

    identifiers = set()
    row_count = 0

    try:
        for batch in parquet_file.iter_batches(
            columns=[id_column],
            batch_size=65536,
        ):
            if batch.num_columns != 1:
                raise SCSTHybridReferenceTaskInputError(
                    f"Failed to read identifier field "
                    f"'{id_column}' from exp_file."
                )

            values = batch.column(0).to_pylist()

            for raw_identifier in values:
                row_count += 1

                if raw_identifier is None:
                    raise SCSTHybridReferenceTaskInputError(
                        f"exp_file contains a null "
                        f"'{id_column}' value at data row "
                        f"{row_count}."
                    )

                identifier = str(raw_identifier).strip()

                if not identifier:
                    raise SCSTHybridReferenceTaskInputError(
                        f"exp_file contains an empty "
                        f"'{id_column}' value at data row "
                        f"{row_count}."
                    )

                if identifier in identifiers:
                    raise SCSTHybridReferenceTaskInputError(
                        f"exp_file contains duplicate "
                        f"'{id_column}' value: {identifier}."
                    )

                identifiers.add(identifier)

    except SCSTHybridReferenceTaskInputError:
        raise

    except (
        ArrowInvalid,
        OSError,
        ValueError,
        KeyError,
    ) as exc:
        raise SCSTHybridReferenceTaskInputError(
            f"Failed to read identifier field "
            f"'{id_column}' from exp_file: {exc}"
        )

    if row_count == 0:
        raise SCSTHybridReferenceTaskInputError(
            "exp_file has no data rows."
        )

    return identifiers


def validate_scst_meta_csv_schema(
    file_path: Path,
    expected_id_column: str,
    group_col: str,
) -> list[str]:
    columns = read_csv_header(file_path)

    actual_first_column = columns[0]

    if actual_first_column != expected_id_column:
        raise SCSTHybridReferenceTaskInputError(
            "The first column of meta_file must be "
            f"'{expected_id_column}', but got "
            f"'{actual_first_column}'."
        )

    if group_col not in columns:
        raise SCSTHybridReferenceTaskInputError(
            f"meta_file is missing group column: "
            f"'{group_col}'."
        )

    if group_col == expected_id_column:
        raise SCSTHybridReferenceTaskInputError(
            f"group_col cannot be the identifier column "
            f"'{expected_id_column}'."
        )

    return columns


def read_scst_meta_identifiers(
    file_path: Path,
    id_column: str,
    group_col: str,
) -> set[str]:
    identifiers = set()
    row_count = 0

    try:
        with file_path.open(
            "r",
            newline="",
            encoding="utf-8-sig",
        ) as file_obj:
            reader = csv.DictReader(file_obj)

            if not reader.fieldnames:
                raise SCSTHybridReferenceTaskInputError(
                    "meta_file is empty or missing header."
                )

            reader.fieldnames = [
                str(column).strip()
                for column in reader.fieldnames
            ]

            for row_number, row in enumerate(reader, start=2):
                row_count += 1

                identifier = str(
                    row.get(id_column, "")
                ).strip()

                group_value = str(
                    row.get(group_col, "")
                ).strip()

                if not identifier:
                    raise SCSTHybridReferenceTaskInputError(
                        f"meta_file contains an empty "
                        f"'{id_column}' value at row "
                        f"{row_number}."
                    )

                if identifier in identifiers:
                    raise SCSTHybridReferenceTaskInputError(
                        f"meta_file contains duplicate "
                        f"'{id_column}' value: "
                        f"{identifier}."
                    )

                if not group_value:
                    raise SCSTHybridReferenceTaskInputError(
                        f"meta_file contains an empty "
                        f"'{group_col}' value at row "
                        f"{row_number}."
                    )

                identifiers.add(identifier)

    except UnicodeDecodeError:
        raise SCSTHybridReferenceTaskInputError(
            f"meta_file must be UTF-8 encoded: "
            f"{file_path.name}."
        )
    except csv.Error as exc:
        raise SCSTHybridReferenceTaskInputError(
            f"Invalid CSV file: {file_path.name}. {exc}"
        )

    if row_count == 0:
        raise SCSTHybridReferenceTaskInputError(
            "meta_file has no data rows."
        )

    return identifiers


def validate_scst_hybrid_reference_file_contents(task) -> dict:
    input_files = validate_scst_hybrid_reference_input_files(task)

    expected_id_column = get_scst_expected_id_column(
        task.data_type
    )

    exp_schema = validate_scst_expression_parquet_schema(
        file_path=input_files["exp_file"],
        expected_id_column=expected_id_column,
    )

    meta_columns = validate_scst_meta_csv_schema(
        file_path=input_files["meta_file"],
        expected_id_column=expected_id_column,
        group_col=task.group_col,
    )

    expression_identifiers = (
        read_scst_parquet_identifier_set(
            file_path=input_files["exp_file"],
            id_column=expected_id_column,
        )
    )

    meta_identifiers = read_scst_meta_identifiers(
        file_path=input_files["meta_file"],
        id_column=expected_id_column,
        group_col=task.group_col,
    )

    missing_in_expression = sorted(
        meta_identifiers - expression_identifiers
    )

    missing_in_meta = sorted(
        expression_identifiers - meta_identifiers
    )

    if missing_in_expression:
        raise SCSTHybridReferenceTaskInputError(
            f"Some '{expected_id_column}' values in "
            "meta_file are not present in exp_file: "
            f"{', '.join(missing_in_expression[:10])}."
        )

    if missing_in_meta:
        raise SCSTHybridReferenceTaskInputError(
            f"Some '{expected_id_column}' values in "
            "exp_file are not present in meta_file: "
            f"{', '.join(missing_in_meta[:10])}."
        )

    return {
        **input_files,
        "id_column": expected_id_column,
        "id_storage": exp_schema["id_storage"],
        "exp_columns": exp_schema["columns"],
        "expression_columns": (
            exp_schema["expression_columns"]
        ),
        "meta_columns": meta_columns,
        "sample_count": len(expression_identifiers),
        "expression_feature_count": len(
            exp_schema["expression_columns"]
        ),
    }


def validate_scst_hybrid_reference_task_params(
    *,
    task_name: str,
    data_type: str,
    tcga_type: str,
    lncrna_type: str,
    group_col: str,
    use_padj: bool,
    logfc_cutoff_mrna: float,
    padj_cutoff_mrna: float,
) -> None:
    validate_safe_name(task_name, "task_name")

    if data_type not in SCST_HYBRID_REFERENCE_VALID_DATA_TYPES:
        raise SCSTHybridReferenceTaskInputError(
            "Invalid field: data_type. Allowed values are: sc, st."
        )

    if tcga_type not in HYBRID_REFERENCE_VALID_TCGA_TYPES:
        raise SCSTHybridReferenceTaskInputError(
            "Invalid field: tcga_type. Allowed values are: "
            f"{', '.join(HYBRID_REFERENCE_VALID_TCGA_TYPES)}."
        )

    if lncrna_type not in HYBRID_REFERENCE_VALID_LNCRNA_TYPES:
        raise SCSTHybridReferenceTaskInputError(
            "Invalid field: lncrna_type. Allowed values are: "
            f"{', '.join(HYBRID_REFERENCE_VALID_LNCRNA_TYPES)}."
        )

    if not group_col:
        raise SCSTHybridReferenceTaskInputError(
            "Missing field: group_col."
        )

    validate_safe_name(group_col, "group_col")

    if not isinstance(use_padj, bool):
        raise SCSTHybridReferenceTaskInputError(
            "Invalid field: use_padj. "
            "Allowed values are true or false."
        )

    if logfc_cutoff_mrna < 0:
        raise SCSTHybridReferenceTaskInputError(
            "logfc_cutoff_mrna must be greater than or equal to 0."
        )

    if padj_cutoff_mrna <= 0 or padj_cutoff_mrna > 1:
        raise SCSTHybridReferenceTaskInputError(
            "padj_cutoff_mrna must be in the range (0, 1]."
        )


def validate_scst_uploaded_file_extensions(files) -> None:
    for field_name, allowed_suffixes in (
        SCST_HYBRID_REFERENCE_ALLOWED_FILE_SUFFIXES.items()
    ):
        uploaded_file = files.get(field_name)

        if uploaded_file is None:
            raise SCSTHybridReferenceTaskInputError(
                f"Missing uploaded file: {field_name}."
            )

        suffix = Path(uploaded_file.name).suffix.lower()

        if suffix not in allowed_suffixes:
            allowed_text = ", ".join(sorted(allowed_suffixes))

            raise SCSTHybridReferenceTaskInputError(
                f"Invalid file type for {field_name}. "
                f"Allowed file extension(s): {allowed_text}."
            )
