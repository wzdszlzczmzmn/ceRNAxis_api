from pathlib import Path
import csv

import pandas as pd

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


class HybridReferenceTaskInputError(ValueError):
    pass


class HybridReferenceTaskPathError(ValueError):
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
