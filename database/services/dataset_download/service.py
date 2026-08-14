import uuid as uuid_lib
import zipfile
from dataclasses import dataclass
from pathlib import Path

from django.conf import settings

from database.models import DatasetMetadata
from database.utils.dataset_annotation_utils.path_utils import build_timedb_group_by_candidates, \
    DatasetAnnotationInputError
from database.utils.dataset_annotation_utils.scst_metadata_utils import read_scst_dataset_group_value_counts
from database.utils.dataset_annotation_utils.scst_path_utils import validate_scst_group_value_path_component, \
    should_skip_scst_group_value_for_results, resolve_scst_dataset_group_annotation_dir, validate_scst_data_type, \
    get_scst_dataset_group_by_fields, get_scst_dataset_meta_file_path
from database.utils.expression_file_utils import (
    EXPRESSION_TYPES_BY_RNA_TYPE,
    get_expression_csv_file_path, get_timedb_expression_file_path, TIMEDB_EXPRESSION_TYPE,
    get_expression_mode_from_metadata, TISCH2_EXPRESSION_TYPE,
)
from database.utils.meta_file_utils import get_tcga_dataset_meta_file, get_timedb_dataset_meta_file, get_large_meta_file
from database.utils.viz_file_utils import (
    get_deg_file_path,
    is_tcga_dataset_metadata, get_timedb_deg_file_path, get_tisch2_deg_file_path,
)


TCGA_ANNOTATION_DEG_RNA_TYPES = [
    "mRNA",
    "miRNA",
    "lncRNA",
    "circRNA",
]


class DatasetDownloadError(Exception):
    pass


class DatasetDownloadNotSupportedError(DatasetDownloadError):
    pass


class DatasetDownloadFileNotFoundError(DatasetDownloadError):
    pass


class DatasetDownloadArchiveError(DatasetDownloadError):
    pass


@dataclass(frozen=True)
class DownloadableDatasetFile:
    path: Path
    arcname: str
    is_directory: bool = False


@dataclass(frozen=True)
class DatasetDownloadResult:
    archive_path: Path
    archive_name: str
    dataset: str


def prepare_dataset_download(dataset: str) -> DatasetDownloadResult:
    """
    Prepare a cached zip archive for dataset download.

    Supported scope:
        - TCGA datasets
        - TIMEDB datasets

    TCGA download rule:
        - {cohort}_meta.csv
        - {dataset}_{expression_type}_exp.csv
        - {dataset}_{expression_type}_deg.csv

    TIMEDB download rule:
        - {dataset}_exp.csv
        - {dataset}_meta.csv
        - {dataset}_deg.csv

    Existing files are packed into the zip archive.
    Missing files are skipped.
    Generated zip is cached under DATASET_DOWNLOAD_ZIP_DIR.
    """

    dataset = validate_dataset_param(dataset)
    metadata = get_dataset_metadata(dataset)

    archive_name = build_dataset_archive_name(dataset)
    archive_dir = get_dataset_data_download_archive_dir()
    archive_path = (archive_dir / archive_name).resolve()

    validate_output_archive_path(
        archive_path=archive_path,
        archive_dir=archive_dir,
    )

    if archive_path.exists() and archive_path.is_file():
        return DatasetDownloadResult(
            archive_path=archive_path,
            archive_name=archive_name,
            dataset=dataset,
        )

    expression_mode = get_expression_mode_from_metadata(metadata)

    if is_tcga_dataset_metadata(metadata):
        downloadable_files = resolve_tcga_dataset_download_files(metadata)
    elif expression_mode == "timedb":
        downloadable_files = resolve_timedb_dataset_download_files(metadata)
    elif expression_mode == "tisch2":
        downloadable_files = resolve_tisch2_dataset_download_files(metadata)
    elif expression_mode == "scTML":
        downloadable_files = resolve_sctml_dataset_download_files(metadata)
    else:
        raise DatasetDownloadNotSupportedError(
            "Dataset download is currently only supported for TCGA, TIMEDB, TISCH2, and scTML datasets. "
            f"Dataset '{dataset}' has expression_mode '{expression_mode}'."
        )

    if not downloadable_files:
        raise DatasetDownloadFileNotFoundError(
            f"No downloadable files found for dataset '{dataset}'."
        )

    create_cached_dataset_archive(
        archive_path=archive_path,
        downloadable_files=downloadable_files,
    )

    return DatasetDownloadResult(
        archive_path=archive_path,
        archive_name=archive_name,
        dataset=dataset,
    )


def prepare_tcga_annotation_download(dataset: str) -> DatasetDownloadResult:
    """
    Prepare a cached zip archive for TCGA dataset annotation download.

    Input dataset must be an mRNA TCGA dataset, for example:
        TCGA_ACC_mRNA

    Annotation files are generated from Module 2-like TCGA annotation results.

    Download rule for TCGA_ACC_mRNA:
        TCGA_ACC_ceRNA_axis.csv
        TCGA_ACC_ceRNA_axis_final.csv
        TCGA_ACC_ceRNA_background.csv
        TCGA_ACC_ceRNA_corr.csv
        TCGA_ACC_ceRNA_network.csv
        TCGA_ACC_CMap.csv
        TCGA_ACC_limma_mRNA.csv
        TCGA_ACC_limma_miRNA.csv
        TCGA_ACC_limma_lncRNA.csv
        TCGA_ACC_limma_circRNA.csv
        TCGA_ACC_map_immune_axis.csv
        TCGA_ACC_mRNA_gsea.csv
        TCGA_ACC_survival_analysis.csv
        TCGA_ACC_sponge_result.csv
        TCGA_ACC_CMdrug_result/

    Existing files/directories are packed.
    Missing files/directories are skipped.
    Generated zip is cached under DATASET_DOWNLOAD_ZIP_DIR.
    """

    dataset = validate_dataset_param(dataset)
    metadata = get_dataset_metadata(dataset)

    validate_tcga_annotation_dataset(metadata)

    annotation_prefix = get_tcga_annotation_prefix(dataset)

    archive_name = build_tcga_annotation_archive_name(annotation_prefix)
    archive_dir = get_dataset_annotation_download_archive_dir()
    archive_path = (archive_dir / archive_name).resolve()

    validate_output_archive_path(
        archive_path=archive_path,
        archive_dir=archive_dir,
    )

    if archive_path.exists() and archive_path.is_file():
        return DatasetDownloadResult(
            archive_path=archive_path,
            archive_name=archive_name,
            dataset=dataset,
        )

    downloadable_files = resolve_tcga_annotation_download_files(
        dataset=dataset,
        annotation_prefix=annotation_prefix,
    )

    if not downloadable_files:
        create_empty_annotation_archive(
            archive_path=archive_path,
            archive_root=build_tcga_annotation_archive_root(
                annotation_prefix
            ),
            dataset=dataset,
        )
    else:
        create_cached_dataset_archive(
            archive_path=archive_path,
            downloadable_files=downloadable_files,
        )

    return DatasetDownloadResult(
        archive_path=archive_path,
        archive_name=archive_name,
        dataset=dataset,
    )


def prepare_timedb_annotation_download(dataset: str) -> DatasetDownloadResult:
    """
    Prepare a cached zip archive for TIMEDB dataset annotation download.

    TIMEDB annotation results are generated by Module 3-like Hybrid Reference logic.

    Download rule for GSE19750:
        GSE19750_ceRNA_axis.csv
        GSE19750_ceRNA_axis_final.csv
        GSE19750_ceRNA_background.csv
        GSE19750_ceRNA_corr.csv
        GSE19750_ceRNA_network.csv
        GSE19750_CMap.csv
        GSE19750_map_immune_axis.csv
        GSE19750_mRNA_gsea.csv
        GSE19750_survival_analysis.csv
        GSE19750_limma_mRNA.csv
        GSE19750_limma_mRNA_intersect.csv
        GSE19750_limma_mRNA_venn.csv
        GSE19750_CMdrug_result/

    Existing files/directories are packed.
    Missing files/directories are skipped.
    Generated zip is cached under DATASET_DOWNLOAD_ZIP_DIR.
    """

    dataset = validate_dataset_param(dataset)
    metadata = get_dataset_metadata(dataset)

    validate_timedb_annotation_dataset(metadata)

    annotation_prefix = dataset

    archive_name = build_timedb_annotation_archive_name(annotation_prefix)
    archive_dir = get_dataset_annotation_download_archive_dir()
    archive_path = (archive_dir / archive_name).resolve()

    validate_output_archive_path(
        archive_path=archive_path,
        archive_dir=archive_dir,
    )

    if archive_path.exists() and archive_path.is_file():
        return DatasetDownloadResult(
            archive_path=archive_path,
            archive_name=archive_name,
            dataset=dataset,
        )

    downloadable_files = resolve_timedb_annotation_download_files(
        annotation_prefix=annotation_prefix,
    )

    if not downloadable_files:
        create_empty_annotation_archive(
            archive_path=archive_path,
            archive_root=build_timedb_annotation_archive_root(
                annotation_prefix
            ),
            dataset=dataset,
        )
    else:
        create_cached_dataset_archive(
            archive_path=archive_path,
            downloadable_files=downloadable_files,
        )

    return DatasetDownloadResult(
        archive_path=archive_path,
        archive_name=archive_name,
        dataset=dataset,
    )


def prepare_scst_annotation_download(
    dataset: str,
) -> DatasetDownloadResult:
    """
    Prepare a cached zip archive for SC/ST Dataset Annotation download.

    Supported sources:
        - TISCH2 -> data_type="sc"
        - scTML  -> data_type="st"

    Annotation results are grouped by configured group_by fields,
    with each group-by stored in its own annotation directory.

    Existing files/directories are packed.
    Missing files/directories are skipped.
    Generated zip is cached under DATASET_DOWNLOAD_ZIP_DIR.
    """

    dataset = validate_dataset_param(dataset)
    metadata = get_dataset_metadata(dataset)

    expression_mode = get_expression_mode_from_metadata(metadata)

    if expression_mode == "tisch2":
        data_type = "sc"
    elif expression_mode == "scTML":
        data_type = "st"
    else:
        raise DatasetDownloadNotSupportedError(
            "SC/ST annotation download is only supported for "
            "TISCH2 and scTML datasets. "
            f"Dataset '{dataset}' has expression_mode "
            f"'{expression_mode}'."
        )

    archive_name = build_scst_annotation_archive_name(
        dataset
    )

    archive_dir = get_dataset_annotation_download_archive_dir()

    archive_path = (
        archive_dir / archive_name
    ).resolve()

    validate_output_archive_path(
        archive_path=archive_path,
        archive_dir=archive_dir,
    )

    if archive_path.exists() and archive_path.is_file():
        return DatasetDownloadResult(
            archive_path=archive_path,
            archive_name=archive_name,
            dataset=dataset,
        )

    downloadable_files = (
        resolve_scst_annotation_download_files(
            dataset_name=dataset,
            data_type=data_type,
        )
    )

    if not downloadable_files:
        create_empty_annotation_archive(
            archive_path=archive_path,
            archive_root=build_scst_annotation_archive_root(
                dataset
            ),
            dataset=dataset,
        )
    else:
        create_cached_dataset_archive(
            archive_path=archive_path,
            downloadable_files=downloadable_files,
        )

    return DatasetDownloadResult(
        archive_path=archive_path,
        archive_name=archive_name,
        dataset=dataset,
    )


def create_empty_annotation_archive(
    *,
    archive_path: Path,
    archive_root: str,
    dataset: str,
) -> None:
    archive_dir = archive_path.parent.resolve()

    temp_archive_path = (
        archive_dir
        / f".{archive_path.name}.{uuid_lib.uuid4().hex}.tmp"
    ).resolve()

    validate_output_archive_path(
        archive_path=temp_archive_path,
        archive_dir=archive_dir,
    )

    try:
        with zipfile.ZipFile(
            temp_archive_path,
            mode="w",
            compression=zipfile.ZIP_DEFLATED,
        ) as zip_file:
            filename = "NO_AVAILABLE_ANNOTATION_FILES.txt"

            arcname = (
                f"{archive_root}/{filename}"
            )

            validate_archive_name(arcname)

            content = (
                "No downloadable annotation files are currently "
                f"available for dataset: {dataset}.\n"
            )

            zip_file.writestr(
                arcname,
                content,
            )

        temp_archive_path.replace(
            archive_path
        )

    except Exception as e:
        if temp_archive_path.exists():
            temp_archive_path.unlink()

        raise DatasetDownloadArchiveError(
            f"Failed to create annotation archive: {str(e)}"
        ) from e


def validate_dataset_param(dataset: str) -> str:
    dataset = str(dataset or "").strip()

    if not dataset:
        raise DatasetDownloadError("Missing required parameter: dataset.")

    if "/" in dataset or "\\" in dataset or ".." in dataset:
        raise DatasetDownloadError("Invalid dataset parameter.")

    return dataset


def get_dataset_metadata(dataset: str) -> DatasetMetadata:
    try:
        return DatasetMetadata.objects.get(dataset=dataset)
    except DatasetMetadata.DoesNotExist as e:
        raise DatasetDownloadError(
            f"Dataset metadata not found for dataset '{dataset}'."
        ) from e


def get_dataset_download_zip_base_dir() -> Path:
    """
    Base output directory for reusable dataset-related zip files.

    settings.py:
        DATASET_DOWNLOAD_ZIP_DIR = Path(BASE_DIR) / "zip_files" / "dataset"
    """

    archive_dir = getattr(
        settings,
        "DATASET_DOWNLOAD_ZIP_DIR",
        Path(settings.BASE_DIR) / "zip_files" / "dataset",
    )

    archive_dir = Path(archive_dir).resolve()
    archive_dir.mkdir(parents=True, exist_ok=True)

    return archive_dir


def get_dataset_data_download_archive_dir() -> Path:
    """
    Cache directory for dataset data zip files.
    """

    archive_dir = (get_dataset_download_zip_base_dir() / "data").resolve()
    archive_dir.mkdir(parents=True, exist_ok=True)

    return archive_dir


def get_dataset_annotation_download_archive_dir() -> Path:
    """
    Cache directory for dataset annotation zip files.
    """

    archive_dir = (get_dataset_download_zip_base_dir() / "annotation").resolve()
    archive_dir.mkdir(parents=True, exist_ok=True)

    return archive_dir


def build_dataset_archive_name(dataset: str) -> str:
    safe_dataset = sanitize_archive_filename_part(dataset)
    return f"{safe_dataset}_dataset.zip"


def build_dataset_archive_root(dataset: str) -> str:
    safe_dataset = sanitize_archive_filename_part(dataset)
    return f"{safe_dataset}_dataset"


def resolve_tcga_dataset_download_files(
    metadata: DatasetMetadata,
) -> list[DownloadableDatasetFile]:
    dataset = metadata.dataset
    rna_type = metadata.gene_bio_type
    archive_root = build_dataset_archive_root(dataset)

    expression_types = EXPRESSION_TYPES_BY_RNA_TYPE.get(rna_type)

    if not expression_types:
        raise DatasetDownloadError(
            f"Invalid RNA type '{rna_type}' for dataset '{dataset}'."
        )

    downloadable_files: list[DownloadableDatasetFile] = []

    meta_file_path = get_tcga_dataset_meta_file(dataset)

    if meta_file_path.exists() and meta_file_path.is_file():
        downloadable_files.append(
            DownloadableDatasetFile(
                path=meta_file_path,
                arcname=f"{archive_root}/{meta_file_path.name}",
            )
        )

    for expression_type in sorted(expression_types):
        expression_file_path = get_expression_csv_file_path(
            dataset=dataset,
            rna_type=rna_type,
            expression_type=expression_type,
        )

        if expression_file_path.exists() and expression_file_path.is_file():
            downloadable_files.append(
                DownloadableDatasetFile(
                    path=expression_file_path,
                    arcname=f"{archive_root}/{expression_file_path.name}",
                )
            )

        deg_file_path = get_deg_file_path(
            dataset=dataset,
            rna_type=rna_type,
            expression_type=expression_type,
        )

        if deg_file_path.exists() and deg_file_path.is_file():
            downloadable_files.append(
                DownloadableDatasetFile(
                    path=deg_file_path,
                    arcname=f"{archive_root}/{deg_file_path.name}",
                )
            )

    return validate_dataset_downloadable_files(downloadable_files)


def resolve_timedb_dataset_download_files(
    metadata: DatasetMetadata,
) -> list[DownloadableDatasetFile]:
    """
    TIMEDB dataset download files.

    Download rule:
        {dataset}_exp.csv
        {dataset}_meta.csv
        {dataset}_deg.csv

    Existing files are packed.
    Missing files are skipped.
    """

    dataset = metadata.dataset
    rna_type = metadata.gene_bio_type
    archive_root = build_dataset_archive_root(dataset)

    downloadable_files: list[DownloadableDatasetFile] = []

    expression_file_path = get_timedb_expression_file_path(
        dataset=dataset,
        rna_type=rna_type,
        file_format="csv",
    )

    if expression_file_path.exists() and expression_file_path.is_file():
        downloadable_files.append(
            DownloadableDatasetFile(
                path=expression_file_path,
                arcname=f"{archive_root}/{expression_file_path.name}",
            )
        )

    meta_file_path = get_timedb_dataset_meta_file(dataset)

    if meta_file_path.exists() and meta_file_path.is_file():
        downloadable_files.append(
            DownloadableDatasetFile(
                path=meta_file_path,
                arcname=f"{archive_root}/{meta_file_path.name}",
            )
        )

    deg_file_path = get_timedb_deg_file_path(
        dataset=dataset,
        rna_type=rna_type,
        expression_type=TIMEDB_EXPRESSION_TYPE,
    )

    if deg_file_path.exists() and deg_file_path.is_file():
        downloadable_files.append(
            DownloadableDatasetFile(
                path=deg_file_path,
                arcname=f"{archive_root}/{deg_file_path.name}",
            )
        )

    return validate_dataset_downloadable_files(downloadable_files)


def resolve_tisch2_dataset_download_files(
    metadata: DatasetMetadata,
) -> list[DownloadableDatasetFile]:
    """
    TISCH2 dataset download files.

    Download rule:
        {dataset}_sparse.h5ad
        {dataset}_meta.csv
        {dataset}_deg.csv

    Existing files are packed.
    Missing files are skipped.
    """

    dataset = metadata.dataset
    rna_type = metadata.gene_bio_type
    expression_mode = get_expression_mode_from_metadata(metadata)
    archive_root = build_dataset_archive_root(dataset)

    downloadable_files: list[DownloadableDatasetFile] = []

    meta_file_path = get_large_meta_file(
        dataset=dataset,
        expression_mode=expression_mode,
    )

    sparse_h5ad_path = (
        meta_file_path.parent / f"{dataset}_sparse.h5ad"
    ).resolve()

    if sparse_h5ad_path.exists() and sparse_h5ad_path.is_file():
        downloadable_files.append(
            DownloadableDatasetFile(
                path=sparse_h5ad_path,
                arcname=f"{archive_root}/{sparse_h5ad_path.name}",
            )
        )

    if meta_file_path.exists() and meta_file_path.is_file():
        downloadable_files.append(
            DownloadableDatasetFile(
                path=meta_file_path,
                arcname=f"{archive_root}/{meta_file_path.name}",
            )
        )

    deg_file_path = get_tisch2_deg_file_path(
        dataset=dataset,
        rna_type=rna_type,
        expression_type=TISCH2_EXPRESSION_TYPE,
    )

    if deg_file_path.exists() and deg_file_path.is_file():
        downloadable_files.append(
            DownloadableDatasetFile(
                path=deg_file_path,
                arcname=f"{archive_root}/{deg_file_path.name}",
            )
        )

    return validate_dataset_downloadable_files(downloadable_files)


def resolve_sctml_dataset_download_files(
    metadata: DatasetMetadata,
) -> list[DownloadableDatasetFile]:
    """
    scTML dataset download files.

    Download rule:
        {dataset}_sparse.h5ad
        {dataset}_meta.csv

    Existing files are packed.
    Missing files are skipped.
    """

    dataset = metadata.dataset
    expression_mode = get_expression_mode_from_metadata(metadata)
    archive_root = build_dataset_archive_root(dataset)

    downloadable_files: list[DownloadableDatasetFile] = []

    meta_file_path = get_large_meta_file(
        dataset=dataset,
        expression_mode=expression_mode,
    )

    sparse_h5ad_path = (
        meta_file_path.parent / f"{dataset}_sparse.h5ad"
    ).resolve()

    if sparse_h5ad_path.exists() and sparse_h5ad_path.is_file():
        downloadable_files.append(
            DownloadableDatasetFile(
                path=sparse_h5ad_path,
                arcname=f"{archive_root}/{sparse_h5ad_path.name}",
            )
        )

    if meta_file_path.exists() and meta_file_path.is_file():
        downloadable_files.append(
            DownloadableDatasetFile(
                path=meta_file_path,
                arcname=f"{archive_root}/{meta_file_path.name}",
            )
        )

    return validate_dataset_downloadable_files(downloadable_files)


def validate_tcga_annotation_dataset(metadata: DatasetMetadata) -> None:
    dataset = metadata.dataset

    if not is_tcga_dataset_metadata(metadata):
        raise DatasetDownloadNotSupportedError(
            "TCGA annotation download is only supported for TCGA datasets. "
            f"Dataset '{dataset}' is not a TCGA dataset."
        )

    if metadata.gene_bio_type != "mRNA":
        raise DatasetDownloadNotSupportedError(
            "TCGA annotation download requires an mRNA dataset. "
            f"Dataset '{dataset}' has gene_bio_type '{metadata.gene_bio_type}'."
        )


def get_tcga_annotation_prefix(dataset: str) -> str:
    """
    Example:
        TCGA_ACC_mRNA -> TCGA_ACC
    """

    parts = str(dataset).split("_")

    if len(parts) < 3:
        raise DatasetDownloadError(
            f"Invalid TCGA mRNA dataset name for annotation download: '{dataset}'."
        )

    return "_".join(parts[:2])


def get_tcga_annotation_dataset_result_dir(annotation_prefix: str) -> Path:
    """
    Source directory for one TCGA annotation dataset.

    Example:
        TCGA_DATASET_ANNOTATIONS_DIR / "TCGA_ACC"
    """

    annotation_prefix = validate_dataset_param(annotation_prefix)

    base_dir = get_tcga_annotation_result_dir()
    result_dir = (base_dir / annotation_prefix).resolve()

    if not is_path_under_dir(
        file_path=result_dir,
        base_dir=base_dir,
    ):
        raise DatasetDownloadArchiveError(
            "Invalid TCGA annotation dataset result directory."
        )

    return result_dir


def get_tcga_annotation_result_dir() -> Path:
    """
    Source directory for TCGA dataset annotation result files.

    Uses:
        settings.TCGA_DATASET_ANNOTATIONS_DIR
    """

    base_dir = Path(settings.TCGA_DATASET_ANNOTATIONS_DIR).resolve()

    if not base_dir.exists() or not base_dir.is_dir():
        raise DatasetDownloadFileNotFoundError(
            "TCGA dataset annotation source directory not found."
        )

    return base_dir


def build_tcga_annotation_archive_name(annotation_prefix: str) -> str:
    safe_prefix = sanitize_archive_filename_part(annotation_prefix)
    return f"{safe_prefix}_annotation.zip"


def build_tcga_annotation_archive_root(annotation_prefix: str) -> str:
    safe_prefix = sanitize_archive_filename_part(annotation_prefix)
    return f"{safe_prefix}_annotation"


def resolve_tcga_annotation_download_files(
    dataset: str,
    annotation_prefix: str,
) -> list[DownloadableDatasetFile]:
    """
    TCGA annotation result files.

    Example:
        dataset: TCGA_ACC_mRNA
        annotation_prefix: TCGA_ACC

    Source directory:
        TCGA_DATASET_ANNOTATIONS_DIR / TCGA_ACC /

    Missing files/directories are skipped.
    """

    validate_dataset_param(dataset)
    annotation_prefix = validate_dataset_param(annotation_prefix)

    annotation_base_dir = get_tcga_annotation_dataset_result_dir(
        annotation_prefix
    )
    archive_root = build_tcga_annotation_archive_root(annotation_prefix)

    result_filenames = [
        f"{annotation_prefix}_ceRNA_axis.csv",
        f"{annotation_prefix}_ceRNA_axis_final.csv",
        f"{annotation_prefix}_ceRNA_background.csv",
        f"{annotation_prefix}_ceRNA_corr.csv",
        f"{annotation_prefix}_ceRNA_network.csv",
        f"{annotation_prefix}_CMap.csv",
        f"{annotation_prefix}_map_immune_axis.csv",
        f"{annotation_prefix}_mRNA_gsea.csv",
        f"{annotation_prefix}_survival_analysis.csv",
        f"{annotation_prefix}_sponge_result.csv",
    ]

    downloadable_files: list[DownloadableDatasetFile] = []

    for filename in result_filenames:
        file_path = (annotation_base_dir / filename).resolve()

        if file_path.exists() and file_path.is_file():
            downloadable_files.append(
                DownloadableDatasetFile(
                    path=file_path,
                    arcname=f"{archive_root}/{filename}",
                )
            )

    for rna_type in TCGA_ANNOTATION_DEG_RNA_TYPES:
        filename = f"{annotation_prefix}_limma_{rna_type}.csv"
        file_path = (annotation_base_dir / filename).resolve()

        if file_path.exists() and file_path.is_file():
            downloadable_files.append(
                DownloadableDatasetFile(
                    path=file_path,
                    arcname=f"{archive_root}/{filename}",
                )
            )

    cm_drug_result_dirname = f"{annotation_prefix}_CMdrug_result"
    cm_drug_result_dir = (
        annotation_base_dir / cm_drug_result_dirname
    ).resolve()

    if cm_drug_result_dir.exists() and cm_drug_result_dir.is_dir():
        downloadable_files.append(
            DownloadableDatasetFile(
                path=cm_drug_result_dir,
                arcname=f"{archive_root}/{cm_drug_result_dirname}/",
                is_directory=True,
            )
        )

    return validate_dataset_downloadable_files(
        downloadable_files=downloadable_files,
        allowed_base_dir=annotation_base_dir,
    )


def validate_timedb_annotation_dataset(metadata: DatasetMetadata) -> None:
    dataset = metadata.dataset
    expression_mode = get_expression_mode_from_metadata(metadata)

    if expression_mode != "timedb":
        raise DatasetDownloadNotSupportedError(
            "TIMEDB annotation download is only supported for TIMEDB datasets. "
            f"Dataset '{dataset}' has expression_mode '{expression_mode}'."
        )

    if metadata.gene_bio_type != "mRNA":
        raise DatasetDownloadNotSupportedError(
            "TIMEDB annotation download requires an mRNA dataset. "
            f"Dataset '{dataset}' has gene_bio_type '{metadata.gene_bio_type}'."
        )


def get_timedb_annotation_dataset_result_dir(annotation_prefix: str) -> Path:
    """
    Source directory for one TIMEDB annotation dataset.

    Example:
        TIMEDB_DATASET_ANNOTATIONS_DIR / "GSE19750"
    """

    annotation_prefix = validate_dataset_param(annotation_prefix)

    base_dir = get_timedb_annotation_result_dir()
    result_dir = (base_dir / annotation_prefix).resolve()

    if not is_path_under_dir(
        file_path=result_dir,
        base_dir=base_dir,
    ):
        raise DatasetDownloadArchiveError(
            "Invalid TIMEDB annotation dataset result directory."
        )

    return result_dir


def get_timedb_annotation_result_dir() -> Path:
    """
    Source directory for TIMEDB dataset annotation result files.

    Uses:
        settings.TIMEDB_DATASET_ANNOTATIONS_DIR
    """

    base_dir = Path(settings.TIMEDB_DATASET_ANNOTATIONS_DIR).resolve()

    if not base_dir.exists() or not base_dir.is_dir():
        raise DatasetDownloadFileNotFoundError(
            "TIMEDB dataset annotation source directory not found."
        )

    return base_dir


def build_timedb_annotation_archive_name(annotation_prefix: str) -> str:
    safe_prefix = sanitize_archive_filename_part(annotation_prefix)
    return f"{safe_prefix}_annotation.zip"


def build_timedb_annotation_archive_root(annotation_prefix: str) -> str:
    safe_prefix = sanitize_archive_filename_part(annotation_prefix)
    return f"{safe_prefix}_annotation"


def resolve_timedb_group_annotation_download_files(
    *,
    annotation_base_dir: Path,
    annotation_prefix: str,
    archive_group_root: str,
) -> list[DownloadableDatasetFile]:
    """
    Resolve downloadable annotation files for one TIMEDB group directory.
    """

    result_filenames = [
        f"{annotation_prefix}_ceRNA_axis.csv",
        f"{annotation_prefix}_ceRNA_axis_final.csv",
        f"{annotation_prefix}_ceRNA_background.csv",
        f"{annotation_prefix}_ceRNA_corr.csv",
        f"{annotation_prefix}_ceRNA_network.csv",
        f"{annotation_prefix}_CMap.csv",
        f"{annotation_prefix}_map_immune_axis.csv",
        f"{annotation_prefix}_mRNA_gsea.csv",
        f"{annotation_prefix}_survival_analysis.csv",
        f"{annotation_prefix}_sponge_result.csv",
    ]

    deg_method = "limma"

    deg_filenames = [
        f"{annotation_prefix}_{deg_method}_mRNA.csv",
        f"{annotation_prefix}_{deg_method}_mRNA_intersect.csv",
        f"{annotation_prefix}_{deg_method}_mRNA_venn.csv",
    ]

    downloadable_files: list[DownloadableDatasetFile] = []

    for filename in result_filenames + deg_filenames:
        file_path = (annotation_base_dir / filename).resolve()

        if file_path.exists() and file_path.is_file():
            downloadable_files.append(
                DownloadableDatasetFile(
                    path=file_path,
                    arcname=f"{archive_group_root}/{filename}",
                )
            )

    cm_drug_result_dirname = f"{annotation_prefix}_CMdrug_result"

    cm_drug_result_dir = (
        annotation_base_dir / cm_drug_result_dirname
    ).resolve()

    if cm_drug_result_dir.exists() and cm_drug_result_dir.is_dir():
        downloadable_files.append(
            DownloadableDatasetFile(
                path=cm_drug_result_dir,
                arcname=(
                    f"{archive_group_root}/"
                    f"{cm_drug_result_dirname}/"
                ),
                is_directory=True,
            )
        )

    return validate_dataset_downloadable_files(
        downloadable_files=downloadable_files,
        allowed_base_dir=annotation_base_dir,
    )



def resolve_timedb_annotation_download_files(
    annotation_prefix: str,
) -> list[DownloadableDatasetFile]:
    """
    TIMEDB annotation result files.

    Group directory mapping:
        other -> {dataset_name}
        grade -> {dataset_name}_grade
        stage -> {dataset_name}_stage

    Each available group directory uses the same result-file rules.
    """

    annotation_prefix = validate_dataset_param(annotation_prefix)

    annotation_root_dir = get_timedb_annotation_result_dir()
    archive_root = build_timedb_annotation_archive_root(annotation_prefix)

    downloadable_files: list[DownloadableDatasetFile] = []

    for candidate in build_timedb_group_by_candidates(annotation_prefix):
        annotation_dir_name = candidate["annotation_dir_name"]
        file_prefix = candidate["file_prefix"]

        annotation_base_dir = (
                annotation_root_dir / annotation_dir_name
        ).resolve()

        if not is_path_under_dir(
                file_path=annotation_base_dir,
                base_dir=annotation_root_dir,
        ):
            raise DatasetDownloadArchiveError(
                "Invalid TIMEDB annotation group directory."
            )

        if (
                not annotation_base_dir.exists()
                or not annotation_base_dir.is_dir()
        ):
            continue

        archive_group_root = (
            f"{archive_root}/{annotation_dir_name}"
        )

        downloadable_files.extend(
            resolve_timedb_group_annotation_download_files(
                annotation_base_dir=annotation_base_dir,
                annotation_prefix=file_prefix,
                archive_group_root=archive_group_root,
            )
        )

    return downloadable_files


def build_scst_annotation_archive_name(
    dataset_name: str,
) -> str:
    safe_dataset = sanitize_archive_filename_part(
        dataset_name
    )

    return f"{safe_dataset}_annotation.zip"


def build_scst_annotation_archive_root(
    dataset_name: str,
) -> str:
    safe_dataset = sanitize_archive_filename_part(
        dataset_name
    )

    return f"{safe_dataset}_annotation"


def resolve_scst_group_annotation_download_files(
    *,
    dataset_name: str,
    annotation_dir: Path,
    group_values,
    archive_group_root: str,
) -> list[DownloadableDatasetFile]:
    """
    Resolve downloadable files for one SC/ST Dataset Annotation
    group-by directory.

    Each valid group_value follows the same output-file naming
    convention as SCSTHybridReferenceTask.
    """

    downloadable_files: list[DownloadableDatasetFile] = []

    for raw_group_value in group_values:
        group_value = str(raw_group_value or "").strip()

        if not group_value:
            continue

        # "/" values are explicitly unsupported by the upstream
        # SC/ST result-generation workflow.
        if should_skip_scst_group_value_for_results(group_value):
            continue

        # Dataset metadata is external data. A bad group value should
        # not make the whole annotation download fail.
        try:
            group_value = validate_scst_group_value_path_component(
                group_value
            )
        except DatasetAnnotationInputError:
            continue

        result_filenames = [
            f"{dataset_name}_ceRNA_corr_{group_value}.csv",
            f"{dataset_name}_ceRNA_{group_value}_axis.csv",
            f"{dataset_name}_ceRNA_{group_value}_axis_final.csv",
            f"{dataset_name}_ceRNA_{group_value}_background.csv",
            f"{dataset_name}_ceRNA_{group_value}_network.csv",
            f"{dataset_name}_CMap_{group_value}.csv",
            f"{dataset_name}_deg_{group_value}.csv",
            f"{dataset_name}_map_immune_axis_{group_value}.csv",
            f"{dataset_name}_mRNA_deg_{group_value}_intersect.csv",
            f"{dataset_name}_mRNA_deg_{group_value}_venn.csv",
            f"{dataset_name}_mRNA_gsea_{group_value}.csv",
            f"{dataset_name}_survival_analysis_{group_value}.csv",
        ]

        for filename in result_filenames:
            file_path = (
                annotation_dir / filename
            ).resolve()

            if file_path.exists() and file_path.is_file():
                downloadable_files.append(
                    DownloadableDatasetFile(
                        path=file_path,
                        arcname=(
                            f"{archive_group_root}/{filename}"
                        ),
                    )
                )

        cm_drug_result_dirname = (
            f"{dataset_name}_CMdrug_result_{group_value}"
        )

        cm_drug_result_dir = (
            annotation_dir / cm_drug_result_dirname
        ).resolve()

        if (
            cm_drug_result_dir.exists()
            and cm_drug_result_dir.is_dir()
        ):
            downloadable_files.append(
                DownloadableDatasetFile(
                    path=cm_drug_result_dir,
                    arcname=(
                        f"{archive_group_root}/"
                        f"{cm_drug_result_dirname}/"
                    ),
                    is_directory=True,
                )
            )

    return validate_dataset_downloadable_files(
        downloadable_files=downloadable_files,
        allowed_base_dir=annotation_dir,
    )


def resolve_scst_annotation_download_files(
    *,
    dataset_name: str,
    data_type: str,
) -> list[DownloadableDatasetFile]:
    """
    Resolve SC/ST Dataset Annotation download files.

    Hierarchy:
        dataset
            -> group_by annotation directory
                -> group_value result files

    Group-by configuration and group values are resolved from the
    existing SC/ST Dataset Annotation metadata utilities.
    """

    dataset_name = validate_dataset_param(dataset_name)
    data_type = validate_scst_data_type(data_type)

    archive_root = build_scst_annotation_archive_root(
        dataset_name
    )

    group_by_fields = get_scst_dataset_group_by_fields(
        dataset_name
    )

    if not group_by_fields:
        return []

    meta_file = get_scst_dataset_meta_file_path(
        dataset_name=dataset_name,
        data_type=data_type,
    )

    (
        group_counts_by_group_by,
        _sample_count,
    ) = read_scst_dataset_group_value_counts(
        meta_file=meta_file,
        group_by_fields=group_by_fields,
    )

    downloadable_files: list[DownloadableDatasetFile] = []

    for group_by in group_by_fields:
        annotation_dir = (
            resolve_scst_dataset_group_annotation_dir(
                dataset_name=dataset_name,
                group_by=group_by,
                data_type=data_type,
            )
        )

        if (
            not annotation_dir.exists()
            or not annotation_dir.is_dir()
        ):
            continue

        group_counts = group_counts_by_group_by.get(
            group_by
        )

        if not group_counts:
            continue

        # Keep the physical group-by directory structure in the zip.
        archive_group_root = (
            f"{archive_root}/{annotation_dir.name}"
        )

        group_downloadable_files = (
            resolve_scst_group_annotation_download_files(
                dataset_name=dataset_name,
                annotation_dir=annotation_dir,
                group_values=group_counts.keys(),
                archive_group_root=archive_group_root,
            )
        )

        downloadable_files.extend(
            group_downloadable_files
        )

    return downloadable_files


def validate_dataset_downloadable_files(
    downloadable_files: list[DownloadableDatasetFile],
    allowed_base_dir: Path | None = None,
) -> list[DownloadableDatasetFile]:
    validated_files = []

    if allowed_base_dir is not None:
        allowed_base_dir = Path(allowed_base_dir).resolve()

    for downloadable_file in downloadable_files:
        file_path = Path(downloadable_file.path).resolve()

        if allowed_base_dir is not None and not is_path_under_dir(
            file_path=file_path,
            base_dir=allowed_base_dir,
        ):
            raise DatasetDownloadArchiveError(
                f"Invalid downloadable path: {file_path.name}"
            )

        if file_path.is_symlink():
            raise DatasetDownloadArchiveError(
                f"Symbolic link is not allowed: {file_path.name}"
            )

        validate_archive_name(downloadable_file.arcname)

        if downloadable_file.is_directory:
            if not file_path.exists():
                continue

            if not file_path.is_dir():
                raise DatasetDownloadArchiveError(
                    f"Expected directory but found file: {file_path.name}"
                )

            validated_files.append(
                DownloadableDatasetFile(
                    path=file_path,
                    arcname=downloadable_file.arcname,
                    is_directory=True,
                )
            )
            continue

        if not file_path.exists() or not file_path.is_file():
            continue

        validated_files.append(
            DownloadableDatasetFile(
                path=file_path,
                arcname=downloadable_file.arcname,
            )
        )

    return validated_files


def write_directory_to_zip(
    zip_file: zipfile.ZipFile,
    directory_path: Path,
    directory_arcname: str,
) -> None:
    directory_path = Path(directory_path).resolve()
    directory_arcname = str(directory_arcname).rstrip("/") + "/"

    validate_archive_name(directory_arcname)

    if not directory_path.exists():
        return

    if not directory_path.is_dir():
        raise DatasetDownloadArchiveError(
            f"Expected directory but found file: {directory_path.name}"
        )

    if directory_path.is_symlink():
        raise DatasetDownloadArchiveError(
            f"Symbolic link is not allowed: {directory_path.name}"
        )

    zip_file.writestr(directory_arcname, "")

    for child_path in directory_path.rglob("*"):
        if child_path.is_symlink():
            raise DatasetDownloadArchiveError(
                f"Symbolic link is not allowed: {child_path.name}"
            )

        child_path = child_path.resolve()

        if not is_path_under_dir(
            file_path=child_path,
            base_dir=directory_path,
        ):
            raise DatasetDownloadArchiveError(
                f"Invalid directory child path: {child_path.name}"
            )

        child_relative_path = child_path.relative_to(directory_path)

        if child_path.is_dir():
            child_arcname = (
                f"{directory_arcname}"
                f"{child_relative_path.as_posix().rstrip('/')}/"
            )
            validate_archive_name(child_arcname)
            zip_file.writestr(child_arcname, "")
            continue

        if child_path.is_file():
            child_arcname = (
                f"{directory_arcname}"
                f"{child_relative_path.as_posix()}"
            )
            validate_archive_name(child_arcname)
            zip_file.write(
                filename=child_path,
                arcname=child_arcname,
            )


def create_cached_dataset_archive(
    archive_path: Path,
    downloadable_files: list[DownloadableDatasetFile],
) -> None:
    archive_dir = archive_path.parent.resolve()

    temp_archive_path = (
        archive_dir / f".{archive_path.name}.{uuid_lib.uuid4().hex}.tmp"
    ).resolve()

    validate_output_archive_path(
        archive_path=temp_archive_path,
        archive_dir=archive_dir,
    )

    try:
        with zipfile.ZipFile(
            temp_archive_path,
            mode="w",
            compression=zipfile.ZIP_DEFLATED,
        ) as zip_file:
            for downloadable_file in downloadable_files:
                file_path = Path(downloadable_file.path).resolve()

                if file_path == archive_path:
                    continue

                if file_path == temp_archive_path:
                    continue

                if file_path.is_symlink():
                    raise DatasetDownloadArchiveError(
                        f"Symbolic link is not allowed: {file_path.name}"
                    )

                if downloadable_file.is_directory:
                    write_directory_to_zip(
                        zip_file=zip_file,
                        directory_path=file_path,
                        directory_arcname=downloadable_file.arcname,
                    )
                    continue

                zip_file.write(
                    filename=file_path,
                    arcname=downloadable_file.arcname,
                )

        temp_archive_path.replace(archive_path)

    except Exception as e:
        if temp_archive_path.exists():
            temp_archive_path.unlink()

        raise DatasetDownloadArchiveError(
            f"Failed to create dataset archive: {str(e)}"
        ) from e


def validate_output_archive_path(
    archive_path: Path,
    archive_dir: Path,
) -> None:
    archive_path = Path(archive_path).resolve()
    archive_dir = Path(archive_dir).resolve()

    if not is_path_under_dir(archive_path, archive_dir):
        raise DatasetDownloadArchiveError(
            "Invalid archive output path."
        )


def is_path_under_dir(file_path: Path, base_dir: Path) -> bool:
    try:
        file_path.resolve().relative_to(base_dir.resolve())
        return True
    except ValueError:
        return False


def validate_archive_name(arcname: str) -> None:
    arcname = str(arcname or "").strip()

    if not arcname:
        raise DatasetDownloadArchiveError(
            "Archive entry name cannot be empty."
        )

    arc_path = Path(arcname)

    if arc_path.is_absolute():
        raise DatasetDownloadArchiveError(
            f"Archive entry name cannot be absolute: {arcname}"
        )

    if ".." in arc_path.parts:
        raise DatasetDownloadArchiveError(
            f"Archive entry name cannot contain '..': {arcname}"
        )


def sanitize_archive_filename_part(value: str) -> str:
    value = str(value or "").strip()

    if not value:
        return "dataset"

    value = value.replace("/", "_").replace("\\", "_").replace("..", "_")
    value = value.strip("._-")

    return value or "dataset"
