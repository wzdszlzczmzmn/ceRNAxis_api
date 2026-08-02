from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from pathlib import Path
from typing import Any

from database.models import (
    AxisDatasetSource,
    AxisGroupType,
    AxisModule,
    AxisResultKind,
)
from database.utils.dataset_annotation_utils.path_utils import (
    get_timedb_json_group_by_fields, TIMEDB_IGNORED_JSON_GROUP_BY_FIELDS,
)

from .import_context import validate_context_spec
from .import_contracts import (
    AxisContextSpec,
    AxisImportJob,
)
from .import_normalization import (
    AxisImportValidationError,
)
from .import_service import import_axis_result_file


DEFAULT_RESULT_FILE_TEMPLATES = {
    AxisResultKind.AXIS_FINAL:
        "{prefix}_ceRNA_axis_final.csv",

    AxisResultKind.SPONGE:
        "{prefix}_sponge_result.csv",
}


OtherGroupByResolver = Callable[
    [str, Path],
    str | None,
]


def iter_result_dirs(
    root_dir: str | Path,
) -> list[Path]:
    """
    Return visible direct child directories in deterministic order.

    Discovery deliberately does not recurse. Both Module2 and Module3 use one
    directory level as part of the dataset/context identity.
    """
    root_dir = Path(root_dir).expanduser().resolve()

    if not root_dir.exists():
        raise AxisImportValidationError(
            f"Result root directory does not exist: {root_dir}"
        )

    if not root_dir.is_dir():
        raise AxisImportValidationError(
            f"Result root path is not a directory: {root_dir}"
        )

    return sorted(
        (
            item.resolve()
            for item in root_dir.iterdir()
            if (
                item.is_dir()
                and not item.name.startswith(".")
            )
        ),
        key=lambda item: item.name.casefold(),
    )


def build_module2_import_jobs(
    *,
    module2_root_dir: str | Path,
    file_templates: Mapping[str, str] | None = None,
    schema_versions: Mapping[str, str] | None = None,
    result_kinds: Iterable[str] | None = None,
) -> tuple[list[AxisImportJob], list[dict]]:
    """
    Discover Module2 result files.

    Directory convention:
        <root>/TCGA_ACC/

    Context convention:
        dataset_name = TCGA_ACC_mRNA
        group_type = none
        group_by = ""

    Default result files:
        TCGA_ACC_ceRNA_axis_final.csv
        TCGA_ACC_sponge_result.csv
    """
    templates = _build_file_templates(
        file_templates=file_templates,
        result_kinds=result_kinds,
    )
    versions = _build_schema_versions(
        schema_versions=schema_versions,
        result_kinds=templates,
    )

    jobs: list[AxisImportJob] = []
    missing_files: list[dict] = []

    for project_dir in iter_result_dirs(
        module2_root_dir
    ):
        prefix = project_dir.name

        context = AxisContextSpec(
            dataset_source=AxisDatasetSource.TCGA,
            module=AxisModule.MODULE2,
            dataset_name=f"{prefix}_mRNA",
            group_type=AxisGroupType.NONE,
            group_by="",
            annotation_dir_name=prefix,
            annotation_file_prefix=prefix,
        )
        validate_context_spec(context)

        _append_context_jobs(
            jobs=jobs,
            missing_files=missing_files,
            context=context,
            result_dir=project_dir,
            templates=templates,
            schema_versions=versions,
        )

    return jobs, missing_files


def build_module3_import_jobs(
    *,
    module3_root_dir: str | Path,
    file_templates: Mapping[str, str] | None = None,
    schema_versions: Mapping[str, str] | None = None,
    result_kinds: Iterable[str] | None = None,
    other_group_by_resolver: (
        OtherGroupByResolver | None
    ) = None,
) -> tuple[list[AxisImportJob], list[dict]]:
    """
    Discover Module3 result files.

    Directory conventions:
        GSE20194
            group_type = other
            group_by = value resolved from TIMEDB JSON

        GSE20194_grade
            group_type = grade
            group_by = Grade

        GSE20194_stage
            group_type = stage
            group_by = Stage

    For all three directory forms, the file prefix is the base dataset name:
        GSE20194_ceRNA_axis_final.csv
        GSE20194_sponge_result.csv
    """
    templates = _build_file_templates(
        file_templates=file_templates,
        result_kinds=result_kinds,
    )
    versions = _build_schema_versions(
        schema_versions=schema_versions,
        result_kinds=templates,
    )
    resolver = (
        other_group_by_resolver
        or resolve_other_group_by
    )

    jobs: list[AxisImportJob] = []
    missing_files: list[dict] = []

    for annotation_dir in iter_result_dirs(
            module3_root_dir
    ):
        context = parse_module3_context(
            annotation_dir=annotation_dir,
            other_group_by_resolver=resolver,
        )

        if context is None:
            continue

        validate_context_spec(context)

        _append_context_jobs(
            jobs=jobs,
            missing_files=missing_files,
            context=context,
            result_dir=annotation_dir,
            templates=templates,
            schema_versions=versions,
        )

    return jobs, missing_files


def resolve_other_group_by(
    dataset_name: str,
    annotation_dir: Path,
) -> str | None:
    """
    Resolve the importable group_by field represented by a Module3
    base annotation directory.

    JSON fields representing tumor grade or tumor stage are ignored because
    they are represented by dedicated ``_grade`` and ``_stage`` directories.

    Returns:
        str:
            The one valid ``other`` group_by field.

        None:
            The JSON contains group_by fields, but all of them are ignored
            grade/stage fields. The base directory should not produce an
            Other Context.

    Raises:
        AxisImportValidationError:
            No usable metadata exists, or more than one importable Other
            group_by field remains.
    """
    raw_fields = get_timedb_json_group_by_fields(
        dataset_name
    )

    non_empty_fields: list[str] = []
    importable_fields: list[str] = []

    ignored_fields = {
        str(field).strip().casefold()
        for field in TIMEDB_IGNORED_JSON_GROUP_BY_FIELDS
    }

    for value in raw_fields or []:
        field = str(value or "").strip()

        if not field:
            continue

        # Deduplicate case-insensitively while preserving the original value.
        if not any(
            existing.casefold() == field.casefold()
            for existing in non_empty_fields
        ):
            non_empty_fields.append(field)

        if field.casefold() in ignored_fields:
            continue

        if not any(
            existing.casefold() == field.casefold()
            for existing in importable_fields
        ):
            importable_fields.append(field)

    if not non_empty_fields:
        raise AxisImportValidationError(
            f"{annotation_dir} does not define a non-empty "
            f"group_by field for dataset={dataset_name!r}."
        )

    if not importable_fields:
        # The JSON only contains c_tumor_grade and/or c_tumor_stage.
        # Dedicated grade/stage directories represent these contexts.
        return None

    if len(importable_fields) != 1:
        raise AxisImportValidationError(
            f"{annotation_dir} resolves to "
            f"{len(importable_fields)} importable Other group_by fields "
            f"for dataset={dataset_name!r}: "
            f"{importable_fields!r}; exactly one is required."
        )

    return importable_fields[0]


def parse_module3_context(
    *,
    annotation_dir: str | Path,
    other_group_by_resolver: (
        OtherGroupByResolver | None
    ) = None,
) -> AxisContextSpec | None:
    annotation_dir = (
        Path(annotation_dir)
        .expanduser()
        .resolve()
    )

    if not annotation_dir.is_dir():
        raise AxisImportValidationError(
            "Module3 annotation path is not a directory: "
            f"{annotation_dir}"
        )

    directory_name = annotation_dir.name

    if directory_name.endswith("_grade"):
        dataset_name = directory_name[
            :-len("_grade")
        ]
        group_type = AxisGroupType.GRADE
        group_by = "Grade"

    elif directory_name.endswith("_stage"):
        dataset_name = directory_name[
            :-len("_stage")
        ]
        group_type = AxisGroupType.STAGE
        group_by = "Stage"



    else:
        dataset_name = directory_name
        group_type = AxisGroupType.OTHER

        resolver = (
                other_group_by_resolver
                or resolve_other_group_by
        )

        group_by = resolver(
            dataset_name,
            annotation_dir,
        )

        if group_by is None:
            return None

    if not dataset_name:
        raise AxisImportValidationError(
            "Unable to derive dataset_name from Module3 "
            f"directory: {annotation_dir}"
        )

    group_by = str(group_by or "").strip()

    context = AxisContextSpec(
        dataset_source=AxisDatasetSource.TIMEDB,
        module=AxisModule.MODULE3,
        dataset_name=dataset_name,
        group_type=group_type,
        group_by=group_by,
        annotation_dir_name=directory_name,
        annotation_file_prefix=dataset_name,
    )

    validate_context_spec(context)

    return context


def import_axis_jobs(
    *,
    jobs: Iterable[AxisImportJob],
    dry_run: bool = False,
    skip_unchanged: bool = True,
    batch_size: int = 1000,
    stop_on_error: bool = False,
) -> dict:
    """
    Execute discovery jobs one artifact at a time.

    Each import_axis_result_file() call owns its own transaction, so a failure
    does not roll back previously completed artifacts.
    """
    jobs = list(jobs)
    results: list[dict] = []

    for job_index, job in enumerate(jobs):
        try:
            result = import_axis_result_file(
                context_spec=job.context,
                result_kind=job.result_kind,
                file_path=job.file_path,
                schema_version=job.schema_version,
                dry_run=dry_run,
                skip_unchanged=skip_unchanged,
                batch_size=batch_size,
            )
        except Exception as exc:
            result = {
                "success": False,
                "dry_run": dry_run,
                "skipped": False,
                "job_index": job_index,
                "dataset_source":
                    job.context.dataset_source,
                "module":
                    job.context.module,
                "dataset_name":
                    job.context.dataset_name,
                "group_type":
                    job.context.group_type,
                "group_by":
                    job.context.group_by,
                "result_kind":
                    job.result_kind,
                "schema_version":
                    job.schema_version,
                "file_path":
                    str(job.file_path),
                "error_type":
                    type(exc).__name__,
                "error":
                    str(exc),
            }

            results.append(result)

            if stop_on_error:
                raise

            continue

        result = {
            "job_index": job_index,
            **result,
        }
        results.append(result)

    imported_count = sum(
        1
        for item in results
        if (
            item.get("success")
            and not item.get("skipped")
            and not item.get("dry_run")
        )
    )
    dry_run_count = sum(
        1
        for item in results
        if (
            item.get("success")
            and item.get("dry_run")
        )
    )
    skipped_count = sum(
        1
        for item in results
        if (
            item.get("success")
            and item.get("skipped")
        )
    )
    failed_count = sum(
        1
        for item in results
        if not item.get("success")
    )

    return {
        "success": failed_count == 0,
        "job_count": len(jobs),
        "imported_count": imported_count,
        "dry_run_count": dry_run_count,
        "skipped_count": skipped_count,
        "failed_count": failed_count,
        "results": results,
    }


def _append_context_jobs(
    *,
    jobs: list[AxisImportJob],
    missing_files: list[dict],
    context: AxisContextSpec,
    result_dir: Path,
    templates: Mapping[str, str],
    schema_versions: Mapping[str, str],
) -> None:
    for result_kind, template in templates.items():
        file_path = _render_result_file_path(
            result_dir=result_dir,
            template=template,
            context=context,
            result_kind=result_kind,
        )

        if file_path.is_file():
            jobs.append(
                AxisImportJob(
                    context=context,
                    result_kind=result_kind,
                    file_path=file_path,
                    schema_version=(
                        schema_versions[result_kind]
                    ),
                )
            )
            continue

        missing_files.append({
            "dataset_source":
                context.dataset_source,
            "module":
                context.module,
            "dataset_name":
                context.dataset_name,
            "group_type":
                context.group_type,
            "group_by":
                context.group_by,
            "annotation_dir_name":
                context.annotation_dir_name,
            "result_kind":
                result_kind,
            "schema_version":
                schema_versions[result_kind],
            "file_path":
                str(file_path),
            "reason":
                "result_file_not_found",
        })


def _render_result_file_path(
    *,
    result_dir: Path,
    template: str,
    context: AxisContextSpec,
    result_kind: str,
) -> Path:
    values: dict[str, Any] = {
        "prefix":
            context.annotation_file_prefix,
        "dataset_name":
            context.dataset_name,
        "annotation_dir_name":
            context.annotation_dir_name,
        "group_type":
            context.group_type,
        "group_by":
            context.group_by,
        "result_kind":
            result_kind,
    }

    try:
        file_name = template.format(**values)
    except KeyError as exc:
        raise AxisImportValidationError(
            f"Unknown placeholder {exc.args[0]!r} in "
            f"file template {template!r}."
        ) from exc
    except (IndexError, ValueError) as exc:
        raise AxisImportValidationError(
            f"Invalid file template: {template!r}."
        ) from exc

    file_name = str(file_name).strip()

    if not file_name:
        raise AxisImportValidationError(
            "A rendered result file name cannot be empty."
        )

    candidate = Path(file_name)

    if candidate.is_absolute():
        raise AxisImportValidationError(
            "Result file templates must render relative "
            f"paths, got: {candidate}"
        )

    result_dir = result_dir.resolve()
    file_path = (
        result_dir / candidate
    ).resolve()

    try:
        file_path.relative_to(result_dir)
    except ValueError as exc:
        raise AxisImportValidationError(
            "Result file template escapes its result "
            f"directory: {template!r}."
        ) from exc

    return file_path


def _build_file_templates(
    *,
    file_templates: Mapping[str, str] | None,
    result_kinds: Iterable[str] | None,
) -> dict[str, str]:
    templates = dict(
        DEFAULT_RESULT_FILE_TEMPLATES
    )

    if file_templates:
        templates.update(file_templates)

    if result_kinds is None:
        selected_result_kinds = tuple(
            templates.keys()
        )
    else:
        selected_result_kinds = tuple(
            dict.fromkeys(result_kinds)
        )

    if not selected_result_kinds:
        raise AxisImportValidationError(
            "At least one result kind must be selected."
        )

    supported_result_kinds = {
        AxisResultKind.AXIS_FINAL,
        AxisResultKind.SPONGE,
    }

    unknown_result_kinds = (
        set(selected_result_kinds)
        - supported_result_kinds
    )

    if unknown_result_kinds:
        raise AxisImportValidationError(
            "Unsupported result kinds: "
            f"{sorted(unknown_result_kinds)!r}."
        )

    selected_templates: dict[str, str] = {}

    for result_kind in selected_result_kinds:
        template = templates.get(result_kind)

        if not isinstance(template, str):
            raise AxisImportValidationError(
                "Missing string file template for "
                f"result_kind={result_kind!r}."
            )

        template = template.strip()

        if not template:
            raise AxisImportValidationError(
                "File template cannot be empty for "
                f"result_kind={result_kind!r}."
            )

        selected_templates[result_kind] = template

    return selected_templates


def _build_schema_versions(
    *,
    schema_versions: Mapping[str, str] | None,
    result_kinds: Iterable[str],
) -> dict[str, str]:
    schema_versions = dict(
        schema_versions or {}
    )

    versions: dict[str, str] = {}

    for result_kind in result_kinds:
        version = str(
            schema_versions.get(
                result_kind,
                "v1",
            )
            or ""
        ).strip()

        if not version:
            raise AxisImportValidationError(
                "schema_version cannot be empty for "
                f"result_kind={result_kind!r}."
            )

        versions[result_kind] = version

    return versions
