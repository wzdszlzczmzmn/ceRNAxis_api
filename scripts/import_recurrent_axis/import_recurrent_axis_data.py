from pathlib import Path
from pprint import pprint

from database.models import (
    AxisResultArtifact,
    AxisResultKind,
    DatasetMetadata,
)

from analysis.services.recurrent_axis.import_discovery import (
    build_module2_import_jobs,
    build_module3_import_jobs,
    import_axis_jobs,
)


# ---------------------------------------------------------------------
# 修改为实际路径
# ---------------------------------------------------------------------

MODULE2_ROOT = Path(
    "E:\\Projects\\ceRNAxis\\data\\module2"
)

MODULE3_ROOT = Path(
    "E:\\Projects\\ceRNAxis\\data\\module3"
)


# ---------------------------------------------------------------------
# 导入配置
# ---------------------------------------------------------------------

RESULT_KINDS = [
    AxisResultKind.AXIS_FINAL,
    AxisResultKind.SPONGE,
]

FILE_TEMPLATES = {
    AxisResultKind.AXIS_FINAL:
        "{prefix}_ceRNA_axis_final.csv",

    AxisResultKind.SPONGE:
        "{prefix}_sponge_result.csv",
}

SCHEMA_VERSIONS = {
    AxisResultKind.AXIS_FINAL: "v1",
    AxisResultKind.SPONGE: "v1",
}

BATCH_SIZE = 1000

# True：发现任何缺失文件时停止。
# False：只导入存在的文件。
STRICT_MISSING_FILES = False


def print_missing_files(
    *,
    title: str,
    missing_files: list[dict],
) -> None:
    print()
    print(f"===== {title}: missing files =====")
    print(f"count: {len(missing_files)}")

    for item in missing_files:
        print(
            f"- dataset={item['dataset_name']}, "
            f"group_type={item['group_type']}, "
            f"group_by={item['group_by']!r}, "
            f"result_kind={item['result_kind']}, "
            f"path={item['file_path']}"
        )


def validate_dataset_metadata(jobs) -> None:
    required_datasets = {
        job.context.dataset_name
        for job in jobs
    }

    existing_datasets = set(
        DatasetMetadata.objects
        .filter(dataset__in=required_datasets)
        .values_list("dataset", flat=True)
    )

    missing_datasets = sorted(
        required_datasets - existing_datasets
    )

    if missing_datasets:
        raise RuntimeError(
            "DatasetMetadata is missing required datasets:\n"
            + "\n".join(
                f"- {dataset_name}"
                for dataset_name in missing_datasets
            )
        )


def print_failed_jobs(result: dict) -> None:
    failed = [
        item
        for item in result["results"]
        if not item.get("success")
    ]

    if not failed:
        return

    print()
    print("===== failed jobs =====")

    for item in failed:
        pprint(item)


def main() -> None:
    # -------------------------------------------------------------
    # 1. Discover Module2
    # -------------------------------------------------------------

    module2_jobs, module2_missing = (
        build_module2_import_jobs(
            module2_root_dir=MODULE2_ROOT,
            file_templates=FILE_TEMPLATES,
            schema_versions=SCHEMA_VERSIONS,
            result_kinds=RESULT_KINDS,
        )
    )

    # -------------------------------------------------------------
    # 2. Discover Module3
    # -------------------------------------------------------------

    module3_jobs, module3_missing = (
        build_module3_import_jobs(
            module3_root_dir=MODULE3_ROOT,
            file_templates=FILE_TEMPLATES,
            schema_versions=SCHEMA_VERSIONS,
            result_kinds=RESULT_KINDS,
        )
    )

    all_jobs = [
        *module2_jobs,
        *module3_jobs,
    ]

    all_missing = [
        *module2_missing,
        *module3_missing,
    ]

    print("===== discovery =====")
    print(f"Module2 jobs: {len(module2_jobs)}")
    print(f"Module3 jobs: {len(module3_jobs)}")
    print(f"Total jobs:   {len(all_jobs)}")

    print_missing_files(
        title="Module2",
        missing_files=module2_missing,
    )

    print_missing_files(
        title="Module3",
        missing_files=module3_missing,
    )

    if STRICT_MISSING_FILES and all_missing:
        raise RuntimeError(
            f"{len(all_missing)} expected result files "
            "are missing."
        )

    if not all_jobs:
        raise RuntimeError(
            "No import jobs were discovered."
        )

    # -------------------------------------------------------------
    # 3. Validate DatasetMetadata before parsing every file
    # -------------------------------------------------------------

    validate_dataset_metadata(all_jobs)

    # -------------------------------------------------------------
    # 4. Dry run
    # -------------------------------------------------------------

    print()
    print("===== dry run =====")

    dry_run_result = import_axis_jobs(
        jobs=all_jobs,
        dry_run=True,
        skip_unchanged=True,
        batch_size=BATCH_SIZE,
        stop_on_error=False,
    )

    pprint({
        key: value
        for key, value in dry_run_result.items()
        if key != "results"
    })

    print_failed_jobs(dry_run_result)

    if dry_run_result["failed_count"] > 0:
        raise RuntimeError(
            "Dry run failed. No database import was started."
        )

    # -------------------------------------------------------------
    # 5. Real import
    # -------------------------------------------------------------

    print()
    print("===== real import =====")

    import_result = import_axis_jobs(
        jobs=all_jobs,
        dry_run=False,
        skip_unchanged=True,
        batch_size=BATCH_SIZE,
        stop_on_error=False,
    )

    pprint({
        key: value
        for key, value in import_result.items()
        if key != "results"
    })

    print_failed_jobs(import_result)

    # -------------------------------------------------------------
    # 6. Basic report
    # -------------------------------------------------------------

    print()
    print("===== active artifacts =====")

    active_artifacts = (
        AxisResultArtifact.objects
        .filter(is_active=True)
        .values(
            "context__dataset_source",
            "context__module",
            "context__dataset_metadata_id",
            "context__group_type",
            "context__group_by",
            "result_kind",
            "row_count",
        )
        .order_by(
            "context__dataset_source",
            "context__dataset_metadata_id",
            "context__group_type",
            "result_kind",
        )
    )

    for artifact in active_artifacts.iterator():
        pprint(artifact)

    if import_result["failed_count"] > 0:
        raise RuntimeError(
            f"{import_result['failed_count']} import jobs failed."
        )

    print()
    print("Import completed successfully.")


main()
