from __future__ import annotations

from pathlib import Path
from pprint import pformat

from django.core.management.base import BaseCommand, CommandError, CommandParser

from analysis.services.recurrent_axis.import_discovery import (
    build_module2_import_jobs,
    build_module3_import_jobs,
    import_axis_jobs,
)
from analysis.services.recurrent_axis.rebuild_summaries import (
    DEFAULT_SUMMARY_VERSION,
    rebuild_recurrent_axis_summaries,
)
from database.models import (
    AxisResultArtifact,
    AxisResultKind,
    DatasetMetadata,
)


DEFAULT_AXIS_FINAL_TEMPLATE = "{prefix}_ceRNA_axis_final.csv"
DEFAULT_SPONGE_TEMPLATE = "{prefix}_sponge_result.csv"
DEFAULT_SCHEMA_VERSION = "v1"
DEFAULT_BATCH_SIZE = 1000


class Command(BaseCommand):
    help = (
        "Discover, validate, and import recurrent Axis results from Module2 "
        "and Module3 result directories."
    )

    def add_arguments(self, parser: CommandParser) -> None:
        parser.add_argument(
            "--module2-root",
            required=True,
            type=Path,
            help=(
                "Module2 result root directory. Its direct child directories "
                "must be dataset/project directories such as TCGA_ACC."
            ),
        )
        parser.add_argument(
            "--module3-root",
            required=True,
            type=Path,
            help=(
                "Module3 result root directory. Its direct child directories "
                "must be TIMEDB annotation directories such as GSE20194, "
                "GSE20194_grade, or GSE20194_stage."
            ),
        )
        parser.add_argument(
            "--dry-run-only",
            action="store_true",
            help=(
                "Run discovery and full file validation without writing to "
                "the database. By default the command performs a dry run and "
                "then imports when validation succeeds."
            ),
        )
        parser.add_argument(
            "--strict-missing-files",
            action="store_true",
            help=(
                "Fail when any expected result file is missing. Without this "
                "option, existing files are imported and missing files are "
                "reported only."
            ),
        )
        parser.add_argument(
            "--batch-size",
            type=int,
            default=DEFAULT_BATCH_SIZE,
            help=f"Database bulk operation batch size. Default: {DEFAULT_BATCH_SIZE}.",
        )
        parser.add_argument(
            "--result-kind",
            action="append",
            choices=[
                AxisResultKind.AXIS_FINAL,
                AxisResultKind.SPONGE,
            ],
            dest="result_kinds",
            help=(
                "Result kind to import. Repeat this option to import multiple "
                "kinds. The default imports both axis_final and sponge."
            ),
        )
        parser.add_argument(
            "--no-skip-unchanged",
            action="store_true",
            help=(
                "Re-import an artifact even when its hash, schema version, "
                "row count, observations, and evidence are unchanged."
            ),
        )
        parser.add_argument(
            "--stop-on-error",
            action="store_true",
            help=(
                "Stop at the first failed import job. By default all jobs are "
                "attempted and failures are summarized at the end."
            ),
        )
        parser.add_argument(
            "--axis-final-template",
            default=DEFAULT_AXIS_FINAL_TEMPLATE,
            help=(
                "Axis Final filename template relative to each result "
                f"directory. Default: {DEFAULT_AXIS_FINAL_TEMPLATE!r}."
            ),
        )
        parser.add_argument(
            "--sponge-template",
            default=DEFAULT_SPONGE_TEMPLATE,
            help=(
                "Sponge filename template relative to each result directory. "
                f"Default: {DEFAULT_SPONGE_TEMPLATE!r}."
            ),
        )
        parser.add_argument(
            "--axis-final-schema-version",
            default=DEFAULT_SCHEMA_VERSION,
            help=(
                "Schema version stored for Axis Final artifacts. "
                f"Default: {DEFAULT_SCHEMA_VERSION!r}."
            ),
        )
        parser.add_argument(
            "--sponge-schema-version",
            default=DEFAULT_SCHEMA_VERSION,
            help=(
                "Schema version stored for Sponge artifacts. "
                f"Default: {DEFAULT_SCHEMA_VERSION!r}."
            ),
        )
        parser.add_argument(
            "--skip-summary-rebuild",
            action="store_true",
            help=(
                "Do not rebuild AxisStructureRecurrentSummary and "
                "AxisFinalRecurrentSummary after a successful import."
            ),
        )
        parser.add_argument(
            "--summary-version",
            type=int,
            default=DEFAULT_SUMMARY_VERSION,
            help=(
                "Version stored in AxisStructureRecurrentSummary. "
                f"Default: {DEFAULT_SUMMARY_VERSION}."
            ),
        )

    def handle(self, *args, **options) -> None:
        module2_root = self._resolve_root_path(
            options["module2_root"],
            option_name="--module2-root",
        )
        module3_root = self._resolve_root_path(
            options["module3_root"],
            option_name="--module3-root",
        )

        batch_size = options["batch_size"]
        if batch_size <= 0:
            raise CommandError("--batch-size must be greater than zero.")

        summary_version = options["summary_version"]
        if summary_version <= 0:
            raise CommandError(
                "--summary-version must be greater than zero."
            )

        result_kinds = options["result_kinds"] or [
            AxisResultKind.AXIS_FINAL,
            AxisResultKind.SPONGE,
        ]
        # Preserve command-line order while removing duplicates.
        result_kinds = list(dict.fromkeys(result_kinds))

        file_templates = {
            AxisResultKind.AXIS_FINAL:
                options["axis_final_template"],
            AxisResultKind.SPONGE:
                options["sponge_template"],
        }
        schema_versions = {
            AxisResultKind.AXIS_FINAL:
                options["axis_final_schema_version"],
            AxisResultKind.SPONGE:
                options["sponge_schema_version"],
        }

        module2_jobs, module2_missing = self._discover_module2(
            root_dir=module2_root,
            file_templates=file_templates,
            schema_versions=schema_versions,
            result_kinds=result_kinds,
        )
        module3_jobs, module3_missing = self._discover_module3(
            root_dir=module3_root,
            file_templates=file_templates,
            schema_versions=schema_versions,
            result_kinds=result_kinds,
        )

        all_jobs = [
            *module2_jobs,
            *module3_jobs,
        ]
        all_missing = [
            *module2_missing,
            *module3_missing,
        ]

        self._print_discovery_summary(
            module2_jobs=module2_jobs,
            module3_jobs=module3_jobs,
            module2_missing=module2_missing,
            module3_missing=module3_missing,
        )

        if options["strict_missing_files"] and all_missing:
            raise CommandError(
                f"{len(all_missing)} expected result files are missing."
            )

        if not all_jobs:
            raise CommandError("No import jobs were discovered.")

        self._validate_dataset_metadata(all_jobs)

        skip_unchanged = not options["no_skip_unchanged"]
        stop_on_error = options["stop_on_error"]

        self.stdout.write("")
        self.stdout.write(self.style.MIGRATE_HEADING("===== dry run ====="))

        dry_run_result = import_axis_jobs(
            jobs=all_jobs,
            dry_run=True,
            skip_unchanged=skip_unchanged,
            batch_size=batch_size,
            stop_on_error=stop_on_error,
        )
        self._print_import_summary(dry_run_result)
        self._print_failed_jobs(dry_run_result)

        if dry_run_result["failed_count"] > 0:
            raise CommandError(
                "Dry run failed. No database import was started."
            )

        if options["dry_run_only"]:
            self.stdout.write("")
            self.stdout.write(
                self.style.SUCCESS(
                    "Dry run completed successfully. No database rows were written."
                )
            )
            return

        self.stdout.write("")
        self.stdout.write(self.style.MIGRATE_HEADING("===== real import ====="))

        import_result = import_axis_jobs(
            jobs=all_jobs,
            dry_run=False,
            skip_unchanged=skip_unchanged,
            batch_size=batch_size,
            stop_on_error=stop_on_error,
        )
        self._print_import_summary(import_result)
        self._print_failed_jobs(import_result)
        self._print_active_artifacts(import_result)

        if import_result["failed_count"] > 0:
            raise CommandError(
                f"{import_result['failed_count']} import jobs failed."
            )

        if not options["skip_summary_rebuild"]:
            self.stdout.write("")
            self.stdout.write(
                self.style.MIGRATE_HEADING(
                    "===== recurrent summary rebuild ====="
                )
            )

            try:
                summary_result = rebuild_recurrent_axis_summaries(
                    batch_size=batch_size,
                    iterator_chunk_size=batch_size,
                    summary_version=summary_version,
                    dry_run=False,
                )
            except Exception as exc:
                raise CommandError(
                    f"Import succeeded, but recurrent summary "
                    f"rebuild failed: {exc}"
                ) from exc

            self.stdout.write(
                pformat(summary_result, sort_dicts=False)
            )

        self.stdout.write("")
        self.stdout.write(
            self.style.SUCCESS(
                "Recurrent Axis import completed successfully."
            )
        )

        if options["skip_summary_rebuild"]:
            self.stdout.write(
                self.style.WARNING(
                    "Summary rebuild was skipped. Run "
                    "rebuild_recurrent_axis_summaries before "
                    "serving recurrent-axis summary queries."
                )
            )

    def _discover_module2(
        self,
        *,
        root_dir: Path,
        file_templates: dict[str, str],
        schema_versions: dict[str, str],
        result_kinds: list[str],
    ):
        try:
            return build_module2_import_jobs(
                module2_root_dir=root_dir,
                file_templates=file_templates,
                schema_versions=schema_versions,
                result_kinds=result_kinds,
            )
        except Exception as exc:
            raise CommandError(
                f"Module2 discovery failed: {exc}"
            ) from exc

    def _discover_module3(
        self,
        *,
        root_dir: Path,
        file_templates: dict[str, str],
        schema_versions: dict[str, str],
        result_kinds: list[str],
    ):
        try:
            return build_module3_import_jobs(
                module3_root_dir=root_dir,
                file_templates=file_templates,
                schema_versions=schema_versions,
                result_kinds=result_kinds,
            )
        except Exception as exc:
            raise CommandError(
                f"Module3 discovery failed: {exc}"
            ) from exc

    @staticmethod
    def _resolve_root_path(
        path: Path,
        *,
        option_name: str,
    ) -> Path:
        resolved = path.expanduser().resolve()

        if not resolved.exists():
            raise CommandError(
                f"{option_name} does not exist: {resolved}"
            )

        if not resolved.is_dir():
            raise CommandError(
                f"{option_name} is not a directory: {resolved}"
            )

        return resolved

    def _validate_dataset_metadata(self, jobs) -> None:
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

        if not missing_datasets:
            self.stdout.write(
                self.style.SUCCESS(
                    "DatasetMetadata validation passed for "
                    f"{len(required_datasets)} datasets."
                )
            )
            return

        details = "\n".join(
            f"- {dataset_name}"
            for dataset_name in missing_datasets
        )
        raise CommandError(
            "DatasetMetadata is missing required datasets:\n"
            f"{details}"
        )

    def _print_discovery_summary(
        self,
        *,
        module2_jobs,
        module3_jobs,
        module2_missing,
        module3_missing,
    ) -> None:
        self.stdout.write(
            self.style.MIGRATE_HEADING("===== discovery =====")
        )
        self.stdout.write(f"Module2 jobs: {len(module2_jobs)}")
        self.stdout.write(f"Module3 jobs: {len(module3_jobs)}")
        self.stdout.write(
            f"Total jobs:   {len(module2_jobs) + len(module3_jobs)}"
        )

        self._print_missing_files(
            title="Module2",
            missing_files=module2_missing,
        )
        self._print_missing_files(
            title="Module3",
            missing_files=module3_missing,
        )

    def _print_missing_files(
        self,
        *,
        title: str,
        missing_files: list[dict],
    ) -> None:
        self.stdout.write("")
        self.stdout.write(
            self.style.MIGRATE_LABEL(
                f"===== {title}: missing files ====="
            )
        )
        self.stdout.write(f"count: {len(missing_files)}")

        for item in missing_files:
            self.stdout.write(
                "- "
                f"dataset={item['dataset_name']}, "
                f"group_type={item['group_type']}, "
                f"group_by={item['group_by']!r}, "
                f"result_kind={item['result_kind']}, "
                f"path={item['file_path']}"
            )

    def _print_import_summary(self, result: dict) -> None:
        summary = {
            key: value
            for key, value in result.items()
            if key != "results"
        }
        self.stdout.write(pformat(summary, sort_dicts=False))

    def _print_failed_jobs(self, result: dict) -> None:
        failed = [
            item
            for item in result["results"]
            if not item.get("success")
        ]

        if not failed:
            return

        self.stdout.write("")
        self.stdout.write(self.style.ERROR("===== failed jobs ====="))

        for item in failed:
            self.stdout.write(pformat(item, sort_dicts=False))

    def _print_active_artifacts(self, import_result: dict) -> None:
        context_ids = {
            item.get("context_id")
            for item in import_result["results"]
            if (
                item.get("success")
                and item.get("context_id") is not None
            )
        }

        self.stdout.write("")
        self.stdout.write(
            self.style.MIGRATE_HEADING(
                "===== active artifacts for imported contexts ====="
            )
        )

        if not context_ids:
            self.stdout.write("No imported context IDs were returned.")
            return

        artifacts = (
            AxisResultArtifact.objects
            .filter(
                is_active=True,
                context_id__in=context_ids,
            )
            .values(
                "id",
                "context__dataset_source",
                "context__module",
                "context__dataset_metadata_id",
                "context__group_type",
                "context__group_by",
                "result_kind",
                "row_count",
                "file_name",
            )
            .order_by(
                "context__dataset_source",
                "context__dataset_metadata_id",
                "context__group_type",
                "context__group_by",
                "result_kind",
            )
        )

        count = 0
        for artifact in artifacts.iterator():
            self.stdout.write(pformat(artifact, sort_dicts=False))
            count += 1

        self.stdout.write(f"active artifact count: {count}")
