from __future__ import annotations

from pprint import pformat

from django.core.management.base import (
    BaseCommand,
    CommandError,
    CommandParser,
)

from analysis.services.recurrent_axis.rebuild_summaries import (
    DEFAULT_BATCH_SIZE,
    DEFAULT_SUMMARY_VERSION,
    rebuild_axis_final_recurrent_summary,
    rebuild_axis_structure_recurrent_summary,
    rebuild_recurrent_axis_summaries,
)


class Command(BaseCommand):
    help = (
        "Rebuild AxisStructureRecurrentSummary and "
        "AxisFinalRecurrentSummary from active recurrent-axis data."
    )

    def add_arguments(self, parser: CommandParser) -> None:
        parser.add_argument(
            "--only",
            choices=["all", "structure", "axis-final"],
            default="all",
            help=(
                "Select the summary to rebuild. "
                "Default: all."
            ),
        )
        parser.add_argument(
            "--batch-size",
            type=int,
            default=DEFAULT_BATCH_SIZE,
            help=(
                "bulk_create batch size. "
                f"Default: {DEFAULT_BATCH_SIZE}."
            ),
        )
        parser.add_argument(
            "--iterator-chunk-size",
            type=int,
            default=None,
            help=(
                "Database iterator chunk size. "
                "Defaults to --batch-size."
            ),
        )
        parser.add_argument(
            "--summary-version",
            type=int,
            default=DEFAULT_SUMMARY_VERSION,
            help=(
                "Version stored in structure summaries. "
                f"Default: {DEFAULT_SUMMARY_VERSION}."
            ),
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help=(
                "Inspect source counts without deleting or "
                "creating summary rows."
            ),
        )

    def handle(self, *args, **options) -> None:
        batch_size = options["batch_size"]
        iterator_chunk_size = (
            options["iterator_chunk_size"]
        )
        summary_version = options["summary_version"]
        dry_run = options["dry_run"]
        only = options["only"]

        if batch_size <= 0:
            raise CommandError(
                "--batch-size must be greater than zero."
            )

        if (
            iterator_chunk_size is not None
            and iterator_chunk_size <= 0
        ):
            raise CommandError(
                "--iterator-chunk-size must be greater "
                "than zero."
            )

        if summary_version <= 0:
            raise CommandError(
                "--summary-version must be greater than zero."
            )

        self.stdout.write(
            self.style.MIGRATE_HEADING(
                "===== recurrent summary rebuild ====="
            )
        )

        try:
            if only == "structure":
                result = (
                    rebuild_axis_structure_recurrent_summary(
                        batch_size=batch_size,
                        iterator_chunk_size=(
                            iterator_chunk_size
                        ),
                        summary_version=summary_version,
                        dry_run=dry_run,
                    )
                )

            elif only == "axis-final":
                result = rebuild_axis_final_recurrent_summary(
                    batch_size=batch_size,
                    iterator_chunk_size=(
                        iterator_chunk_size
                    ),
                    dry_run=dry_run,
                )

            else:
                result = rebuild_recurrent_axis_summaries(
                    batch_size=batch_size,
                    iterator_chunk_size=(
                        iterator_chunk_size
                    ),
                    summary_version=summary_version,
                    dry_run=dry_run,
                )

        except Exception as exc:
            raise CommandError(
                f"Recurrent summary rebuild failed: {exc}"
            ) from exc

        self.stdout.write(
            pformat(result, sort_dicts=False)
        )
        self.stdout.write("")

        if dry_run:
            self.stdout.write(
                self.style.SUCCESS(
                    "Summary dry run completed; no rows "
                    "were modified."
                )
            )
        else:
            self.stdout.write(
                self.style.SUCCESS(
                    "Recurrent summaries rebuilt successfully."
                )
            )
