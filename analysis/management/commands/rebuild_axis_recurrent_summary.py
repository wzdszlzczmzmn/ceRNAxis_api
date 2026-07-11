from django.core.management.base import (
    BaseCommand,
    CommandError,
)

from analysis.services.axis_final import (
    rebuild_axis_recurrent_summary,
)


class Command(BaseCommand):
    help = (
        "Rebuild AxisRecurrentSummary from "
        "AxisSignatureProjectIndex."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--batch-size",
            type=int,
            default=5000,
            help=(
                "bulk_create batch size. "
                "Default: 5000."
            ),
        )

        parser.add_argument(
            "--iterator-chunk-size",
            type=int,
            default=None,
            help=(
                "Database iterator chunk size. "
                "Defaults to batch-size."
            ),
        )

        parser.add_argument(
            "--dry-run",
            action="store_true",
            help=(
                "Count source rows and signatures without "
                "modifying AxisRecurrentSummary."
            ),
        )

        parser.add_argument(
            "--include-inactive-projects",
            action="store_true",
            help=(
                "Include AxisSignatureProjectIndex rows "
                "belonging to inactive projects."
            ),
        )

    def handle(self, *args, **options):
        try:
            result = rebuild_axis_recurrent_summary(
                batch_size=options["batch_size"],
                iterator_chunk_size=(
                    options["iterator_chunk_size"]
                ),
                active_projects_only=not options[
                    "include_inactive_projects"
                ],
                dry_run=options["dry_run"],
            )
        except Exception as exc:
            raise CommandError(str(exc)) from exc

        if result["dry_run"]:
            self.stdout.write(
                self.style.WARNING(
                    "Axis recurrent summary dry run finished."
                )
            )
        else:
            self.stdout.write(
                self.style.SUCCESS(
                    "Axis recurrent summary rebuild finished."
                )
            )

        self.stdout.write(
            f"dry_run: {result['dry_run']}"
        )
        self.stdout.write(
            "source_index_row_count: "
            f"{result['source_index_row_count']}"
        )
        self.stdout.write(
            "source_axis_signature_count: "
            f"{result['source_axis_signature_count']}"
        )
        self.stdout.write(
            "deleted_summary_count: "
            f"{result['deleted_summary_count']}"
        )
        self.stdout.write(
            "created_summary_count: "
            f"{result['created_summary_count']}"
        )
        self.stdout.write(
            "final_summary_count: "
            f"{result['final_summary_count']}"
        )
        self.stdout.write(
            f"batch_size: {result['batch_size']}"
        )
        self.stdout.write(
            "iterator_chunk_size: "
            f"{result['iterator_chunk_size']}"
        )
        self.stdout.write(
            "active_projects_only: "
            f"{result['active_projects_only']}"
        )
