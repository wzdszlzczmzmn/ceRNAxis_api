from django.core.management.base import BaseCommand

from analysis.services.axis_final import rebuild_axis_signature_project_index


class Command(BaseCommand):
    help = (
        "Rebuild AxisSignatureProjectIndex from "
        "DatasetAxisFinalOccurrence."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--batch-size",
            type=int,
            default=5000,
            help="bulk_create batch size. Default: 5000.",
        )

        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Only count source rows without writing database.",
        )

        parser.add_argument(
            "--include-inactive-projects",
            action="store_true",
            help="Also index occurrences from inactive projects.",
        )

        parser.add_argument(
            "--no-clear",
            action="store_true",
            help=(
                "Do not clear existing AxisSignatureProjectIndex rows "
                "before rebuilding. Usually not recommended."
            ),
        )

    def handle(self, *args, **options):
        result = rebuild_axis_signature_project_index(
            batch_size=options["batch_size"],
            dry_run=options["dry_run"],
            active_projects_only=not options["include_inactive_projects"],
            clear_existing=not options["no_clear"],
        )

        if result["success"]:
            self.stdout.write(
                self.style.SUCCESS(
                    "AxisSignatureProjectIndex rebuild finished."
                )
            )
        else:
            self.stdout.write(
                self.style.ERROR(
                    "AxisSignatureProjectIndex rebuild failed."
                )
            )

        self.stdout.write(f"dry_run: {result['dry_run']}")
        self.stdout.write(
            f"source_occurrence_count: {result['source_occurrence_count']}"
        )
        self.stdout.write(
            f"deleted_index_count: {result['deleted_index_count']}"
        )
        self.stdout.write(
            f"created_index_count: {result['created_index_count']}"
        )

        if not result["dry_run"]:
            self.stdout.write(
                f"final_index_count: {result['final_index_count']}"
            )

        self.stdout.write(f"batch_size: {result['batch_size']}")
        self.stdout.write(
            f"active_projects_only: {result['active_projects_only']}"
        )
        self.stdout.write(f"clear_existing: {result['clear_existing']}")
