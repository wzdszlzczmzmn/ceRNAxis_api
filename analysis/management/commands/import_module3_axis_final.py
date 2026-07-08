from django.core.management.base import BaseCommand

from analysis.services.axis_final import import_module3_axis_final_projects


class Command(BaseCommand):
    help = "Import Module3 TIMEDB axis_final results into database."

    def add_arguments(self, parser):
        parser.add_argument(
            "--root",
            required=True,
            help="Module3 result root directory.",
        )

        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Read and validate files without writing database.",
        )

        parser.add_argument(
            "--show-results",
            action="store_true",
            help="Print per-project import results.",
        )

    def handle(self, *args, **options):
        result = import_module3_axis_final_projects(
            module3_root_dir=options["root"],
            dry_run=options["dry_run"],
        )

        if result["success"]:
            self.stdout.write(
                self.style.SUCCESS(
                    "Module3 axis_final import finished."
                )
            )
        else:
            self.stdout.write(
                self.style.ERROR(
                    "Module3 axis_final import finished with errors."
                )
            )

        self.stdout.write(f"root: {result['module3_root_dir']}")
        self.stdout.write(f"source: {result['source']}")
        self.stdout.write(f"module: {result['module']}")
        self.stdout.write(f"annotation_dir_count: {result['annotation_dir_count']}")
        self.stdout.write(f"project_count: {result['project_count']}")
        self.stdout.write(f"imported_count: {result['imported_count']}")
        self.stdout.write(f"dry_run_count: {result['dry_run_count']}")
        self.stdout.write(f"skipped_count: {result['skipped_count']}")
        self.stdout.write(f"failed_count: {result['failed_count']}")

        if options["show_results"]:
            for item in result["results"]:
                if item.get("success"):
                    self.stdout.write(self.style.SUCCESS(str(item)))
                elif item.get("skipped"):
                    self.stdout.write(self.style.WARNING(str(item)))
                else:
                    self.stdout.write(self.style.ERROR(str(item)))

        if not result["success"]:
            raise SystemExit(1)
