from django.core.management.base import BaseCommand, CommandError, CommandParser
from django.core.management.color import no_style
from django.db import DEFAULT_DB_ALIAS, DatabaseError, connections

from database.models import (
    AxisContextPresence,
    AxisDatasetContext,
    AxisFinalEvidence,
    AxisFinalRecurrentSummary,
    AxisObservation,
    AxisResultArtifact,
    AxisStructureRecurrentSummary,
    CanonicalAxis,
    SpongeEvidence,
)


# Children come first for database backends that clear tables with DELETE.
RECURRENT_AXIS_MODELS = (
    AxisFinalEvidence,
    SpongeEvidence,
    AxisObservation,
    AxisContextPresence,
    AxisResultArtifact,
    AxisStructureRecurrentSummary,
    AxisFinalRecurrentSummary,
    AxisDatasetContext,
    CanonicalAxis,
)


class Command(BaseCommand):
    help = (
        "Clear all current and historical data in the nine recurrent-axis "
        "tables. By default, report row counts. Pass --confirm to delete."
    )

    def add_arguments(self, parser: CommandParser) -> None:
        mode = parser.add_mutually_exclusive_group()
        mode.add_argument(
            "--dry-run",
            action="store_true",
            help="Report how many rows would be deleted (default).",
        )
        mode.add_argument(
            "--confirm",
            action="store_true",
            help="Delete all rows in the nine recurrent-axis tables.",
        )
        parser.add_argument(
            "--database",
            choices=tuple(connections),
            default=DEFAULT_DB_ALIAS,
            help="Django database alias to target. Defaults to 'default'.",
        )

    def handle(self, *args, **options) -> None:
        database = options["database"]
        if options["dry_run"] or not options["confirm"]:
            self.stdout.write(f"Dry run on database '{database}':")
            for model in RECURRENT_AXIS_MODELS:
                count = model.objects.using(database).count()
                self.stdout.write(f"  {model._meta.db_table}: {count}")
            return

        connection = connections[database]
        if not connection.features.can_rollback_ddl:
            raise CommandError(
                f"Database '{database}' does not support transactional table clearing."
            )

        tables = [model._meta.db_table for model in RECURRENT_AXIS_MODELS]
        try:
            statements = connection.ops.sql_flush(
                no_style(), tables, reset_sequences=False, allow_cascade=False
            )
            connection.ops.execute_sql_flush(statements)
        except DatabaseError as exc:
            raise CommandError(
                f"Could not clear recurrent-axis data on database '{database}': {exc}"
            ) from exc

        self.stdout.write(
            self.style.SUCCESS(
                f"Cleared {len(tables)} recurrent-axis tables on database '{database}'."
            )
        )
