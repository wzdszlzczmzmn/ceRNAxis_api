from django.core.management.base import BaseCommand
from django.db import transaction

from database.models import (
    RNAInteraction,
    InteractionDatabase,
    FilterField,
    FilterOption,
)


class Command(BaseCommand):
    help = "Initialize FilterField and FilterOption from database data"

    TABLE_NAME = "rna_interaction"

    FIELD_CONFIGS = [
        {
            "field_name": "species",
            "field_label": "Species",
            "sort_order": 1,
        },
        {
            "field_name": "database",
            "field_label": "Database",
            "sort_order": 2,
        },
        {
            "field_name": "type",
            "field_label": "Type",
            "sort_order": 3,
        },
    ]

    def handle(self, *args, **options):
        with transaction.atomic():
            field_map = self.create_filter_fields()
            self.create_species_options(field_map["species"])
            self.create_database_options(field_map["database"])
            self.create_type_options(field_map["type"])

        self.stdout.write(self.style.SUCCESS("Filter fields and options initialized."))

    def create_filter_fields(self):
        field_map = {}

        for config in self.FIELD_CONFIGS:
            field, created = FilterField.objects.update_or_create(
                table_name=self.TABLE_NAME,
                field_name=config["field_name"],
                defaults={
                    "field_label": config["field_label"],
                    "sort_order": config["sort_order"],
                    "is_active": True,
                },
            )

            field_map[config["field_name"]] = field

            self.stdout.write(
                f"FilterField {'created' if created else 'updated'}: "
                f"{field.table_name}.{field.field_name}"
            )

        return field_map

    def create_species_options(self, field):
        values = (
            RNAInteraction.objects
            .exclude(species__isnull=True)
            .exclude(species="")
            .values_list("species", flat=True)
            .distinct()
            .order_by("species")
        )

        self.bulk_upsert_options(field, values)

    def create_database_options(self, field):
        values = (
            InteractionDatabase.objects
            .exclude(name__isnull=True)
            .exclude(name="")
            .values_list("name", flat=True)
            .distinct()
            .order_by("name")
        )

        self.bulk_upsert_options(field, values)

    def create_type_options(self, field):
        values = (
            RNAInteraction.objects
            .exclude(interaction_type__isnull=True)
            .exclude(interaction_type="")
            .values_list("interaction_type", flat=True)
            .distinct()
            .order_by("interaction_type")
        )

        self.bulk_upsert_options(field, values)

    def bulk_upsert_options(self, field, values):
        values = list(values)

        existing_values = set(
            FilterOption.objects
            .filter(field=field)
            .values_list("value", flat=True)
        )

        objs = []

        for index, value in enumerate(values, start=1):
            value = value.strip()

            if not value:
                continue

            if value in existing_values:
                FilterOption.objects.filter(
                    field=field,
                    value=value,
                ).update(
                    sort_order=index,
                    is_active=True,
                )
                continue

            objs.append(
                FilterOption(
                    field=field,
                    value=value,
                    sort_order=index,
                    is_active=True,
                )
            )

        FilterOption.objects.bulk_create(
            objs,
            batch_size=1000,
            ignore_conflicts=True,
        )

        self.stdout.write(
            f"FilterOption initialized for {field.field_name}: "
            f"{len(values)} values, {len(objs)} newly created"
        )
