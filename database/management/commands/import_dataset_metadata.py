import pandas as pd

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from database.models import DatasetMetadata


class Command(BaseCommand):
    help = "Import Dataset Metadata v2 sheet from xlsx file."

    REQUIRED_COLUMNS = [
        "dataset",
        "c_programme",
        "c_obs_type",
        "c_reference",
        "c_cancer_type",
        "c_cancer_type_full_name",
        "c_gene_bio_type",
        "c_workflow",
        "n_sample_nums",
        "n_cell_nums",
        "n_spot_nums",
    ]

    TEXT_COLUMNS = [
        "dataset",
        "c_programme",
        "c_obs_type",
        "c_reference",
        "c_cancer_type",
        "c_cancer_type_full_name",
        "c_gene_bio_type",
        "c_workflow",
    ]

    COUNT_COLUMNS = [
        "n_sample_nums",
        "n_cell_nums",
        "n_spot_nums",
    ]

    VALID_GENE_BIO_TYPES = {
        "miRNA",
        "mRNA",
        "lncRNA",
        "circRNA",
    }

    VALID_OBS_TYPES = {
        "sample",
        "cell",
        "spot",
    }

    def add_arguments(self, parser):
        parser.add_argument(
            "xlsx_path",
            type=str,
            help="Path to the xlsx file.",
        )
        parser.add_argument(
            "--sheet",
            type=str,
            default="Dataset Metadata v2",
            help="Sheet name to import.",
        )
        parser.add_argument(
            "--clear",
            action="store_true",
            help="Clear existing DatasetMetadata records before import.",
        )

    def handle(self, *args, **options):
        xlsx_path = options["xlsx_path"]
        sheet_name = options["sheet"]
        clear = options["clear"]

        try:
            df = pd.read_excel(xlsx_path, sheet_name=sheet_name)
        except Exception as e:
            raise CommandError(f"Failed to read xlsx file: {e}")

        missing_columns = [
            col for col in self.REQUIRED_COLUMNS
            if col not in df.columns
        ]

        if missing_columns:
            raise CommandError(
                f"Missing required columns: {', '.join(missing_columns)}"
            )

        df = df[self.REQUIRED_COLUMNS].copy()
        raw_count = len(df)

        # dataset 为空的行没有业务意义
        df = df.dropna(subset=["dataset"]).copy()

        for col in self.TEXT_COLUMNS:
            df[col] = df[col].fillna("").astype(str).str.strip()

        # 清洗后 dataset 为空字符串的行继续丢弃
        df = df[df["dataset"] != ""].copy()

        invalid_obs_types = sorted(
            set(df["c_obs_type"]) - self.VALID_OBS_TYPES
        )
        if invalid_obs_types:
            raise CommandError(
                "Invalid c_obs_type values: "
                + ", ".join(invalid_obs_types)
            )

        invalid_gene_types = sorted(
            set(df["c_gene_bio_type"]) - self.VALID_GENE_BIO_TYPES
        )
        if invalid_gene_types:
            raise CommandError(
                "Invalid c_gene_bio_type values: "
                + ", ".join(invalid_gene_types)
            )

        duplicated_datasets = (
            df[df.duplicated(subset=["dataset"], keep=False)]["dataset"]
            .drop_duplicates()
            .tolist()
        )

        if duplicated_datasets:
            raise CommandError(
                "Duplicated dataset values in input file: "
                + ", ".join(duplicated_datasets[:20])
            )

        for col in self.COUNT_COLUMNS:
            df[col] = pd.to_numeric(
                df[col],
                errors="coerce",
            ).fillna(0)

            negative_rows = df[df[col] < 0]
            if not negative_rows.empty:
                raise CommandError(
                    f"{col} contains negative values: "
                    + ", ".join(
                        negative_rows["dataset"].head(10).tolist()
                    )
                )

            df[col] = df[col].astype(int)

        objects = []

        for _, row in df.iterrows():
            objects.append(
                DatasetMetadata(
                    dataset=row["dataset"],
                    programme=row["c_programme"],
                    obs_type=row["c_obs_type"],
                    reference=row["c_reference"],
                    cancer_type=row["c_cancer_type"],
                    cancer_type_full_name=row["c_cancer_type_full_name"],
                    gene_bio_type=row["c_gene_bio_type"],
                    workflow=row["c_workflow"],
                    sample_nums=row["n_sample_nums"],
                    cell_nums=row["n_cell_nums"],
                    spot_nums=row["n_spot_nums"],
                )
            )

        with transaction.atomic():
            if clear:
                deleted_count, _ = DatasetMetadata.objects.all().delete()
                self.stdout.write(
                    self.style.WARNING(
                        f"Deleted existing records: {deleted_count}"
                    )
                )

            DatasetMetadata.objects.bulk_create(
                objects,
                update_conflicts=True,
                update_fields=[
                    "programme",
                    "obs_type",
                    "reference",
                    "cancer_type",
                    "cancer_type_full_name",
                    "gene_bio_type",
                    "workflow",
                    "sample_nums",
                    "cell_nums",
                    "spot_nums",
                ],
                unique_fields=["dataset"],
            )

        self.stdout.write(
            self.style.SUCCESS(
                (
                    f"Successfully imported {len(objects)} dataset metadata "
                    f"records. Raw rows: {raw_count}, "
                    f"ignored rows: {raw_count - len(df)}."
                )
            )
        )
