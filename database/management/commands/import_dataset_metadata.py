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
    ]

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
            col for col in self.REQUIRED_COLUMNS if col not in df.columns
        ]

        if missing_columns:
            raise CommandError(
                f"Missing required columns: {', '.join(missing_columns)}"
            )

        df = df[self.REQUIRED_COLUMNS].copy()

        # 清理空行：dataset 为空的行直接忽略
        df = df.dropna(subset=["dataset"])

        # 将 NaN 转为空字符串，避免 CharField 写入 nan
        text_columns = [
            "dataset",
            "c_programme",
            "c_obs_type",
            "c_reference",
            "c_cancer_type",
            "c_cancer_type_full_name",
            "c_gene_bio_type",
            "c_workflow",
        ]

        for col in text_columns:
            df[col] = df[col].fillna("").astype(str).str.strip()

        df["n_sample_nums"] = df["n_sample_nums"].fillna(0).astype(int)

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
                ],
                unique_fields=["dataset"],
            )

        self.stdout.write(
            self.style.SUCCESS(
                f"Successfully imported {len(objects)} dataset metadata records."
            )
        )
