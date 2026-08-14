from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from django.conf import settings
from django.core.management.base import (
    BaseCommand,
    CommandError,
)


AXIS_FINAL_FILE_PATTERN = "*_ceRNA_axis_final.csv"


@dataclass(frozen=True, slots=True)
class AxisFinalFileStats:
    file_count: int = 0
    empty_file_count: int = 0
    data_row_count: int = 0


def count_csv_data_rows(
    file_path: Path,
) -> tuple[int, bool]:
    """
    Count data rows in one Axis Final CSV.

    Returns:
        (
            data_row_count,
            is_empty_file,
        )

    Rules:
        - zero-byte file:
            0 data rows, no header expected;
        - non-empty file:
            first CSV row is treated as header;
        - blank rows after the header are ignored.

    This function is strictly read-only.
    """

    file_size = file_path.stat().st_size

    if file_size == 0:
        return 0, True

    with file_path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as file_obj:
        reader = csv.reader(file_obj)

        # First row is the header.
        header = next(reader, None)

        # Defensive handling for a non-zero file that nevertheless
        # produces no CSV rows.
        if header is None:
            return 0, True

        data_row_count = 0

        for row in reader:
            # Ignore completely blank CSV rows.
            if not row:
                continue

            if not any(
                str(value).strip()
                for value in row
            ):
                continue

            data_row_count += 1

    return data_row_count, False


def collect_axis_final_stats(
    root_dir: Path,
) -> AxisFinalFileStats:
    """
    Recursively scan one annotation root and count Axis Final files
    and their data rows.

    This function never modifies the filesystem.
    """

    root_dir = Path(root_dir).expanduser().resolve()

    if not root_dir.exists():
        raise CommandError(
            f"Annotation root does not exist: {root_dir}"
        )

    if not root_dir.is_dir():
        raise CommandError(
            f"Annotation root is not a directory: {root_dir}"
        )

    file_count = 0
    empty_file_count = 0
    data_row_count = 0

    for file_path in sorted(
        root_dir.rglob(
            AXIS_FINAL_FILE_PATTERN
        ),
        key=lambda path: str(path).casefold(),
    ):
        if not file_path.is_file():
            continue

        file_count += 1

        try:
            (
                row_count,
                is_empty,
            ) = count_csv_data_rows(
                file_path
            )
        except (
            OSError,
            UnicodeDecodeError,
            csv.Error,
        ) as exc:
            raise CommandError(
                "Failed to read Axis Final file: "
                f"{file_path}. "
                f"{type(exc).__name__}: {exc}"
            ) from exc

        data_row_count += row_count

        if is_empty:
            empty_file_count += 1

    return AxisFinalFileStats(
        file_count=file_count,
        empty_file_count=empty_file_count,
        data_row_count=data_row_count,
    )


class Command(BaseCommand):
    help = (
        "Read-only statistics for TCGA and TIMEDB Dataset "
        "Annotation Axis Final files."
    )

    def handle(self, *args, **options) -> None:
        sources = [
            (
                "TCGA",
                getattr(
                    settings,
                    "TCGA_DATASET_ANNOTATIONS_DIR",
                    None,
                ),
            ),
            (
                "TIMEDB",
                getattr(
                    settings,
                    "TIMEDB_DATASET_ANNOTATIONS_DIR",
                    None,
                ),
            ),
        ]

        total_file_count = 0
        total_empty_file_count = 0
        total_data_row_count = 0

        self.stdout.write(
            self.style.MIGRATE_HEADING(
                "===== Axis Final Statistics ====="
            )
        )

        for source_name, root_value in sources:
            if not root_value:
                raise CommandError(
                    f"{source_name} annotation root "
                    "is not configured."
                )

            root_dir = Path(
                root_value
            ).expanduser().resolve()

            stats = collect_axis_final_stats(
                root_dir
            )

            self.stdout.write("")
            self.stdout.write(
                self.style.MIGRATE_LABEL(
                    f"===== {source_name} ====="
                )
            )

            self.stdout.write(
                f"root:              {root_dir}"
            )
            self.stdout.write(
                "axis final files:  "
                f"{stats.file_count}"
            )
            self.stdout.write(
                "empty files:       "
                f"{stats.empty_file_count}"
            )
            self.stdout.write(
                "data rows:         "
                f"{stats.data_row_count}"
            )

            total_file_count += (
                stats.file_count
            )
            total_empty_file_count += (
                stats.empty_file_count
            )
            total_data_row_count += (
                stats.data_row_count
            )

        self.stdout.write("")
        self.stdout.write(
            self.style.MIGRATE_HEADING(
                "===== Total ====="
            )
        )

        self.stdout.write(
            f"axis final files:  {total_file_count}"
        )
        self.stdout.write(
            f"empty files:       {total_empty_file_count}"
        )
        self.stdout.write(
            f"data rows:         {total_data_row_count}"
        )
