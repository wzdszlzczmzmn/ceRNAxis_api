from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from django.conf import settings
from django.core.management.base import (
    BaseCommand,
    CommandError,
)

from database.utils.dataset_annotation_utils.scst_path_utils import (
    get_scst_dataset_group_by_fields,
    get_scst_group_by_meta_column,
)


META_FILE_SUFFIX = "_meta.csv"


@dataclass(frozen=True, slots=True)
class GroupColumnStats:
    group_by: str
    meta_column: str

    row_count: int
    empty_count: int

    first_empty_rows: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class DatasetMetaStats:
    source: str
    dataset_name: str
    file_path: Path

    row_count: int

    group_stats: tuple[GroupColumnStats, ...]

    missing_columns: tuple[str, ...]

    @property
    def empty_value_count(self) -> int:
        return sum(
            item.empty_count
            for item in self.group_stats
        )

    @property
    def has_empty_values(self) -> bool:
        return self.empty_value_count > 0

    @property
    def has_missing_columns(self) -> bool:
        return bool(self.missing_columns)


def normalize_meta_value(value) -> str:
    """
    Normalize one metadata value for empty-value detection.

    This command treats the following as empty:
        None
        ""
        whitespace-only text

    CSV textual values such as "NA", "None", or "nan" are not
    automatically treated as empty here.
    """
    if value is None:
        return ""

    return str(value).strip()


def get_dataset_name_from_meta_file(
    file_path: Path,
) -> str:
    file_name = file_path.name

    if not file_name.endswith(META_FILE_SUFFIX):
        raise ValueError(
            f"Invalid SC/ST metadata filename: {file_name}"
        )

    dataset_name = file_name[
        :-len(META_FILE_SUFFIX)
    ].strip()

    if not dataset_name:
        raise ValueError(
            f"Unable to resolve dataset name from: {file_name}"
        )

    return dataset_name


def iter_meta_files(
    root_dir: Path,
) -> list[Path]:
    """
    Return all direct/recursive *_meta.csv files.

    Strictly read-only.
    """
    root_dir = (
        Path(root_dir)
        .expanduser()
        .resolve()
    )

    if not root_dir.exists():
        raise CommandError(
            f"SC/ST metadata root does not exist: {root_dir}"
        )

    if not root_dir.is_dir():
        raise CommandError(
            f"SC/ST metadata root is not a directory: {root_dir}"
        )

    return sorted(
        (
            path.resolve()
            for path in root_dir.rglob(
                f"*{META_FILE_SUFFIX}"
            )
            if path.is_file()
        ),
        key=lambda path: str(path).casefold(),
    )


def inspect_dataset_meta_file(
    *,
    source: str,
    file_path: Path,
    max_example_rows: int = 10,
) -> DatasetMetaStats:
    """
    Inspect all configured group columns for one dataset metadata CSV.

    No filesystem or database mutation is performed.
    """

    dataset_name = (
        get_dataset_name_from_meta_file(
            file_path
        )
    )

    group_by_fields = (
        get_scst_dataset_group_by_fields(
            dataset_name
        )
    )

    # Dataset has no configured annotation grouping.
    if not group_by_fields:
        return DatasetMetaStats(
            source=source,
            dataset_name=dataset_name,
            file_path=file_path,
            row_count=0,
            group_stats=(),
            missing_columns=(),
        )

    group_columns = {
        group_by: get_scst_group_by_meta_column(
            group_by
        )
        for group_by in group_by_fields
    }

    try:
        with file_path.open(
            "r",
            encoding="utf-8-sig",
            newline="",
        ) as file_obj:
            reader = csv.DictReader(
                file_obj
            )

            if not reader.fieldnames:
                raise CommandError(
                    "SC/ST metadata file is empty or "
                    "does not contain a header: "
                    f"{file_path}"
                )

            fieldnames = [
                str(column).strip()
                for column in reader.fieldnames
            ]

            reader.fieldnames = fieldnames

            missing_columns = tuple(
                meta_column
                for meta_column
                in group_columns.values()
                if meta_column not in fieldnames
            )

            empty_count_by_group = {
                group_by: 0
                for group_by in group_by_fields
            }

            first_empty_rows_by_group = {
                group_by: []
                for group_by in group_by_fields
            }

            row_count = 0

            for row_number, row in enumerate(
                reader,
                start=2,
            ):
                row_count += 1

                for (
                    group_by,
                    meta_column,
                ) in group_columns.items():
                    # Missing physical column is reported separately.
                    if meta_column not in fieldnames:
                        continue

                    value = normalize_meta_value(
                        row.get(meta_column)
                    )

                    if value:
                        continue

                    empty_count_by_group[
                        group_by
                    ] += 1

                    examples = (
                        first_empty_rows_by_group[
                            group_by
                        ]
                    )

                    if (
                        len(examples)
                        < max_example_rows
                    ):
                        examples.append(
                            row_number
                        )

    except UnicodeDecodeError as exc:
        raise CommandError(
            "SC/ST metadata file must be UTF-8 encoded: "
            f"{file_path}"
        ) from exc

    except csv.Error as exc:
        raise CommandError(
            "Failed to parse SC/ST metadata CSV: "
            f"{file_path}: {exc}"
        ) from exc

    group_stats = tuple(
        GroupColumnStats(
            group_by=group_by,
            meta_column=meta_column,
            row_count=row_count,
            empty_count=(
                empty_count_by_group[
                    group_by
                ]
            ),
            first_empty_rows=tuple(
                first_empty_rows_by_group[
                    group_by
                ]
            ),
        )
        for (
            group_by,
            meta_column,
        ) in group_columns.items()
        if meta_column in fieldnames
    )

    return DatasetMetaStats(
        source=source,
        dataset_name=dataset_name,
        file_path=file_path,
        row_count=row_count,
        group_stats=group_stats,
        missing_columns=missing_columns,
    )


class Command(BaseCommand):
    help = (
        "Read-only validation of configured SC/ST metadata "
        "group columns for empty values."
    )

    def add_arguments(
        self,
        parser,
    ) -> None:
        parser.add_argument(
            "--show-valid",
            action="store_true",
            help=(
                "Also print datasets whose configured group "
                "columns contain no empty values."
            ),
        )

        parser.add_argument(
            "--max-example-rows",
            type=int,
            default=10,
            help=(
                "Maximum number of CSV row numbers displayed "
                "for each group column containing empty values. "
                "Default: 10."
            ),
        )

    def handle(
        self,
        *args,
        **options,
    ) -> None:
        show_valid = options["show_valid"]
        max_example_rows = (
            options["max_example_rows"]
        )

        if max_example_rows <= 0:
            raise CommandError(
                "--max-example-rows must be "
                "greater than zero."
            )

        sources = [
            (
                "SC",
                getattr(
                    settings,
                    "TISCH2_DATASET_BASE_DIR",
                    None,
                ),
            ),
            (
                "ST",
                getattr(
                    settings,
                    "SCTML_DATASET_BASE_DIR",
                    None,
                ),
            ),
        ]

        total_dataset_count = 0
        configured_dataset_count = 0
        invalid_dataset_count = 0
        empty_value_dataset_count = 0
        missing_column_dataset_count = 0

        total_empty_value_count = 0

        self.stdout.write(
            self.style.MIGRATE_HEADING(
                "===== SC/ST Metadata Group Empty Value Check ====="
            )
        )

        for source, root_value in sources:
            if not root_value:
                raise CommandError(
                    f"{source} dataset base directory "
                    "is not configured."
                )

            root_dir = (
                Path(root_value)
                .expanduser()
                .resolve()
            )

            meta_files = iter_meta_files(
                root_dir
            )

            self.stdout.write("")
            self.stdout.write(
                self.style.MIGRATE_HEADING(
                    f"===== {source} ====="
                )
            )

            self.stdout.write(
                f"root:       {root_dir}"
            )
            self.stdout.write(
                f"meta files: {len(meta_files)}"
            )

            for file_path in meta_files:
                total_dataset_count += 1

                try:
                    stats = (
                        inspect_dataset_meta_file(
                            source=source,
                            file_path=file_path,
                            max_example_rows=(
                                max_example_rows
                            ),
                        )
                    )
                except Exception as exc:
                    invalid_dataset_count += 1

                    self.stdout.write("")
                    self.stdout.write(
                        self.style.ERROR(
                            f"[ERROR] "
                            f"{file_path.name}: {exc}"
                        )
                    )
                    continue

                if not stats.group_stats:
                    # No configured group_by for this dataset.
                    continue

                configured_dataset_count += 1

                if stats.has_missing_columns:
                    missing_column_dataset_count += 1

                    self.stdout.write("")
                    self.stdout.write(
                        self.style.ERROR(
                            f"[MISSING COLUMN] "
                            f"{stats.dataset_name}"
                        )
                    )

                    self.stdout.write(
                        f"  file: {stats.file_path}"
                    )
                    self.stdout.write(
                        "  missing columns: "
                        + ", ".join(
                            stats.missing_columns
                        )
                    )

                if stats.has_empty_values:
                    empty_value_dataset_count += 1
                    total_empty_value_count += (
                        stats.empty_value_count
                    )

                    self.stdout.write("")
                    self.stdout.write(
                        self.style.WARNING(
                            f"[EMPTY VALUES] "
                            f"{stats.dataset_name}"
                        )
                    )

                    self.stdout.write(
                        f"  file: {stats.file_path}"
                    )
                    self.stdout.write(
                        f"  rows: {stats.row_count}"
                    )

                    for group_stat in (
                        stats.group_stats
                    ):
                        if (
                            group_stat.empty_count
                            == 0
                        ):
                            continue

                        self.stdout.write(
                            "  - "
                            f"group_by="
                            f"{group_stat.group_by!r}, "
                            f"column="
                            f"{group_stat.meta_column!r}, "
                            f"empty="
                            f"{group_stat.empty_count}"
                        )

                        if (
                            group_stat.first_empty_rows
                        ):
                            row_text = ", ".join(
                                str(row_number)
                                for row_number
                                in (
                                    group_stat
                                    .first_empty_rows
                                )
                            )

                            self.stdout.write(
                                "    example CSV rows: "
                                f"{row_text}"
                            )

                elif (
                    show_valid
                    and not stats.has_missing_columns
                ):
                    self.stdout.write(
                        self.style.SUCCESS(
                            f"[OK] "
                            f"{stats.dataset_name} "
                            f"({stats.row_count} rows)"
                        )
                    )

        self.stdout.write("")
        self.stdout.write(
            self.style.MIGRATE_HEADING(
                "===== Summary ====="
            )
        )

        self.stdout.write(
            "metadata files scanned:       "
            f"{total_dataset_count}"
        )
        self.stdout.write(
            "datasets with group config:   "
            f"{configured_dataset_count}"
        )
        self.stdout.write(
            "datasets with empty values:   "
            f"{empty_value_dataset_count}"
        )
        self.stdout.write(
            "total empty group values:     "
            f"{total_empty_value_count}"
        )
        self.stdout.write(
            "datasets missing group column:"
            f" {missing_column_dataset_count}"
        )
        self.stdout.write(
            "metadata read/parse errors:   "
            f"{invalid_dataset_count}"
        )

        if (
            empty_value_dataset_count == 0
            and missing_column_dataset_count == 0
            and invalid_dataset_count == 0
        ):
            self.stdout.write("")
            self.stdout.write(
                self.style.SUCCESS(
                    "No SC/ST configured group-column "
                    "empty values were found."
                )
            )
