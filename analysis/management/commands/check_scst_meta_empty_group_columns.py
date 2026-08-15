from __future__ import annotations

import csv
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
    Recursively return all *_meta.csv files.

    Strictly read-only.
    """
    root_dir = (
        Path(root_dir)
        .expanduser()
        .resolve()
    )

    if not root_dir.exists():
        raise CommandError(
            f"Metadata root does not exist: {root_dir}"
        )

    if not root_dir.is_dir():
        raise CommandError(
            f"Metadata root is not a directory: {root_dir}"
        )

    return sorted(
        (
            file_path.resolve()
            for file_path in root_dir.rglob(
                f"*{META_FILE_SUFFIX}"
            )
            if file_path.is_file()
        ),
        key=lambda path: str(path).casefold(),
    )


def is_empty_group_value(value) -> bool:
    """
    Match the current SC/ST grouping semantics:
    empty/whitespace-only values do not form a group.
    """
    if value is None:
        return True

    return not str(value).strip()


def find_all_empty_group_columns(
    *,
    file_path: Path,
) -> tuple[
    int,
    list[dict],
]:
    """
    Return:
        (
            data_row_count,
            [
                {
                    "group_by": ...,
                    "meta_column": ...,
                },
                ...
            ],
        )

    Only configured group-by columns are inspected.
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

    if not group_by_fields:
        return 0, []

    meta_column_by_group_by = {
        group_by: (
            get_scst_group_by_meta_column(
                group_by
            )
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
                    "Metadata file is empty or "
                    "missing a header: "
                    f"{file_path}"
                )

            fieldnames = [
                str(column).strip()
                for column in reader.fieldnames
            ]

            reader.fieldnames = fieldnames

            missing_columns = [
                (
                    group_by,
                    meta_column,
                )
                for (
                    group_by,
                    meta_column,
                ) in (
                    meta_column_by_group_by.items()
                )
                if meta_column not in fieldnames
            ]

            if missing_columns:
                missing_text = ", ".join(
                    (
                        f"{group_by} -> "
                        f"{meta_column}"
                    )
                    for (
                        group_by,
                        meta_column,
                    ) in missing_columns
                )

                raise CommandError(
                    "Metadata file is missing "
                    "configured group column(s): "
                    f"{missing_text}. "
                    f"File: {file_path}"
                )

            has_nonempty_value = {
                group_by: False
                for group_by in group_by_fields
            }

            row_count = 0

            for row in reader:
                row_count += 1

                for group_by in group_by_fields:
                    # Once a non-empty value has been seen,
                    # this column can no longer be "all empty".
                    if has_nonempty_value[group_by]:
                        continue

                    meta_column = (
                        meta_column_by_group_by[
                            group_by
                        ]
                    )

                    value = row.get(
                        meta_column
                    )

                    if not is_empty_group_value(
                        value
                    ):
                        has_nonempty_value[
                            group_by
                        ] = True

    except UnicodeDecodeError as exc:
        raise CommandError(
            "Metadata file must be UTF-8 encoded: "
            f"{file_path}"
        ) from exc

    except csv.Error as exc:
        raise CommandError(
            "Invalid metadata CSV file: "
            f"{file_path}: {exc}"
        ) from exc

    if row_count == 0:
        raise CommandError(
            "Metadata file contains no data rows: "
            f"{file_path}"
        )

    all_empty_columns = [
        {
            "group_by": group_by,
            "meta_column": (
                meta_column_by_group_by[
                    group_by
                ]
            ),
        }
        for group_by in group_by_fields
        if not has_nonempty_value[
            group_by
        ]
    ]

    return (
        row_count,
        all_empty_columns,
    )


class Command(BaseCommand):
    help = (
        "Read-only scan of SC/ST metadata files for "
        "configured group columns whose values are all empty."
    )

    def handle(
        self,
        *args,
        **options,
    ) -> None:
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

        scanned_meta_count = 0
        configured_meta_count = 0
        affected_dataset_count = 0
        all_empty_column_count = 0
        error_count = 0

        self.stdout.write(
            self.style.MIGRATE_HEADING(
                "===== SC/ST All-Empty Group Column Check ====="
            )
        )

        for (
            source,
            root_value,
        ) in sources:
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

            meta_files = (
                iter_meta_files(
                    root_dir
                )
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
                scanned_meta_count += 1

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

                if not group_by_fields:
                    continue

                configured_meta_count += 1

                try:
                    (
                        row_count,
                        all_empty_columns,
                    ) = find_all_empty_group_columns(
                        file_path=file_path
                    )
                except Exception as exc:
                    error_count += 1

                    self.stdout.write("")
                    self.stdout.write(
                        self.style.ERROR(
                            "[ERROR] "
                            f"{dataset_name}: "
                            f"{exc}"
                        )
                    )
                    continue

                if not all_empty_columns:
                    continue

                affected_dataset_count += 1
                all_empty_column_count += (
                    len(all_empty_columns)
                )

                self.stdout.write("")
                self.stdout.write(
                    self.style.WARNING(
                        "[ALL EMPTY GROUP COLUMN] "
                        f"{dataset_name}"
                    )
                )

                self.stdout.write(
                    f"  source: {source}"
                )
                self.stdout.write(
                    f"  file:   {file_path}"
                )
                self.stdout.write(
                    f"  rows:   {row_count}"
                )

                for item in (
                    all_empty_columns
                ):
                    self.stdout.write(
                        "  - "
                        f"group_by="
                        f"{item['group_by']!r}, "
                        f"meta_column="
                        f"{item['meta_column']!r}"
                    )

        self.stdout.write("")
        self.stdout.write(
            self.style.MIGRATE_HEADING(
                "===== Summary ====="
            )
        )

        self.stdout.write(
            "metadata files scanned:       "
            f"{scanned_meta_count}"
        )
        self.stdout.write(
            "datasets with group config:   "
            f"{configured_meta_count}"
        )
        self.stdout.write(
            "datasets with all-empty group:"
            f" {affected_dataset_count}"
        )
        self.stdout.write(
            "all-empty group columns:      "
            f"{all_empty_column_count}"
        )
        self.stdout.write(
            "metadata errors:              "
            f"{error_count}"
        )

        if (
            affected_dataset_count == 0
            and error_count == 0
        ):
            self.stdout.write("")
            self.stdout.write(
                self.style.SUCCESS(
                    "No configured SC/ST group column "
                    "was entirely empty."
                )
            )
