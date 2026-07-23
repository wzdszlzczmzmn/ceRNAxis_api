from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq


def parse_pandas_metadata(metadata: dict[bytes, bytes] | None) -> dict[str, Any] | None:
    """解析 Parquet 中由 pandas 写入的 metadata。"""
    if not metadata:
        return None

    raw_metadata = metadata.get(b"pandas")
    if not raw_metadata:
        return None

    try:
        return json.loads(raw_metadata.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        return None


def inspect_parquet(
    file_path: Path,
    preview_rows: int = 5,
    preview_columns: int = 10,
) -> None:
    if not file_path.exists():
        raise FileNotFoundError(f"File does not exist: {file_path}")

    if not file_path.is_file():
        raise ValueError(f"Path is not a file: {file_path}")

    if file_path.suffix.lower() != ".parquet":
        raise ValueError(
            f"Expected a .parquet file, got: {file_path.suffix or '<no suffix>'}"
        )

    parquet_file = pq.ParquetFile(file_path)
    schema = parquet_file.schema_arrow
    metadata = parquet_file.metadata
    column_names = schema.names

    print("=" * 80)
    print("FILE")
    print("=" * 80)
    print(f"Path:             {file_path.resolve()}")
    print(f"Size:             {file_path.stat().st_size:,} bytes")
    print(f"Rows:             {metadata.num_rows:,}")
    print(f"Columns:          {metadata.num_columns:,}")
    print(f"Row groups:       {metadata.num_row_groups:,}")
    print(f"Created by:       {metadata.created_by}")
    print()

    print("=" * 80)
    print("PHYSICAL COLUMN ORDER")
    print("=" * 80)

    for index, field in enumerate(schema):
        marker = "  <-- first physical column" if index == 0 else ""
        print(
            f"[{index:>4}] "
            f"name={field.name!r}, "
            f"type={field.type}, "
            f"nullable={field.nullable}"
            f"{marker}"
        )

    print()

    pandas_metadata = parse_pandas_metadata(schema.metadata)

    print("=" * 80)
    print("PANDAS METADATA")
    print("=" * 80)

    if pandas_metadata is None:
        print("No readable pandas metadata was found.")
        index_columns = []
    else:
        index_columns = pandas_metadata.get("index_columns", [])
        print(f"index_columns: {index_columns}")

        pandas_columns = pandas_metadata.get("columns", [])
        print("Recorded pandas columns:")

        for column in pandas_columns:
            print(
                f"  name={column.get('name')!r}, "
                f"field_name={column.get('field_name')!r}, "
                f"pandas_type={column.get('pandas_type')!r}, "
                f"numpy_type={column.get('numpy_type')!r}"
            )

    print()

    print("=" * 80)
    print("PANDAS PREVIEW")
    print("=" * 80)

    selected_columns = column_names[:preview_columns]

    table = parquet_file.read(columns=selected_columns)
    dataframe = table.to_pandas()

    print(f"DataFrame shape for preview: {dataframe.shape}")
    print(f"DataFrame index name:        {dataframe.index.name!r}")
    print(f"DataFrame index type:        {type(dataframe.index).__name__}")
    print(f"Previewed columns:           {list(dataframe.columns)}")
    print()
    print(dataframe.head(preview_rows).to_string())

    print()

    print("=" * 80)
    print("IDENTIFIER DIAGNOSIS")
    print("=" * 80)

    identifier_candidates = ["cell_id", "spot_id", "sample_id"]

    normal_identifier_columns = [
        name for name in identifier_candidates
        if name in column_names
    ]

    if normal_identifier_columns:
        print(
            "Identifier column(s) found in physical schema: "
            f"{normal_identifier_columns}"
        )
    else:
        print(
            "No cell_id, spot_id, or sample_id column was found "
            "in the physical schema."
        )

    if dataframe.index.name in identifier_candidates:
        print(
            f"Identifier appears to be stored as the pandas index: "
            f"{dataframe.index.name!r}"
        )

    if index_columns:
        print(
            "Parquet pandas metadata declares index column(s): "
            f"{index_columns}"
        )

    first_column = column_names[0] if column_names else None
    print(f"First physical column: {first_column!r}")

    if first_column not in identifier_candidates:
        print(
            "WARNING: the first physical column is not cell_id, spot_id, "
            "or sample_id."
        )

        if dataframe.index.name in identifier_candidates or index_columns:
            print(
                "Likely cause: the identifier was written as a pandas index, "
                "not as a regular Parquet column."
            )
        else:
            print(
                "Possible causes:\n"
                "1. The matrix orientation is genes x cells/spots.\n"
                "2. The identifier column was removed before writing.\n"
                "3. The identifier uses another name.\n"
                "4. The first gene column is genuinely the first stored column."
            )

    print()

    print("=" * 80)
    print("ORIENTATION HEURISTIC")
    print("=" * 80)

    gene_like_first_columns = [
        name for name in column_names[:10]
        if isinstance(name, str)
        and (
            name.startswith("ENSG")
            or name.startswith("ENSMUSG")
            or name.upper() == name
        )
    ]

    if gene_like_first_columns:
        print(
            "The first columns look gene-like: "
            f"{gene_like_first_columns[:10]}"
        )
        print(
            "This is consistent with rows representing cells/spots and "
            "columns representing genes, but the row identifier may be stored "
            "as the DataFrame index."
        )
    else:
        print(
            "No clear conclusion could be made from the first column names."
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Inspect the structure and metadata of a Parquet expression file."
    )
    parser.add_argument(
        "parquet_file",
        type=Path,
        help="Path to the Parquet file.",
    )
    parser.add_argument(
        "--rows",
        type=int,
        default=5,
        help="Number of preview rows. Default: 5.",
    )
    parser.add_argument(
        "--columns",
        type=int,
        default=10,
        help="Number of preview columns. Default: 10.",
    )

    args = parser.parse_args()

    try:
        inspect_parquet(
            file_path=args.parquet_file,
            preview_rows=args.rows,
            preview_columns=args.columns,
        )
    except Exception as exc:
        raise SystemExit(f"Inspection failed: {exc}") from exc


if __name__ == "__main__":
    main()
