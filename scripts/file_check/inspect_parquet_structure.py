import argparse
import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


def print_section(title: str):
    print()
    print("=" * 80)
    print(title)
    print("=" * 80)


def inspect_parquet_structure(
    parquet_file: str,
    nrows: int = 5,
    ncols: int = 5,
):
    parquet_path = Path(parquet_file)

    if not parquet_path.exists() or not parquet_path.is_file():
        raise FileNotFoundError(f"Parquet file not found: {parquet_path}")

    parquet = pq.ParquetFile(parquet_path)
    schema = parquet.schema_arrow
    columns = schema.names

    print_section("Basic parquet info")
    print(f"File: {parquet_path}")
    print(f"Total rows: {parquet.metadata.num_rows}")
    print(f"Total columns in schema: {len(columns)}")
    print(f"Row groups: {parquet.num_row_groups}")

    print_section("Schema columns")
    print(f"First {min(ncols, len(columns))} columns:")
    for i, col in enumerate(columns[:ncols]):
        print(f"  [{i}] {col}")

    print()
    print(f"Last {min(ncols, len(columns))} columns:")
    for i, col in enumerate(columns[-ncols:], start=max(len(columns) - ncols, 0)):
        print(f"  [{i}] {col}")

    print_section("Arrow schema")
    print(schema)

    print_section("Parquet key-value metadata")
    metadata = schema.metadata or {}

    if not metadata:
        print("No key-value metadata found.")
    else:
        for key, value in metadata.items():
            key_text = key.decode("utf-8", errors="replace")
            value_text = value.decode("utf-8", errors="replace")

            if key_text == "pandas":
                print("pandas metadata found.")
                try:
                    pandas_meta = json.loads(value_text)
                    print(json.dumps(pandas_meta, indent=2, ensure_ascii=False)[:5000])
                    if len(json.dumps(pandas_meta)) > 5000:
                        print("... <truncated>")
                except json.JSONDecodeError:
                    print(value_text[:5000])
            else:
                print(f"{key_text}: {value_text[:1000]}")

    if not columns:
        print_section("Read sample data")
        print("No columns available.")
        return

    selected_columns = columns[: min(ncols, len(columns))]

    print_section("Read sample with pandas")
    print(f"Reading columns: {selected_columns}")
    print(f"Reading rows: {nrows}")

    df = pd.read_parquet(
        parquet_path,
        columns=selected_columns,
        engine="pyarrow",
    ).head(nrows)

    print()
    print("DataFrame shape:")
    print(df.shape)

    print()
    print("Index info:")
    print(f"  index type: {type(df.index)}")
    print(f"  index name: {df.index.name}")
    print(f"  index names: {getattr(df.index, 'names', None)}")
    print(f"  first index values: {list(df.index[:nrows])}")

    print()
    print("DataFrame columns:")
    for i, col in enumerate(df.columns.tolist()):
        print(f"  [{i}] {col}")

    print()
    print("DataFrame head:")
    print(df)

    print_section("After reset_index")
    reset_df = df.reset_index()

    print("reset_index columns:")
    for i, col in enumerate(reset_df.columns.tolist()):
        print(f"  [{i}] {col}")

    print()
    print("reset_index head:")
    print(reset_df)


def main():
    parser = argparse.ArgumentParser(
        description="Inspect parquet columns, pandas index, and reset_index structure."
    )

    parser.add_argument(
        "parquet_file",
        help="Path to parquet file.",
    )

    parser.add_argument(
        "--nrows",
        type=int,
        default=5,
        help="Number of rows to preview. Default: 5.",
    )

    parser.add_argument(
        "--ncols",
        type=int,
        default=5,
        help="Number of schema columns to preview and read. Default: 5.",
    )

    args = parser.parse_args()

    inspect_parquet_structure(
        parquet_file=args.parquet_file,
        nrows=args.nrows,
        ncols=args.ncols,
    )


if __name__ == "__main__":
    main()
