import argparse
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


def check_parquet_zero_ratio(
    parquet_file: str,
    nrows: int = 10000,
    include_first_column: bool = False,
):
    parquet_path = Path(parquet_file)

    if not parquet_path.exists() or not parquet_path.is_file():
        raise FileNotFoundError(f"Parquet file not found: {parquet_path}")

    parquet = pq.ParquetFile(parquet_path)

    columns = parquet.schema.names

    if not columns:
        raise ValueError("Parquet file has no columns.")

    first_column = columns[0]

    print("=" * 80)
    print("Parquet file summary")
    print("=" * 80)
    print(f"File: {parquet_path}")
    print(f"Total rows: {parquet.metadata.num_rows}")
    print(f"Total columns: {len(columns)}")
    print(f"First column: {first_column}")
    print()

    read_nrows = min(nrows, parquet.metadata.num_rows)

    if read_nrows == 0:
        print("File has zero rows. Skip zero-ratio calculation.")
        return

    batches = []
    remaining = read_nrows

    for batch in parquet.iter_batches(batch_size=read_nrows):
        take_rows = min(batch.num_rows, remaining)
        batches.append(batch.slice(0, take_rows))
        remaining -= take_rows

        if remaining <= 0:
            break

    df = pd.concat(
        [batch.to_pandas() for batch in batches],
        ignore_index=True,
    )

    if include_first_column:
        value_df = df
        checked_columns = columns
    else:
        value_df = df.iloc[:, 1:]
        checked_columns = columns[1:]

    if value_df.empty:
        print("No value columns available for zero-ratio calculation.")
        return

    numeric_df = value_df.apply(pd.to_numeric, errors="coerce")

    total_cells = numeric_df.size
    valid_cells = numeric_df.notna().sum().sum()
    zero_cells = (numeric_df == 0).sum().sum()

    zero_ratio_total = zero_cells / total_cells if total_cells else 0
    zero_ratio_valid = zero_cells / valid_cells if valid_cells else 0

    print("=" * 80)
    print(f"Zero ratio in first {read_nrows} rows")
    print("=" * 80)
    print(f"Include first column: {include_first_column}")
    print(f"Checked columns: {len(checked_columns)}")
    print(f"Total checked cells: {total_cells}")
    print(f"Valid numeric cells: {valid_cells}")
    print(f"Zero cells: {zero_cells}")
    print(f"Zero ratio among all checked cells: {zero_ratio_total:.6f}")
    print(f"Zero ratio among valid numeric cells: {zero_ratio_valid:.6f}")


def main():
    parser = argparse.ArgumentParser(
        description="Check first column and zero ratio in the first N rows of a parquet file."
    )

    parser.add_argument(
        "parquet_file",
        help="Path to the parquet file.",
    )

    parser.add_argument(
        "--nrows",
        type=int,
        default=10000,
        help="Number of rows to check. Default: 10000.",
    )

    parser.add_argument(
        "--include-first-column",
        action="store_true",
        help="Include the first column in zero-ratio calculation.",
    )

    args = parser.parse_args()

    check_parquet_zero_ratio(
        parquet_file=args.parquet_file,
        nrows=args.nrows,
        include_first_column=args.include_first_column,
    )


if __name__ == "__main__":
    main()
