from pathlib import Path
from typing import List

import pandas as pd
import typer


app = typer.Typer(help="Convert RNA CSV files to Parquet format.")


DEFAULT_RNA_FILES = [
    "mrna.csv",
    "mirna.csv",
    "lncrna.csv",
    "circrna.csv",
]


@app.command()
def convert(
    input_dir: Path = typer.Argument(
        ...,
        help="Directory containing mrna.csv, mirna.csv, lncrna.csv and circrna.csv",
    ),
    output_dir: Path | None = typer.Option(
        None,
        "--output-dir",
        "-o",
        help="Directory to save parquet files. Defaults to input_dir.",
    ),
    files: List[str] = typer.Option(
        DEFAULT_RNA_FILES,
        "--file",
        "-f",
        help="CSV filename to convert. Can be specified multiple times.",
    ),
    encoding: str = typer.Option(
        "utf-8",
        "--encoding",
        help="CSV file encoding.",
    ),
    index: bool = typer.Option(
        False,
        "--index/--no-index",
        help="Whether to write DataFrame index into parquet.",
    ),
    overwrite: bool = typer.Option(
        True,
        "--overwrite/--no-overwrite",
        help="Whether to overwrite existing parquet files.",
    ),
):
    """
    Convert RNA CSV files in a directory to parquet files.
    """

    input_dir = input_dir.resolve()

    if output_dir is None:
        output_dir = input_dir
    else:
        output_dir = output_dir.resolve()

    if not input_dir.exists():
        raise typer.BadParameter(f"Input directory does not exist: {input_dir}")

    if not input_dir.is_dir():
        raise typer.BadParameter(f"Input path is not a directory: {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    converted_count = 0
    skipped_count = 0
    missing_count = 0

    for filename in files:
        csv_path = input_dir / filename

        if not csv_path.exists():
            typer.echo(f"[MISS] {csv_path}")
            missing_count += 1
            continue

        if not csv_path.is_file():
            typer.echo(f"[SKIP] Not a file: {csv_path}")
            skipped_count += 1
            continue

        parquet_name = csv_path.with_suffix(".parquet").name
        parquet_path = output_dir / parquet_name

        if parquet_path.exists() and not overwrite:
            typer.echo(f"[SKIP] Exists: {parquet_path}")
            skipped_count += 1
            continue

        try:
            df = pd.read_csv(csv_path, encoding=encoding)
            df.to_parquet(parquet_path, index=index)
        except Exception as exc:
            typer.echo(f"[ERROR] Failed to convert {csv_path}: {exc}")
            raise typer.Exit(code=1)

        typer.echo(f"[OK] {csv_path.name} -> {parquet_path}")
        converted_count += 1

    typer.echo("-" * 60)
    typer.echo(f"Converted: {converted_count}")
    typer.echo(f"Skipped:   {skipped_count}")
    typer.echo(f"Missing:   {missing_count}")


if __name__ == "__main__":
    app()
