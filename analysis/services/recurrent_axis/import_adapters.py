from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

import pandas as pd

from database.models import AxisResultKind

from .import_contracts import (
    NormalizedAxisRow,
    ParsedAxisArtifact,
)
from .import_normalization import (
    AxisImportValidationError,
    build_axis_key,
    build_axis_signature,
    calculate_file_sha256,
    derive_canonical_axis_type,
    normalize_optional_float,
    normalize_text,
    to_json_value,
)


BASE_COLUMNS = {
    "axis_id",
    "axis_type",
    "mRNA",
    "miRNA",
    "lncRNA",
    "circRNA",
}


class BaseAxisResultAdapter(ABC):
    result_kind: str
    schema_version = "v1"

    required_columns = {
        "mRNA",
        "miRNA",
    }

    optional_columns = {
        "axis_id",
        "axis_type",
        "lncRNA",
        "circRNA",
    }

    @property
    @abstractmethod
    def evidence_columns(self) -> set[str]:
        raise NotImplementedError

    @abstractmethod
    def build_evidence(
        self,
        *,
        row_index: int,
        row: pd.Series,
    ) -> dict[str, Any]:
        raise NotImplementedError

    def validate_result_schema(
        self,
        df: pd.DataFrame,
    ) -> None:
        pass

    def get_empty_result_columns(self) -> list[str]:
        """
        Return a complete synthetic schema for a valid zero-row result.

        Some upstream workflows use a zero-byte file to mean that the
        analysis completed successfully but produced no Axis rows. The
        importer converts that representation into an empty DataFrame with
        the adapter's expected columns, allowing the file to be persisted as
        an active artifact with row_count=0.
        """
        return sorted(
            BASE_COLUMNS
            | self.required_columns
            | self.optional_columns
            | self.evidence_columns
        )

    def read_dataframe(
        self,
        file_path: Path,
    ) -> pd.DataFrame:
        try:
            file_size = file_path.stat().st_size
        except OSError as exc:
            raise AxisImportValidationError(
                f"Unable to inspect CSV file: {file_path}"
            ) from exc

        if file_size == 0:
            return pd.DataFrame(
                columns=self.get_empty_result_columns(),
                dtype=object,
            )

        try:
            df = pd.read_csv(
                file_path,
                dtype=object,
            )
        except pd.errors.EmptyDataError as exc:
            raise AxisImportValidationError(
                "CSV file contains no parseable header: "
                f"{file_path}. "
                "A zero-byte file is supported as an explicit zero-row "
                "result, but a non-zero blank or whitespace-only file is "
                "invalid."
            ) from exc
        except Exception as exc:
            raise AxisImportValidationError(
                f"Unable to read CSV file: {file_path}. "
                f"Cause: {type(exc).__name__}: {exc}"
            ) from exc

        columns = [
            str(column).strip()
            for column in df.columns
        ]

        if len(columns) != len(set(columns)):
            raise AxisImportValidationError(
                "CSV contains duplicate column names "
                f"after trimming: {file_path}"
            )

        df.columns = columns

        return df

    def parse_file(
        self,
        *,
        file_path: Path,
        schema_version: str | None = None,
    ) -> ParsedAxisArtifact:
        file_path = Path(file_path).resolve()

        if not file_path.is_file():
            raise AxisImportValidationError(
                f"Result file is unavailable: {file_path}"
            )

        df = self.read_dataframe(file_path)

        missing_columns = sorted(
            self.required_columns - set(df.columns)
        )

        if missing_columns:
            raise AxisImportValidationError(
                "Missing required columns: "
                + ", ".join(missing_columns)
            )

        self.validate_result_schema(df)

        known_columns = (
            BASE_COLUMNS
            | self.required_columns
            | self.optional_columns
            | self.evidence_columns
        )

        rows = []
        first_row_by_axis_key = {}

        for row_index, row in df.iterrows():
            normalized_row = self.normalize_row(
                row_index=int(row_index),
                row=row,
                known_columns=known_columns,
            )

            previous_row_index = (
                first_row_by_axis_key.get(
                    normalized_row.axis_key
                )
            )

            if previous_row_index is not None:
                raise AxisImportValidationError(
                    "Duplicate structural Axis in one artifact: "
                    f"rows {previous_row_index} and {row_index}; "
                    f"signature={normalized_row.axis_signature!r}."
                )

            first_row_by_axis_key[
                normalized_row.axis_key
            ] = int(row_index)

            rows.append(normalized_row)

        return ParsedAxisArtifact(
            file_path=file_path,
            file_sha256=calculate_file_sha256(file_path),
            result_kind=self.result_kind,
            schema_version=(
                schema_version or self.schema_version
            ),
            rows=tuple(rows),
        )

    def normalize_row(
        self,
        *,
        row_index: int,
        row: pd.Series,
        known_columns: set[str],
    ) -> NormalizedAxisRow:
        miRNA = normalize_text(row.get("miRNA"))
        mRNA = normalize_text(row.get("mRNA"))
        lncRNA = normalize_text(row.get("lncRNA"))
        circRNA = normalize_text(row.get("circRNA"))

        axis_type = derive_canonical_axis_type(
            miRNA=miRNA,
            mRNA=mRNA,
            lncRNA=lncRNA,
            circRNA=circRNA,
            row_index=row_index,
        )

        axis_signature = build_axis_signature(
            miRNA=miRNA,
            mRNA=mRNA,
            lncRNA=lncRNA,
            circRNA=circRNA,
        )

        extra_data = {
            str(column): to_json_value(row.get(column))
            for column in row.index
            if str(column) not in known_columns
        }

        return NormalizedAxisRow(
            row_index=row_index,
            source_axis_id=normalize_text(
                row.get("axis_id")
            ),
            source_axis_type=normalize_text(
                row.get("axis_type")
            ),
            axis_key=build_axis_key(axis_signature),
            axis_signature=axis_signature,
            axis_type=axis_type,
            miRNA=miRNA,
            mRNA=mRNA,
            lncRNA=lncRNA,
            circRNA=circRNA,
            evidence=self.build_evidence(
                row_index=row_index,
                row=row,
            ),
            extra_data=extra_data,
        )


class AxisFinalResultAdapter(
    BaseAxisResultAdapter
):
    result_kind = AxisResultKind.AXIS_FINAL

    @property
    def evidence_columns(self) -> set[str]:
        return {
            "axis_regulation",
            "mRNA_log2FC",
            "mRNA_regulation",
            "miRNA_log2FC",
            "miRNA_regulation",
            "lncRNA_log2FC",
            "lncRNA_regulation",
            "circRNA_log2FC",
            "circRNA_regulation",
        }

    def validate_result_schema(
        self,
        df: pd.DataFrame,
    ) -> None:
        if not (
            self.evidence_columns & set(df.columns)
        ):
            raise AxisImportValidationError(
                "The file does not contain any "
                "Axis Final evidence columns."
            )

    def build_evidence(
        self,
        *,
        row_index: int,
        row: pd.Series,
    ) -> dict[str, Any]:
        return {
            "axis_regulation": normalize_text(
                row.get("axis_regulation")
            ),

            "mRNA_log2FC": normalize_optional_float(
                row.get("mRNA_log2FC"),
                field_name="mRNA_log2FC",
                row_index=row_index,
            ),
            "mRNA_regulation": normalize_text(
                row.get("mRNA_regulation")
            ),

            "miRNA_log2FC": normalize_optional_float(
                row.get("miRNA_log2FC"),
                field_name="miRNA_log2FC",
                row_index=row_index,
            ),
            "miRNA_regulation": normalize_text(
                row.get("miRNA_regulation")
            ),

            "lncRNA_log2FC": normalize_optional_float(
                row.get("lncRNA_log2FC"),
                field_name="lncRNA_log2FC",
                row_index=row_index,
            ),
            "lncRNA_regulation": normalize_text(
                row.get("lncRNA_regulation")
            ),

            "circRNA_log2FC": normalize_optional_float(
                row.get("circRNA_log2FC"),
                field_name="circRNA_log2FC",
                row_index=row_index,
            ),
            "circRNA_regulation": normalize_text(
                row.get("circRNA_regulation")
            ),
        }


class SpongeResultAdapter(
    BaseAxisResultAdapter
):
    result_kind = AxisResultKind.SPONGE

    required_columns = (
        BaseAxisResultAdapter.required_columns
        | {
            "cor",
            "pcor",
            "mscor",
        }
    )

    @property
    def evidence_columns(self) -> set[str]:
        return {
            "cor",
            "pcor",
            "mscor",
        }

    def build_evidence(
        self,
        *,
        row_index: int,
        row: pd.Series,
    ) -> dict[str, Any]:
        return {
            "cor": normalize_optional_float(
                row.get("cor"),
                field_name="cor",
                row_index=row_index,
            ),
            "pcor": normalize_optional_float(
                row.get("pcor"),
                field_name="pcor",
                row_index=row_index,
            ),
            "mscor": normalize_optional_float(
                row.get("mscor"),
                field_name="mscor",
                row_index=row_index,
            ),
        }


ADAPTERS = {
    AxisResultKind.AXIS_FINAL:
        AxisFinalResultAdapter(),

    AxisResultKind.SPONGE:
        SpongeResultAdapter(),
}


def get_axis_result_adapter(
    result_kind: str,
) -> BaseAxisResultAdapter:
    try:
        return ADAPTERS[result_kind]
    except KeyError as exc:
        raise AxisImportValidationError(
            f"Unsupported Axis result kind: "
            f"{result_kind!r}."
        ) from exc
