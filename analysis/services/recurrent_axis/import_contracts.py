from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping


@dataclass(frozen=True, slots=True)
class AxisContextSpec:
    dataset_source: str
    module: str
    dataset_name: str

    group_type: str
    group_by: str = ""
    group_value: str = ""

    annotation_dir_name: str = ""
    annotation_file_prefix: str = ""


@dataclass(frozen=True, slots=True)
class NormalizedAxisRow:
    row_index: int

    source_axis_id: str
    source_axis_type: str

    axis_key: str
    axis_signature: str
    axis_type: str

    miRNA: str
    mRNA: str
    lncRNA: str = ""
    circRNA: str = ""

    evidence: Mapping[str, Any] = field(default_factory=dict)
    extra_data: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class ParsedAxisArtifact:
    file_path: Path
    file_sha256: str

    result_kind: str
    schema_version: str

    rows: tuple[NormalizedAxisRow, ...]

    @property
    def row_count(self) -> int:
        return len(self.rows)


@dataclass(frozen=True, slots=True)
class AxisImportJob:
    context: AxisContextSpec
    result_kind: str
    file_path: Path
    schema_version: str = "v1"
