from collections import defaultdict

from database.models import AxisSignatureProjectIndex

from analysis.services.axis_final.axis_signature import (
    add_axis_signature_to_records,
)


DEFAULT_MATCH_QUERY_CHUNK_SIZE = 2000


class AxisFinalProjectMatchError(ValueError):
    pass


def get_unique_axis_signatures_from_records(records: list[dict]) -> list[str]:
    signatures = {
        str(record.get("axis_signature") or "").strip()
        for record in records
        if str(record.get("axis_signature") or "").strip()
    }

    return sorted(signatures)


def serialize_axis_project_index_row(row: AxisSignatureProjectIndex) -> dict:
    """
    Serialize one AxisSignatureProjectIndex row for axis_final response.

    This is intentionally lightweight.
    Detailed log2FC/regulation values can be queried later by occurrence_id.
    """
    return {
        "project_id": row.project_id,
        "occurrence_id": row.occurrence_id,

        "source": row.source,
        "module": row.module,
        "dataset_name": row.dataset_name,
        "group_type": row.group_type,
        "group_by": row.group_by,

        "axis_id": row.axis_id,
        "axis_type": row.axis_type,
        "axis_regulation": row.axis_regulation,
    }


def query_project_matches_by_axis_signatures(
    *,
    axis_signatures: list[str],
    chunk_size: int = DEFAULT_MATCH_QUERY_CHUNK_SIZE,
) -> dict[str, list[dict]]:
    """
    Query AxisSignatureProjectIndex by axis_signature.

    Return:
        {
            axis_signature: [project_match, ...]
        }
    """
    if not axis_signatures:
        return {}

    matches_by_signature = defaultdict(list)

    qs = (
        AxisSignatureProjectIndex.objects
        .filter(axis_signature__in=axis_signatures)
        .order_by(
            "axis_signature",
            "source",
            "module",
            "dataset_name",
            "group_type",
            "group_by",
        )
        .only(
            "axis_signature",
            "project_id",
            "occurrence_id",
            "source",
            "module",
            "dataset_name",
            "group_type",
            "group_by",
            "axis_id",
            "axis_type",
            "axis_regulation",
        )
    )

    for row in qs.iterator(chunk_size=chunk_size):
        matches_by_signature[row.axis_signature].append(
            serialize_axis_project_index_row(row)
        )

    return dict(matches_by_signature)


def build_axis_project_match_summary(
    *,
    records: list[dict],
    matches_by_signature: dict[str, list[dict]],
) -> dict:
    matched_axis_count = 0
    matched_project_keys = set()
    project_hit_count_map = defaultdict(int)

    for record in records:
        signature = record.get("axis_signature", "")
        matches = matches_by_signature.get(signature, [])

        if matches:
            matched_axis_count += 1

        for match in matches:
            project_key = (
                match["source"],
                match["module"],
                match["dataset_name"],
                match["group_type"],
                match["group_by"],
            )

            matched_project_keys.add(project_key)
            project_hit_count_map[project_key] += 1

    project_hits = []

    for project_key, matched_axis_count_for_project in project_hit_count_map.items():
        source, module, dataset_name, group_type, group_by = project_key

        project_hits.append(
            {
                "source": source,
                "module": module,
                "dataset_name": dataset_name,
                "group_type": group_type,
                "group_by": group_by,
                "matched_axis_count": matched_axis_count_for_project,
            }
        )

    project_hits = sorted(
        project_hits,
        key=lambda item: (
            item["source"],
            item["module"],
            item["dataset_name"],
            item["group_type"],
            item["group_by"],
        ),
    )

    total_axis_count = len(records)

    return {
        "total_axis_count": total_axis_count,
        "matched_axis_count": matched_axis_count,
        "unmatched_axis_count": total_axis_count - matched_axis_count,
        "matched_project_count": len(matched_project_keys),
        "project_hits": project_hits,
    }


def attach_project_matches_to_axis_records(
    *,
    records: list[dict],
    max_matches_per_axis: int | None = None,
) -> dict:
    """
    Attach reference dataset project matches to serialized user axis_final records.

    User axis_final records are not written to database.

    Return:
        {
            "records": enriched_records,
            "summary": {...}
        }
    """
    if records is None:
        records = []

    records_with_signature = add_axis_signature_to_records(records)

    axis_signatures = get_unique_axis_signatures_from_records(
        records_with_signature
    )

    matches_by_signature = query_project_matches_by_axis_signatures(
        axis_signatures=axis_signatures,
    )

    summary = build_axis_project_match_summary(
        records=records_with_signature,
        matches_by_signature=matches_by_signature,
    )

    enriched_records = []

    for record in records_with_signature:
        signature = record.get("axis_signature", "")
        matches = matches_by_signature.get(signature, [])

        match_count = len(matches)

        if max_matches_per_axis is not None:
            matches_for_response = matches[:max_matches_per_axis]
        else:
            matches_for_response = matches

        enriched_records.append(
            {
                **record,
                "dataset_project_match_count": match_count,
                "dataset_project_matches": matches_for_response,
            }
        )

    return {
        "records": enriched_records,
        "summary": summary,
    }


def enrich_axis_final_response_with_project_matches(
    *,
    response_data: dict,
    max_matches_per_axis: int | None = None,
) -> dict:
    """
    Enrich build_axis_final_response_from_dataframe() response data.

    Adds:
        - axis_signature to every row
        - dataset_project_match_count to every row
        - dataset_project_matches to every row
        - axis_project_match_summary to top-level response
    """
    if not isinstance(response_data, dict):
        raise AxisFinalProjectMatchError(
            "response_data must be a dict."
        )

    records = response_data.get("results", [])

    match_result = attach_project_matches_to_axis_records(
        records=records,
        max_matches_per_axis=max_matches_per_axis,
    )

    return {
        **response_data,
        "results": match_result["records"],
        "axis_project_match_summary": match_result["summary"],
        "axis_project_match_enabled": True,
    }
