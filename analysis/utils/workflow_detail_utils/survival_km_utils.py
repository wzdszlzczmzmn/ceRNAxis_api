import math

import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test


DEFAULT_SURVIVAL_TIME_COL = "n_os"
DEFAULT_SURVIVAL_EVENT_COL = "c_os_status"
DEFAULT_SURVIVAL_GROUP_COL = "ceRNA_cluster"

DEFAULT_SURVIVAL_GROUPS = [
    "Cluster_1",
    "Cluster_2",
]

SURVIVAL_REQUIRED_COLUMNS = {
    DEFAULT_SURVIVAL_TIME_COL,
    DEFAULT_SURVIVAL_EVENT_COL,
    DEFAULT_SURVIVAL_GROUP_COL,
}


class SurvivalKMInputError(ValueError):
    pass


def safe_float_or_none(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None

    if pd.isna(result) or not math.isfinite(result):
        return None

    return result


def validate_survival_dataframe_columns(
    df: pd.DataFrame,
    required_columns: set[str] | None = None,
) -> None:
    required_columns = required_columns or SURVIVAL_REQUIRED_COLUMNS
    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise SurvivalKMInputError(
            "Survival analysis file is missing required column(s): "
            f"{', '.join(sorted(missing_columns))}."
        )


def normalize_survival_dataframe(
    df: pd.DataFrame,
    time_col: str = DEFAULT_SURVIVAL_TIME_COL,
    event_col: str = DEFAULT_SURVIVAL_EVENT_COL,
    group_col: str = DEFAULT_SURVIVAL_GROUP_COL,
    valid_groups: list[str] | None = None,
) -> pd.DataFrame:
    valid_groups = valid_groups or DEFAULT_SURVIVAL_GROUPS

    normalized_df = df.copy()

    normalized_df[time_col] = pd.to_numeric(
        normalized_df[time_col],
        errors="coerce",
    )

    normalized_df[event_col] = pd.to_numeric(
        normalized_df[event_col],
        errors="coerce",
    )

    normalized_df[group_col] = (
        normalized_df[group_col]
        .astype(str)
        .str.strip()
    )

    normalized_df = normalized_df.dropna(
        subset=[
            time_col,
            event_col,
            group_col,
        ]
    ).copy()

    normalized_df = normalized_df[
        normalized_df[group_col].isin(valid_groups)
    ].copy()

    normalized_df = normalized_df[
        normalized_df[time_col] >= 0
    ].copy()

    normalized_df = normalized_df[
        normalized_df[event_col].isin([0, 1])
    ].copy()

    normalized_df[event_col] = normalized_df[event_col].astype(int)

    return normalized_df


def dataframe_to_km_points(
    survival_df: pd.DataFrame,
    time_col: str,
    value_col: str,
) -> list[dict]:
    points = []

    for _, row in survival_df.iterrows():
        time_value = safe_float_or_none(row.get(time_col))
        survival_value = safe_float_or_none(row.get(value_col))

        if time_value is None or survival_value is None:
            continue

        points.append(
            {
                "time": time_value,
                "survival": survival_value,
            }
        )

    return points


def build_single_km_group_data(
    group_name: str,
    group_df: pd.DataFrame,
    time_col: str = DEFAULT_SURVIVAL_TIME_COL,
    event_col: str = DEFAULT_SURVIVAL_EVENT_COL,
) -> dict | None:
    if group_df.empty:
        return None

    kmf = KaplanMeierFitter()

    kmf.fit(
        durations=group_df[time_col],
        event_observed=group_df[event_col],
        label=group_name,
    )

    survival_df = kmf.survival_function_.reset_index()
    survival_time_col = survival_df.columns[0]
    survival_value_col = survival_df.columns[1]

    ci_df = kmf.confidence_interval_survival_function_.reset_index()
    ci_time_col = ci_df.columns[0]
    ci_lower_col = ci_df.columns[1]
    ci_upper_col = ci_df.columns[2]

    event_count = int(group_df[event_col].sum())
    sample_count = int(group_df.shape[0])

    return {
        "name": group_name,
        "n": sample_count,
        "event_count": event_count,
        "censored_count": sample_count - event_count,
        "points": dataframe_to_km_points(
            survival_df=survival_df,
            time_col=survival_time_col,
            value_col=survival_value_col,
        ),
        "ci_lower": dataframe_to_km_points(
            survival_df=ci_df,
            time_col=ci_time_col,
            value_col=ci_lower_col,
        ),
        "ci_upper": dataframe_to_km_points(
            survival_df=ci_df,
            time_col=ci_time_col,
            value_col=ci_upper_col,
        ),
    }


def calculate_logrank_p_value(
    surv_df: pd.DataFrame,
    group_a: str,
    group_b: str,
    time_col: str = DEFAULT_SURVIVAL_TIME_COL,
    event_col: str = DEFAULT_SURVIVAL_EVENT_COL,
    group_col: str = DEFAULT_SURVIVAL_GROUP_COL,
) -> float | None:
    group_a_df = surv_df[surv_df[group_col] == group_a]
    group_b_df = surv_df[surv_df[group_col] == group_b]

    if group_a_df.empty or group_b_df.empty:
        return None

    try:
        result = logrank_test(
            group_a_df[time_col],
            group_b_df[time_col],
            event_observed_A=group_a_df[event_col],
            event_observed_B=group_b_df[event_col],
        )
    except Exception:
        return None

    return safe_float_or_none(result.p_value)


def build_survival_km_data_from_dataframe(
    *,
    task,
    survival_file_name: str,
    df: pd.DataFrame,
    title: str,
    time_col: str = DEFAULT_SURVIVAL_TIME_COL,
    event_col: str = DEFAULT_SURVIVAL_EVENT_COL,
    group_col: str = DEFAULT_SURVIVAL_GROUP_COL,
    valid_groups: list[str] | None = None,
    x_label: str = "Time days",
    y_label: str = "Overall survival probability",
) -> dict:
    return build_survival_km_data_from_dataframe_common(
        survival_file_name=survival_file_name,
        df=df,
        title=title,
        base_response={
            "uuid": str(task.uuid),
            "task_name": task.task_name,
        },
        time_col=time_col,
        event_col=event_col,
        group_col=group_col,
        valid_groups=valid_groups,
        x_label=x_label,
        y_label=y_label,
    )


def build_survival_km_data_from_dataframe_common(
    *,
    survival_file_name: str,
    df: pd.DataFrame,
    title: str,
    base_response: dict | None = None,
    time_col: str = DEFAULT_SURVIVAL_TIME_COL,
    event_col: str = DEFAULT_SURVIVAL_EVENT_COL,
    group_col: str = DEFAULT_SURVIVAL_GROUP_COL,
    valid_groups: list[str] | None = None,
    x_label: str = "Time days",
    y_label: str = "Overall survival probability",
) -> dict:
    valid_groups = valid_groups or DEFAULT_SURVIVAL_GROUPS

    validate_survival_dataframe_columns(
        df=df,
        required_columns={time_col, event_col, group_col},
    )

    surv_df = normalize_survival_dataframe(
        df=df,
        time_col=time_col,
        event_col=event_col,
        group_col=group_col,
        valid_groups=valid_groups,
    )

    raw_count = int(df.shape[0])
    cleaned_count = int(surv_df.shape[0])
    dropped_count = raw_count - cleaned_count

    groups = []

    for group_name in valid_groups:
        group_df = surv_df[
            surv_df[group_col] == group_name
        ].copy()

        group_data = build_single_km_group_data(
            group_name=group_name,
            group_df=group_df,
            time_col=time_col,
            event_col=event_col,
        )

        if group_data is not None:
            groups.append(group_data)

    logrank_p = None

    if len(valid_groups) >= 2:
        logrank_p = calculate_logrank_p_value(
            surv_df=surv_df,
            group_a=valid_groups[0],
            group_b=valid_groups[1],
            time_col=time_col,
            event_col=event_col,
            group_col=group_col,
        )

    response_data = {
        "survival_file": survival_file_name,
        "title": title,
        "x_label": x_label,
        "y_label": y_label,
        "summary": {
            "raw_count": raw_count,
            "cleaned_count": cleaned_count,
            "dropped_count": dropped_count,
            "group_count": len(groups),
            "logrank_p": logrank_p,
        },
        "groups": groups,
    }

    if base_response:
        response_data = {
            **base_response,
            **response_data,
        }

    return response_data
