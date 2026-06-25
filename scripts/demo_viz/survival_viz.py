import os
import argparse
import pandas as pd
import matplotlib

matplotlib.use("TkAgg")

import matplotlib.pyplot as plt

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test


def plot_cerna_survival(
        survival_data_file,
        title="ceRNA axis-based survival analysis"
):
    surv_df = pd.read_csv(survival_data_file)

    surv_df = surv_df.dropna(subset=["n_os", "c_os_status", "ceRNA_cluster"]).copy()

    surv_df["n_os"] = pd.to_numeric(surv_df["n_os"], errors="coerce")
    surv_df["c_os_status"] = pd.to_numeric(surv_df["c_os_status"], errors="coerce")

    surv_df = surv_df.dropna(subset=["n_os", "c_os_status"]).copy()
    surv_df["c_os_status"] = surv_df["c_os_status"].astype(int)

    kmf = KaplanMeierFitter()

    plt.figure(figsize=(6, 5))

    for group_name in ["Cluster_1", "Cluster_2"]:
        group_df = surv_df[surv_df["ceRNA_cluster"] == group_name]

        kmf.fit(
            durations=group_df["n_os"],
            event_observed=group_df["c_os_status"],
            label=f"{group_name} (n={group_df.shape[0]})"
        )

        kmf.plot_survival_function(ci_show=True)

    plt.title(f"{title}")

    plt.xlabel("Time days")
    plt.ylabel("Overall survival probability")
    plt.grid(alpha=0.3)
    plt.tight_layout()

    plt.show()


plot_cerna_survival(
    survival_data_file="E:\\Projects\\ceRNAxis\\ceRNAxis_api\\workspace\\3603a460-1972-4771-b3e0-21e05f00d9c8\\output"
                       "\\demo_task_survival_analysis.csv",
)
