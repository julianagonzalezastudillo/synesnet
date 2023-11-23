"""
=================================
            SYNESNET
=================================
This module is design to perform a Student t-test for independent samples between synesthets and control subjects.
In the case of coreness, the analysis is performed on a specific subgroup of nodes that exhibit a significant difference
between synesthetic and control subjects in terms of strength. Importantly, the coreness False Discovery Rate (FDR)
correction was exclusively applied to this subgroup.
"""

import os.path
import numpy as np
from scipy import stats
import pandas as pd
import statsmodels.stats.multitest as smt
from statsmodels.stats.multitest import multipletests
import statsmodels
from tools import load_net_metrics, load_xyz, load_node_names
from config import DATA_DIR, P_VAL


# CONSTANTS
SELECTION = True  # True: select base on strength significance
metric_list = ["strength", "coreness"]

# nodes nodes positions and names
xyz = load_xyz()
n_name, n_name_full = load_node_names()
lh_ind = [
    index for index, element in enumerate(np.array(n_name)) if element.endswith("_L")
]
rh_ind = [
    index for index, element in enumerate(np.array(n_name)) if element.endswith("_R")
]

idx_select = slice(None)
df_select = []
for net_key in metric_list:
    # Load net metrics for all subjects
    metric = (
        "coreness_norm_by_rand"
        if net_key == "coreness"
        else net_key
    )
    Xnet_syn, Xnet_ctr = load_net_metrics(metric, idx_select=idx_select)

    # Perform t-test, Student's t-test for independent samples
    t_val, p_val = stats.ttest_ind(Xnet_syn, Xnet_ctr)
    t_val[np.isnan(t_val)] = 0  # Replace NaNs with 0

    # Perform FDR correction
    p_val = p_val[idx_select]
    rejected, corrected_p_values, _, _ = multipletests(p_val, method="fdr_bh")
    statsmodels.stats.multitest.fdrcorrection(
        p_val, alpha=P_VAL, method="indep", is_sorted=False
    )

    # Corrected by hand FDR (Benjamini-Hochberg)
    rank = np.arange(len(p_val)) + 1
    rank_idx = np.argsort(p_val)
    N = len(p_val)
    p_val_corrected_ = np.sort(p_val) * N / rank
    p_val_corrected = np.zeros(len(p_val))
    p_val_corrected[rank_idx] = p_val_corrected_

    if net_key == "coreness":
        Xnet_syn, Xnet_ctr = load_net_metrics(net_key, idx_select=idx_select)

    # Create DataFrame for results
    df = pd.DataFrame(
        {
            "node": np.array(n_name)[idx_select],
            "node_complete": np.array(n_name_full)[idx_select],
            "node_idx": np.arange(len(n_name))[idx_select],
            "metric": np.array([net_key] * len(n_name))[idx_select],
            "metric_synes": Xnet_syn.mean(axis=0)[idx_select],
            "metric_ctr": Xnet_ctr.mean(axis=0)[idx_select],
            "t-val": np.array(t_val)[idx_select],
            "p-val": p_val,
            "p-val_corrected": p_val_corrected,
        }
    )

    # get indexes of significant nodes and append the DataFrame to the list
    if SELECTION and net_key == "strength":
        idx_select = np.array(df["node_idx"][(df["p-val"] < P_VAL)])
        df_select.append(df.loc[idx_select])
    elif SELECTION and net_key != "strength":
        idx_select = np.array(df["node_idx"][(df["p-val_corrected"] < P_VAL)])
        df_select.append(df)

# Concatenate the list of DataFrames into a single DataFrame and save to .csv file
df_stats = pd.concat(df_select, axis=0, ignore_index=True)
stats_file = os.path.join(DATA_DIR, "results", "stats_results.csv")
# df_stats.to_csv(stats_file, index=False)
