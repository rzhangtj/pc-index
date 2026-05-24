from pathlib import Path
import sys

import numpy as np
import pandas as pd


EXPRESSION_FILE = "data/expression_cp10k.csv"
METADATA_FILE = "data/metadata.csv"
OUTPUT_FILE = "output/supplementary_table1_top100_union.csv"

CELL_ID_COL = "cell_id"
SAMPLE_COL = "sample_id"
CONDITION_COL = "condition"
CLASS_COL = "class"
SUBCLASS_COL = "subclass"

CONTROL_LABEL = "Uninjured"
INJURY_LABEL = "1 WkPI"
MICROGLIA_CLASS_LABEL = "Microglia"
MICROGLIA_SUBCLASS_LABEL = "Microglia/Hematopoietic"


SCRIPT_DIR = Path(__file__).resolve().parent
for path in [SCRIPT_DIR, *SCRIPT_DIR.parents]:
    if (path / "pc_index.py").exists():
        sys.path.insert(0, str(path))
        break

from pc_index import compute_pc_index


def read_inputs():
    expression = pd.read_csv(SCRIPT_DIR / EXPRESSION_FILE)
    metadata = pd.read_csv(SCRIPT_DIR / METADATA_FILE)

    expression = expression.set_index(CELL_ID_COL)
    metadata = metadata.set_index(CELL_ID_COL)
    expression.index = expression.index.astype(str)
    metadata.index = metadata.index.astype(str)

    shared_cells = expression.index.intersection(metadata.index)
    if shared_cells.empty:
        raise ValueError("No matching cell IDs between expression and metadata.")

    expression = expression.loc[shared_cells].apply(pd.to_numeric, errors="coerce").fillna(0)
    metadata = metadata.loc[shared_cells]
    return expression, metadata


def select_cells(expression, metadata):
    keep = (
        metadata[CLASS_COL].eq(MICROGLIA_CLASS_LABEL)
        & metadata[SUBCLASS_COL].eq(MICROGLIA_SUBCLASS_LABEL)
        & metadata[CONDITION_COL].isin([CONTROL_LABEL, INJURY_LABEL])
    )
    if not keep.any():
        raise ValueError("No cells remain after microglia and condition filtering.")
    return expression.loc[keep], metadata.loc[keep]


def summarize_samples(expression, metadata):
    records = []
    for (condition, sample), sample_meta in metadata.groupby([CONDITION_COL, SAMPLE_COL]):
        sample_expression = expression.loc[sample_meta.index]
        records.append(
            pd.DataFrame(
                {
                    "gene": sample_expression.columns,
                    "condition": condition,
                    "sample_id": sample,
                    "pc_index": sample_expression.apply(
                        lambda x: compute_pc_index(x.to_numpy()), axis=0
                    ).to_numpy(),
                    "mean_cp10k": sample_expression.mean(axis=0).to_numpy(),
                    "pct_positive": (100 * sample_expression.gt(0).mean(axis=0)).to_numpy(),
                }
            )
        )
    return pd.concat(records, ignore_index=True)


def calculate_gene_metrics(sample_metrics):
    means = (
        sample_metrics.groupby(["gene", "condition"])[["pc_index", "mean_cp10k", "pct_positive"]]
        .mean()
        .unstack("condition")
    )

    for label in [CONTROL_LABEL, INJURY_LABEL]:
        if label not in means.columns.get_level_values("condition"):
            raise ValueError(f"Missing condition: {label}")

    table = pd.DataFrame(index=means.index)
    table["pc_index_uninjured"] = means[("pc_index", CONTROL_LABEL)]
    table["pc_index_1wkpi"] = means[("pc_index", INJURY_LABEL)]
    table["delta_pc_index"] = table["pc_index_1wkpi"] - table["pc_index_uninjured"]
    table["mean_cp10k_uninjured"] = means[("mean_cp10k", CONTROL_LABEL)]
    table["mean_cp10k_1wkpi"] = means[("mean_cp10k", INJURY_LABEL)]
    table["delta_mean_cp10k"] = table["mean_cp10k_1wkpi"] - table["mean_cp10k_uninjured"]
    table["pct_positive_uninjured"] = means[("pct_positive", CONTROL_LABEL)]
    table["pct_positive_1wkpi"] = means[("pct_positive", INJURY_LABEL)]
    table["delta_pct_positive"] = (
        table["pct_positive_1wkpi"] - table["pct_positive_uninjured"]
    )

    table = table.reset_index()
    table["rank_delta_pc_index"] = table["delta_pc_index"].rank(
        method="min", ascending=False
    ).astype(int)
    table["rank_delta_mean_cp10k"] = table["delta_mean_cp10k"].rank(
        method="min", ascending=False
    ).astype(int)
    table["rank_delta_pct_positive"] = table["delta_pct_positive"].rank(
        method="min", ascending=False
    ).astype(int)
    return table


def add_top100_membership(table):
    n_top = min(100, len(table))
    top_pc = set(table.nsmallest(n_top, "rank_delta_pc_index")["gene"])
    top_mean = set(table.nsmallest(n_top, "rank_delta_mean_cp10k")["gene"])
    top_pct = set(table.nsmallest(n_top, "rank_delta_pct_positive")["gene"])

    table = table.copy()
    table["top100_delta_pc_index"] = table["gene"].isin(top_pc)
    table["top100_delta_mean_cp10k"] = table["gene"].isin(top_mean)
    table["top100_delta_pct_positive"] = table["gene"].isin(top_pct)
    return table


def membership_group(row):
    pc = row["top100_delta_pc_index"]
    mean = row["top100_delta_mean_cp10k"]
    pct = row["top100_delta_pct_positive"]
    if pc and mean and pct:
        return "all_three"
    if pc and mean:
        return "pc_index_and_mean"
    if pc and pct:
        return "pc_index_and_pct_positive"
    if mean and pct:
        return "mean_and_pct_positive"
    if pc:
        return "pc_index_only"
    if mean:
        return "mean_only"
    return "pct_positive_only"


def make_supplementary_table(table):
    table = add_top100_membership(table)
    table = table[
        table["top100_delta_pc_index"]
        | table["top100_delta_mean_cp10k"]
        | table["top100_delta_pct_positive"]
    ].copy()
    table["membership_group"] = table.apply(membership_group, axis=1)

    columns = [
        "gene",
        "pc_index_uninjured",
        "pc_index_1wkpi",
        "delta_pc_index",
        "mean_cp10k_uninjured",
        "mean_cp10k_1wkpi",
        "delta_mean_cp10k",
        "pct_positive_uninjured",
        "pct_positive_1wkpi",
        "delta_pct_positive",
        "rank_delta_pc_index",
        "rank_delta_mean_cp10k",
        "rank_delta_pct_positive",
        "top100_delta_pc_index",
        "top100_delta_mean_cp10k",
        "top100_delta_pct_positive",
        "membership_group",
    ]
    return table[columns].sort_values(["membership_group", "rank_delta_pc_index", "gene"])


def print_overlap_summary(table):
    counts = table["membership_group"].value_counts()
    labels = [
        ("total union size", len(table)),
        ("number shared by all three", counts.get("all_three", 0)),
        ("PC-index only", counts.get("pc_index_only", 0)),
        ("mean only", counts.get("mean_only", 0)),
        ("percent-positive only", counts.get("pct_positive_only", 0)),
        ("PC-index + mean only", counts.get("pc_index_and_mean", 0)),
        ("PC-index + percent-positive only", counts.get("pc_index_and_pct_positive", 0)),
        ("mean + percent-positive only", counts.get("mean_and_pct_positive", 0)),
    ]
    for label, value in labels:
        print(f"{label}: {value}")


def main():
    expression, metadata = read_inputs()
    expression, metadata = select_cells(expression, metadata)
    sample_metrics = summarize_samples(expression, metadata)
    gene_metrics = calculate_gene_metrics(sample_metrics)
    supplementary_table = make_supplementary_table(gene_metrics)

    output_path = SCRIPT_DIR / OUTPUT_FILE
    output_path.parent.mkdir(exist_ok=True)
    supplementary_table.to_csv(output_path, index=False)
    print_overlap_summary(supplementary_table)
    print(f"wrote: {output_path}")


if __name__ == "__main__":
    main()
