#!/usr/bin/env python3
"""
Compute PC-index, mean CP10K, and percent-positive for each gene
from a cell-by-gene CP10K matrix.

Definition used here:
PC-index = max_k min(sorted_expression_k, 100 * k / n_cells)

This is the continuous version of the prevalence-CP10K intersection.
For example, a PC-index of 21 means roughly that at least 21% of cells
express the gene at >= 21 CP10K.

Input:
- CSV or TSV
- default: rows = cells, columns = genes
- values = CP10K

Output:
- CSV with one row per gene and columns:
    gene, n_cells, mean_cp10k, pct_positive, pc_index
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def compute_pc_index(values: np.ndarray) -> float:
    """
    Continuous PC-index:
    max over ranks of min(expression_at_rank, prevalence_at_rank_percent)

    Parameters
    ----------
    values : np.ndarray
        1D array of CP10K values across cells for one gene.

    Returns
    -------
    float
        Continuous PC-index.
    """
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]

    if x.size == 0:
        return np.nan

    x = np.sort(x)[::-1]  # descending
    prevalence_pct = 100.0 * np.arange(1, x.size + 1) / x.size
    return float(np.max(np.minimum(x, prevalence_pct)))


def summarize_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize each gene (column) in a cell-by-gene CP10K matrix.
    """
    results = []

    n_cells = len(df)

    for gene in df.columns:
        values = pd.to_numeric(df[gene], errors="coerce").to_numpy(dtype=float)

        mean_cp10k = float(np.nanmean(values)) if values.size else np.nan
        pct_positive = float(np.mean(values > 0) * 100.0) if values.size else np.nan
        pc_index = compute_pc_index(values)

        results.append(
            {
                "gene": gene,
                "n_cells": n_cells,
                "mean_cp10k": mean_cp10k,
                "pct_positive": pct_positive,
                "pc_index": pc_index,
            }
        )

    out = pd.DataFrame(results)
    out = out.sort_values("pc_index", ascending=False).reset_index(drop=True)
    return out


def read_table(path: Path, sep: str | None = None) -> pd.DataFrame:
    """
    Read CSV or TSV automatically unless --sep is specified.
    """
    if sep is None:
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep, index_col=0)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute PC-index, mean CP10K, and percent-positive for each gene."
    )
    parser.add_argument("input", type=Path, help="Input CSV/TSV file")
    parser.add_argument("output", type=Path, help="Output CSV file")
    parser.add_argument(
        "--sep",
        default=None,
        help="Input delimiter. Default: auto-detect by file extension.",
    )
    parser.add_argument(
        "--transpose",
        action="store_true",
        help="Use this if input is gene-by-cell instead of cell-by-gene.",
    )

    args = parser.parse_args()

    df = read_table(args.input, sep=args.sep)

    if args.transpose:
        df = df.T

    summary = summarize_matrix(df)
    summary.to_csv(args.output, index=False)

    print(f"Done. Wrote {len(summary)} genes to: {args.output}")


if __name__ == "__main__":
    main()
