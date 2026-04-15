# PC-index

**PC-index** (Prevalence–CP10K index) is a metric for summarizing gene expression in single-cell RNA-seq data by integrating expression magnitude and cellular prevalence into a single value.

A PC-index of *j* means that at least *j%* of cells express the gene at **≥ j CP10K**.

## Definition

For a given gene, let the CP10K expression values across cells be sorted in decreasing order:

x(1) ≥ x(2) ≥ ... ≥ x(N)

For rank *k*, the corresponding prevalence is:

100 × k / N

The **continuous PC-index** is defined as:

PC-index = max over k of min(x(k), 100 × k / N)

This corresponds to the intersection between the ranked expression curve and the prevalence axis expressed as a percentage of cells.

## Input

A **cell-by-gene matrix** of CP10K-normalized expression values:

- rows = cells
- columns = genes
- values = CP10K
- first column = cell identifiers

Supported formats:
- CSV
- TSV

## Output

A gene-level summary table with:

- `gene`
- `n_cells`
- `mean_cp10k`
- `pct_positive`
- `pc_index`

## Python usage

Run from the command line:

```bash
python pc_index.py input_matrix.csv output_summary.csv
