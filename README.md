# PC-index

**PC-index** (Prevalence–CP10K index) summarizes gene expression in single-cell RNA-seq by integrating expression magnitude and cellular prevalence into a single value.

A PC-index of *j* means that at least *j%* of cells express the gene at **≥ j CP10K**.

---

## Definition

For a given gene, let expression values across cells be sorted in decreasing order:

x(1) ≥ x(2) ≥ ... ≥ x(N)

At rank *k*, the corresponding prevalence is:

100 × k / N

The **continuous PC-index** is:

PC-index = max over k of min(x(k), 100 × k / N)

---

## Input

A **cell-by-gene matrix** of CP10K-normalized expression values:

- rows = cells  
- columns = genes  
- values = CP10K  
- first column = cell identifiers  

Supported formats: CSV, TSV

---

## Output

A gene-level table with:

- `gene`
- `n_cells`
- `mean_cp10k`
- `pct_positive`
- `pc_index`

---

## Usage

Command line:

```bash
python pc_index.py input.csv output.csv
Rscript pc_index.R input.csv output.csv
