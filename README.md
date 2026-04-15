# PC-index

PC-index (Prevalence–CP10K index) is a metric for summarizing gene expression in single-cell RNA-seq data by integrating expression magnitude and cellular prevalence into a single value.

## Definition

For each gene, the PC-index is defined as the largest value *j* such that at least *j%* of cells have expression ≥ *j* (CP10K).

## Input

A **cell-by-gene matrix** of CP10K-normalized expression values:

- rows = cells  
- columns = genes  
- values = CP10K (counts per 10,000)  
- first column = cell identifiers (optional)

## Usage

```bash
python code/summarize_pcindex_matrix.py input_matrix.csv output_summary.csv
```

## Output

A gene-level summary table including:

- number of cells  
- mean CP10K  
- percentage of expressing cells  
- PC-index  

## Example

```bash
python code/summarize_pcindex_matrix.py example_data/small_matrix.csv example_output/summary.csv
```

## Notes

- Input must be **CP10K-normalized** (not raw counts, not log-transformed).  
- The PC-index emphasizes genes with **jointly strong and widespread expression**.  

## Citation

Manuscript in preparation.
