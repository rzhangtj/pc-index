# PC-index

**PC-index** (Prevalence-CP10K index) is a gene-level summary metric for single-cell RNA-seq data. It integrates expression magnitude and cellular prevalence into a single interpretable value.

A PC-index of *j* means that at least *j*% of cells express the gene at **>= *j* CP10K**.

## Definition

For a given gene, let CP10K-normalized expression values across *N* cells be sorted in decreasing order:

```text
x(1) >= x(2) >= ... >= x(N)
```

At rank *k*, the corresponding cellular prevalence is:

```text
100 * k / N
```

The continuous PC-index is defined as:

```text
PC-index = max over k of min(x(k), 100 * k / N)
```

Genes with no detectable expression are assigned PC-index = 0.

## Core PC-index calculation

The input should be a cell-by-gene matrix of CP10K-normalized expression values:

- rows = cells
- columns = genes
- values = CP10K
- first column = cell identifiers

Supported formats are CSV and TSV.

The output is a gene-level table with:

- `gene`
- `n_cells`
- `mean_cp10k`
- `pct_positive`
- `pc_index`

Run with Python:

```bash
python code/pc_index.py input.csv output.csv
```

Run with R:

```bash
Rscript code/pc_index.R input.csv output.csv
```

The PC-index is intended as a compact summary of expression magnitude and cellular prevalence. It complements, rather than replaces, standard summaries such as mean expression and percentage of expressing cells. Input expression values should be normalized consistently across cells before calculating PC-index. In the manuscript, CP10K-normalized values were used.

## GSE172167 spinal cord injury analysis

This repository also includes code used to generate **Supplementary Table 1** for the PC-index manuscript.

Supplementary Table 1 lists the union of genes ranked in the top 100 by injury-associated change in:

- Delta PC-index
- Delta mean CP10K
- Delta % positive cells

The analysis uses nuclei from GSE172167 annotated as `class = Microglia` and `subclass = Microglia/Hematopoietic`. The comparison is **1 WkPI spinal cord injury minus uninjured control**.

Raw data are publicly available from GEO accession **GSE172167**. Full processed expression data are not included in this repository. Small real-data example files are provided only to illustrate input format.

## GSE172167 input files

The Supplementary Table 1 script expects two input files:

- `example_data/expression_cp10k.csv`
- `example_data/metadata.csv`

The expression file should have nuclei/cells as rows, genes as columns, and `cell_id` as the first column.

Example expression file:

```text
cell_id,GeneA,GeneB,GeneC
cell_001,0,2.5,10.1
cell_002,1.2,0,8.4
cell_003,0,3.1,0
```

The metadata file should have one row per nucleus/cell and include:

- `cell_id`
- `sample_id`
- `condition`
- `class`
- `subclass`

Example metadata file:

```text
cell_id,sample_id,condition,class,subclass
cell_001,S1,Uninjured,Microglia,Microglia/Hematopoietic
cell_002,S1,Uninjured,Microglia,Microglia/Hematopoietic
cell_003,S2,1 WkPI,Microglia,Microglia/Hematopoietic
```

Default condition labels are `Uninjured` and `1 WkPI`. Column names, condition labels, and annotation labels can be edited near the top of:

```text
code/make_supplementary_table1.py
```

## Run the GSE172167 Supplementary Table 1 script

From the repository root:

```bash
python code/make_supplementary_table1.py
```

The output table is written to:

```text
output/supplementary_table1_top100_union.csv
```

The script calculates, for each gene:

- PC-index in uninjured control
- PC-index at 1 WkPI
- 1 WkPI-minus-uninjured change in PC-index
- mean CP10K in uninjured control
- mean CP10K at 1 WkPI
- 1 WkPI-minus-uninjured change in mean CP10K
- percentage-positive cells in uninjured control
- percentage-positive cells at 1 WkPI
- 1 WkPI-minus-uninjured change in percentage-positive cells
- ranks by each change metric
- top-100 membership for each change metric
- top-100 membership group

## Repository structure

```text
public/
  code/
    pc_index.py
    pc_index.R
    make_supplementary_table1.py

  example_data/
    expression_cp10k.csv
    metadata.csv

  
  README.md
  requirements.txt
```

## Requirements

Python dependencies:

```text
numpy
pandas
scipy
```

The core PC-index script requires `numpy` and `pandas`. The GSE172167 Supplementary Table 1 script may also use `scipy`.

## Citation

Zhang R. **PC-index: a metric integrating expression magnitude and cellular prevalence in single-cell RNA-seq.**
