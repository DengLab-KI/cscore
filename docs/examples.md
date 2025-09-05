# C-score Examples

This page provides comprehensive examples of how to use `cscore` with the datasets provided in the `testdata/` directory.

## Available Test Datasets

The repository includes three main test datasets:

1. **GSE237099** - Bulk RNA-seq data (unloading/reloading muscle study)
2. **GSE236519** - Single-cell RNA-seq data (neuronal cell types)
3. **simulation** - Simulated C-score data for testing

---

## 1. Bulk RNA-seq Example (GSE237099)

This example demonstrates C-score analysis using bulk RNA-seq data from a muscle unloading/reloading study.

### Basic Usage

```bash
cscore -i testdata/GSE237099 \
  -a GSE237099_1_unloading_reloading_Reloading_vs_Control_deseq2.txt \
  -b GSE237099_1_unloading_reloading_Unloading_vs_Control_deseq2.txt \
  -o output_bulk_example.tsv \
  -n ensembl_gene_id \
  -e log2FoldChange \
  -f padj
```

### What this example does:
- Compares gene expression between **Reloading vs Control** and **Unloading vs Control**
- Uses `ensembl_gene_id` as the gene identifier column
- Uses `log2FoldChange` as the effect size column
- Uses `padj` as the adjusted p-value column
- Identifies genes with common vs. divergent responses to unloading and reloading

### Expected Output:
The output file will contain columns including:
- `score`: C-score (positive = common direction; negative = divergent)
- `p`: permutation p-value for commonness/divergence
- `q_value`: Benjamini–Hochberg adjusted p-value
- `convergence`: "high" when commonness dominates, "low" otherwise
- All original columns from both input files with `_comp1`/`_comp2` suffixes

---

## 2. Single-cell RNA-seq Examples (GSE236519)

This dataset contains differential expression results for different neuronal cell types comparing knockout (KO) vs. knockdown (KD) conditions.

### Example 1: Deep Layer Neurons

```bash
cscore -i testdata/GSE236519 \
  -a KO_Deep_Layer_neurons.txt \
  -b KD_Deep_Layer_neurons.txt \
  -o output_deep_layer_neurons.tsv \
  -n gene \
  -e logFC \
  -f FDR
```

### Example 2: Interneurons

```bash
cscore -i testdata/GSE236519 \
  -a KO_Interneurons.txt \
  -b KD_Interneurons.txt \
  -o output_interneurons.tsv \
  -n gene \
  -e logFC \
  -f FDR
```

### Example 3: Superficial Layer Neurons

```bash
cscore -i testdata/GSE236519 \
  -a KO_Sup_Layer_neurons.txt \
  -b KD_Sup_Layer_neurons.txt \
  -o output_sup_layer_neurons.tsv \
  -n gene \
  -e logFC \
  -f FDR
```

### Alternative: Using Filtered Data

You can also use the filtered versions of the data:

```bash
cscore -i testdata/GSE236519 \
  -a filtered_kO_kD_Deep_Layer_neurons.txt \
  -b filtered_kO_kD_Interneurons.txt \
  -o output_filtered_comparison.tsv \
  -n gene \
  -e logFC \
  -f FDR
```

### What these examples show:
- How KO (knockout) and KD (knockdown) affect different neuronal cell types
- Which genes show consistent vs. divergent responses across different perturbation methods
- Cell-type-specific effects of genetic perturbations

---

## 3. Cross-Dataset Comparisons

You can also compare across different datasets or conditions:

### Example: Compare Bulk vs Single-cell Approaches

```bash
# First, prepare compatible data formats if needed, then run:
cscore -i testdata \
  -a GSE237099/GSE237099_1_unloading_reloading_Reloading_vs_Control_deseq2.txt \
  -b GSE236519/KO_Deep_Layer_neurons.txt \
  -o output_cross_dataset.tsv \
  -n gene \
  -e log2FoldChange \
  -f padj
```

Note: Cross-dataset comparisons require careful consideration of gene identifier compatibility and experimental context.

---

## 4. Working with Pre-computed Results

The repository also includes a pre-computed C-score result file:

```bash
# View pre-computed results
head testdata/GSE237099/GSE237099_1_unloading_reloading_cscore.txt
```

This file contains the output of running C-score analysis, including:
- Gene identifiers and metadata from both comparisons
- C-scores and significance values
- Convergence classifications

---

## 5. Advanced Options

### Adding Parallel Processing

```bash
cscore -i testdata/GSE237099 \
  -a GSE237099_1_unloading_reloading_Reloading_vs_Control_deseq2.txt \
  -b GSE237099_1_unloading_reloading_Unloading_vs_Control_deseq2.txt \
  -o output_parallel.tsv \
  -n ensembl_gene_id \
  -e log2FoldChange \
  -f padj \
  -w 4 \
  -s 42
```

### Using Custom Column Names

If your data has different column names, you can specify them:

```bash
cscore -i testdata/GSE236519 \
  -a KO_Deep_Layer_neurons.txt \
  -b KD_Deep_Layer_neurons.txt \
  -o output_custom_cols.tsv \
  -n gene \
  -e logFC \
  -f FDR
```

---

## 6. Data Format Requirements

### Input File Format

Your input TSV files should contain:
- **Gene identifier column**: First column by default, or specify with `-n`
- **Effect size column**: `avg_log2FC` by default, or specify with `-e`
- **Adjusted p-value column**: `p_val_adj` by default, or specify with `-f`

### Supported Column Name Variants

C-score automatically detects common column name variants:

**Effect size columns:**
- `avg_log2FC` (default for single-cell)
- `log2FoldChange` (common in bulk RNA-seq)
- `logFC`
- `avg_logFC`

**Adjusted p-value columns:**
- `p_val_adj` (default for single-cell)
- `padj` (common in bulk RNA-seq)
- `FDR`
- `q_value`

---

## 7. Simulation Data

The `simulation/` directory contains simulated C-score data for testing and validation:

```bash
# View simulation parameters and results
ls testdata/simulation/
# sim_hist.png - histogram of simulated scores
# sim_hist.tsv - histogram data
# sim_scores.tsv - raw simulated scores
```

These files can be used to understand the expected distribution of C-scores under different scenarios.

---

## Output Interpretation

### C-score Values
- **Positive scores**: Genes show consistent direction of change (both up or both down)
- **Negative scores**: Genes show divergent responses (up in one condition, down in another)
- **Magnitude**: Larger absolute values indicate stronger commonness/divergence

### Significance Testing
- **p-value**: Permutation-based significance of the observed C-score
- **q_value**: Multiple testing corrected p-value (Benjamini-Hochberg)
- **convergence**: Classification as "high" (common) or "low" (divergent) commonness

### Practical Interpretation
- Focus on genes with significant q-values (e.g., < 0.05)
- High convergence + positive score = consistent upregulation
- High convergence + negative score = consistent downregulation  
- Low convergence = divergent responses between conditions

---

## Getting Help

For additional options and help:

```bash
cscore -h
```

For more detailed documentation about the C-score methodology, see the [main README](../README.md).
