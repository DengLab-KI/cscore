# Examples

This page provides comprehensive examples of how to use `cscore` with simulations and real datasets.

## Simulated scenarios

**[→ Open Interactive Notebook](cscore_showcase.ipynb)**

The simulation demonstrates three fundamental scenarios:
1. **Independent responses**: No relationship
2. **Common responses**: Genes behave similarly in both comparisons
3. **Divergent responses**: Genes behave oppositely

---

## Real Data Examples

The following examples use datasets provided in the `testdata/` directory.

1. **GSE237099** - Bulk RNA-seq data (unloading/reloading muscle study)
2. **GSE236519** - Single-cell RNA-seq data (neuronal cell types)

---

## Bulk RNA-seq Example (GSE237099)

This example demonstrates C-score analysis using bulk RNA-seq data from a muscle unloading/reloading study.

```bash
cscore -i testdata/GSE237099 \
  -a GSE237099_1_unloading_reloading_Reloading_vs_Control_deseq2.txt \
  -b GSE237099_1_unloading_reloading_Unloading_vs_Control_deseq2.txt \
  -o output_bulk_example.tsv \
  -n ensembl_gene_id \
  -e log2FoldChange \
  -f padj \
  -w 16 ## parallel workers
```

### Expected Output:
The output file will contain columns including:
- `score`: C-score (positive = common direction; negative = divergent)
- `p`: permutation p-value for commonness/divergence
- `q_value`: Benjamini–Hochberg adjusted p-value
- `convergence`: "high" when commonness dominates, "low" otherwise
- All original columns from both input files with `_comp1`/`_comp2` suffixes

---

## 3. Single-cell RNA-seq Examples (GSE236519)

This dataset contains differential expression results for different neuronal cell types comparing knockout (KO) vs. knockdown (KD) conditions.

```bash
cscore -i testdata/GSE236519 \
  -a KO_Deep_Layer_neurons.txt \
  -b KD_Deep_Layer_neurons.txt \
  -o output_deep_layer_neurons.tsv \
  -n gene \
  -e logFC \
  -f FDR
  -w 16
```

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


---

## Getting Help

For additional options and help:

```bash
cscore -h
```
