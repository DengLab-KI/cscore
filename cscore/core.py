from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Iterable, Tuple, Optional

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.utils import shuffle as sk_shuffle
import sys


def weight(fdr_values: np.ndarray) -> np.ndarray:
    fdr_values = np.asarray(fdr_values, dtype=float)
    clipped = np.clip(fdr_values, 1e-300, 1.0)
    weights = np.where(clipped < 0.05, 1.0, np.log10(clipped) / np.log10(0.05))
    return weights


def ratio(fc_comp1: np.ndarray, fc_comp2: np.ndarray) -> np.ndarray:
    fc1 = np.asarray(fc_comp1, dtype=float)
    fc2 = np.asarray(fc_comp2, dtype=float)
    same_direction = fc1 * fc2 > 0
    max_abs = np.maximum(np.abs(fc1), np.abs(fc2))
    diff_abs = np.abs(fc1 - fc2)
    pos_part = max_abs / (diff_abs + 1.0)
    neg_part = -diff_abs / (max_abs + 1.0)
    return np.where(same_direction, pos_part, neg_part)


def compute_score(comp1_np: np.ndarray, comp2_np: np.ndarray) -> np.ndarray:
    fc1 = comp1_np[:, 0]
    fc2 = comp2_np[:, 0]
    w1 = weight(comp1_np[:, 1])
    w2 = weight(comp2_np[:, 1])
    magnitude = np.abs(fc1 * w1) + np.abs(fc2 * w2)
    return magnitude * ratio(fc1, fc2)


def _permute_rows_sklearn(array: np.ndarray, random_state: int) -> np.ndarray:
    return sk_shuffle(array, random_state=random_state)


def compute_pvalues(
    comp1_np: np.ndarray,
    comp2_np: np.ndarray,
    observed_scores: np.ndarray,
    num_permutations: int,
    seed: Optional[int] = None,
    workers: int = 1,
) -> Tuple[np.ndarray, list[str]]:
    """Compute permutation p-values.

    Uses scikit-learn shuffling. If workers > 1, runs batched permutations in parallel
    and aggregates tail counts to avoid storing full permutation matrices.
    """
    num_genes = observed_scores.shape[0]

    def _batch_count(batch_index: int, batch_size: int) -> np.ndarray:
        start_seed = (seed or 0) + batch_index * batch_size
        counts = np.zeros(num_genes, dtype=np.int64)
        for j in range(batch_size):
            step_seed = start_seed + j
            perm1 = _permute_rows_sklearn(comp1_np, random_state=step_seed)
            perm2 = _permute_rows_sklearn(comp2_np, random_state=step_seed + 40000)
            perm_scores = compute_score(perm1, perm2)
            counts += (perm_scores > observed_scores).astype(np.int64)
        return counts

    if workers is None or workers < 1:
        workers = 1

    if workers == 1:
        greater_counts = _batch_count(0, num_permutations)
    else:
        # Choose a batch size that balances overhead and memory
        batch_size = 512 if num_permutations >= 4096 else max(1, num_permutations // (workers * 2))
        num_batches, rem = divmod(num_permutations, batch_size)
        batches = [batch_size] * num_batches + ([rem] if rem else [])
        counts_list = Parallel(n_jobs=workers, prefer="threads")(
            delayed(_batch_count)(i, bsz) for i, bsz in enumerate(batches)
        )
        greater_counts = np.sum(counts_list, axis=0, dtype=np.int64)

    ps_tmp = greater_counts.astype(float) / float(num_permutations)
    sense = ["high" if p < 0.5 else "low" for p in ps_tmp]
    pvals = np.where(ps_tmp < 0.5, ps_tmp, 1.0 - ps_tmp)
    return pvals, sense


def _bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = q
    return np.clip(out, 0.0, 1.0)


@dataclass
class CScoreConfig:
    input_folder: str
    comp1_file: str
    comp2_file: str
    output_file: str
    mode: str = "gene"
    effect: Optional[str] = "avg_log2FC"
    gname: Optional[str] = None
    fdr: Optional[str] = "p_val_adj"
    gtf_file: Optional[str] = None
    workers: int = max(1, (os.cpu_count() or 1))
    seed: Optional[int] = 1234


def run_cscore(cfg: CScoreConfig) -> None:
    comp1_path = os.path.join(cfg.input_folder, cfg.comp1_file)
    comp2_path = os.path.join(cfg.input_folder, cfg.comp2_file)

    comp1 = pd.read_csv(comp1_path, sep="\t", decimal=".").dropna()
    comp2 = pd.read_csv(comp2_path, sep="\t", decimal=".").dropna()

    # Normalize column names: strip whitespace and BOMs
    def _normalize_cols(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        df.columns = [str(c).replace("\ufeff", "").strip() for c in df.columns]
        return df

    comp1 = _normalize_cols(comp1)
    comp2 = _normalize_cols(comp2)

    # Determine key columns
    comp1_key = cfg.gname if (cfg.gname in comp1.columns if cfg.gname else False) else comp1.columns[0]
    comp2_key = cfg.gname if (cfg.gname in comp2.columns if cfg.gname else False) else comp2.columns[0]

    if (cfg.gname is None) or (comp1_key != (cfg.gname or comp1_key)) or (comp2_key != (cfg.gname or comp2_key)):
        print(
            f"[cscore] Using key column '{comp1_key}' for comp1 and '{comp2_key}' for comp2.",
            file=sys.stderr,
        )

    merge_key = cfg.gname if (cfg.gname and comp1_key == cfg.gname and comp2_key == cfg.gname) else "__gene_key__"
    if merge_key not in comp1.columns:
        comp1[merge_key] = comp1[comp1_key]
    if merge_key not in comp2.columns:
        comp2[merge_key] = comp2[comp2_key]

    rows_keep = np.intersect1d(comp1[merge_key], comp2[merge_key])
    comp1 = comp1[comp1[merge_key].isin(rows_keep)].sort_values(merge_key).reset_index(drop=True)
    comp2 = comp2[comp2[merge_key].isin(rows_keep)].sort_values(merge_key).reset_index(drop=True)

    def _resolve_column(df: pd.DataFrame, requested: Optional[str], kind: str) -> str:
        """Resolve column by rules:
        - If user specified a name, use it (case-insensitive, trimmed); error if absent.
        - If not specified, auto-detect from common names; otherwise fall back heuristics.
        """
        def _norm(s: str) -> str:
            return s.replace("\ufeff", "").strip().lower()

        lower_to_actual = {_norm(c): c for c in df.columns}

        if requested is not None:
            key = _norm(requested)
            if key in lower_to_actual:
                actual = lower_to_actual[key]
                return actual
            raise KeyError(
                f"Column '{requested}' not found in file. Available: {list(df.columns)}"
            )

        # Auto-detect when not requested
        if kind == "effect":
            candidates = ["avg_log2FC", "log2FoldChange", "logFC", "avg_logFC", "log2fc", "l2fc"]
        else:
            candidates = ["p_val_adj", "padj", "p_adj", "FDR", "qval", "q_value", "qvalue", "pvalue", "p_val", "p"]

        for cand in candidates:
            key = _norm(cand)
            if key in lower_to_actual:
                return lower_to_actual[key]

        # Heuristic fallback
        if kind == "effect":
            for col in df.columns:
                if col == merge_key:
                    continue
                if pd.api.types.is_numeric_dtype(df[col]):
                    print(f"[cscore] Falling back to effect column '{col}'.", file=sys.stderr)
                    return col
        else:
            ordered = [c for c in df.columns if _norm(c) in {"p_val_adj", "padj", "p_adj", "fdr", "q_value", "qval", "qvalue"}]
            if ordered:
                col = ordered[0]
                print(f"[cscore] Falling back to FDR column '{col}'.", file=sys.stderr)
                return col
            for col in df.columns:
                if _norm(col) in {"p", "pvalue", "p_val"}:
                    print(f"[cscore] Falling back to FDR column '{col}'.", file=sys.stderr)
                    return col
        raise KeyError(f"Could not resolve a valid column for {kind} in file.")

    comp1_effect = _resolve_column(comp1, cfg.effect, kind="effect")
    comp1_fdr = _resolve_column(comp1, cfg.fdr, kind="fdr")
    comp2_effect = _resolve_column(comp2, cfg.effect, kind="effect")
    comp2_fdr = _resolve_column(comp2, cfg.fdr, kind="fdr")

    print(
        f"[cscore] comp1 columns -> key='{merge_key}', effect='{comp1_effect}', fdr='{comp1_fdr}'",
        file=sys.stderr,
    )
    print(
        f"[cscore] comp2 columns -> key='{merge_key}', effect='{comp2_effect}', fdr='{comp2_fdr}'",
        file=sys.stderr,
    )

    comp1_np = comp1[[comp1_effect, comp1_fdr]].to_numpy(dtype=float)
    comp2_np = comp2[[comp2_effect, comp2_fdr]].to_numpy(dtype=float)

    scores = compute_score(comp1_np, comp2_np)
    keep_mask = scores != 0
    comp1 = comp1.loc[keep_mask]
    comp2 = comp2.loc[keep_mask]
    comp1_np = comp1_np[keep_mask]
    comp2_np = comp2_np[keep_mask]
    scores = scores[keep_mask]

    # permutations:  k = n^2 if n < 200 else 40000
    n = comp1_np.shape[0]
    num_permutations = n * n if n < 200 else 40000
    pvals, sense = compute_pvalues(
        comp1_np, comp2_np, scores, num_permutations=num_permutations, seed=cfg.seed, workers=cfg.workers
    )

    merged = pd.merge(comp1, comp2, on=merge_key, suffixes=("_comp1", "_comp2"))
    merged = merged.sort_values(merge_key).reset_index(drop=True)
    merged["score"] = scores
    merged["p"] = pvals
    merged["convergence"] = sense

    if cfg.mode == "gene" and cfg.gtf_file:
        gencode = pd.read_table(
            cfg.gtf_file,
            comment="#",
            sep="\t",
            names=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"],
        )
        gencode_genes = gencode[(gencode.feature == "transcript")][["seqname", "start", "end", "attribute"]].copy().reset_index().drop(
            "index", axis=1
        )

        def gene_info(attributes: str) -> Tuple[str, str]:
            fields = attributes.split(";")
            gene_name = [f for f in fields if "gene_name" in f][0].split(" ")[2].strip('"')
            gene_type = [f for f in fields if "gene_type" in f][0].split(" ")[2]
            return gene_name, gene_type

        gencode_genes["gene_name"], gencode_genes["gene_type"] = zip(*gencode_genes.attribute.apply(gene_info))
        pc_genes = gencode_genes.query("gene_type=='\"protein_coding\"'")
        pc_gene_set = set(pc_genes.gene_name)
        merged["coding"] = merged[merge_key].isin(pc_gene_set)

    # BH correction on p-values
    merged["q_value"] = _bh_fdr(merged["p"].to_numpy())

    merged.to_csv(cfg.output_file, index=False, float_format="%.64f", sep="\t")


