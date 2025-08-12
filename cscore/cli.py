from __future__ import annotations

import argparse
import os

from .core import CScoreConfig, run_cscore


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cscore", description="Compute C-score and permutation-based significance")
    parser.add_argument("-i", "--input_folder", required=True, help="path containing the comparison stats files")
    parser.add_argument("-a", "--comp1_file", required=True, help="first comparison file (TSV)")
    parser.add_argument("-b", "--comp2_file", required=True, help="second comparison file (TSV)")
    parser.add_argument("-o", "--output_file", required=True, help="output TSV path")
    parser.add_argument("-m", "--mode", choices=["gene", "pathway"], default="gene")
    parser.add_argument("-e", "--effect", default="avg_log2FC", help="effect size column; accepts common variants like log2FoldChange")
    parser.add_argument("-n", "--gname", default=None, help="key column name; if omitted or missing, uses the first column in each file and warns with the column names")
    parser.add_argument("-f", "--fdr", default="p_val_adj", help="FDR/p-adjusted column; accepts common variants like padj")
    parser.add_argument("-g", "--gtf", dest="gtf_file", default=None, help="gene annotation GTF for protein-coding annotation")
    parser.add_argument("-w", "--workers", type=int, default=None, help="number of workers (currently used for future parallelism)")
    parser.add_argument("-s", "--seed", type=int, default=1234, help="random seed for permutations")
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    default_workers = args.workers if args.workers is not None else max(1, (os.cpu_count() or 1))

    cfg = CScoreConfig(
        input_folder=args.input_folder,
        comp1_file=args.comp1_file,
        comp2_file=args.comp2_file,
        output_file=args.output_file,
        mode=args.mode,
        effect=args.effect,
        gname=args.gname,
        fdr=args.fdr,
        gtf_file=args.gtf_file,
        workers=default_workers,
        seed=args.seed,
    )

    run_cscore(cfg)


if __name__ == "__main__":
    main()


