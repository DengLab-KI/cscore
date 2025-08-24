#!/usr/bin/env python3
import argparse, os, sys
import numpy as np
import pandas as pd
from cscore.core import compute_score

def simulate(n=10000, mu_fc1=0.0, mu_fc2=0.0, sd_fc=1.0, fdr_alpha=0.05, seed=1234):
    rng = np.random.default_rng(seed)
    fc1 = rng.normal(mu_fc1, sd_fc, size=n)
    fc2 = rng.normal(mu_fc2, sd_fc, size=n)
    mix = rng.uniform(0,1,size=n) < 0.2
    fdr1 = np.where(mix, rng.uniform(1e-6, fdr_alpha, size=n), rng.uniform(fdr_alpha, 1.0, size=n))
    fdr2 = np.where(mix, rng.uniform(1e-6, fdr_alpha, size=n), rng.uniform(fdr_alpha, 1.0, size=n))
    comp1 = np.column_stack([fc1, fdr1])
    comp2 = np.column_stack([fc2, fdr2])
    scores = compute_score(comp1, comp2)
    return scores

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=100000)
    ap.add_argument("--mu-fc1", type=float, default=0.0)
    ap.add_argument("--mu-fc2", type=float, default=0.0)
    ap.add_argument("--sd-fc", type=float, default=1.0)
    ap.add_argument("--fdr-alpha", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--out", type=str, default="/mnt/run/jh/my_github/cscore/misc/sim_scores.tsv")
    ap.add_argument("--hist", type=str, default="/mnt/run/jh/my_github/cscore/misc/sim_hist.tsv")
    ap.add_argument("--png", type=str, default="/mnt/run/jh/my_github/cscore/misc/sim_hist.png")
    args = ap.parse_args()
    scores = simulate(args.n, args.mu_fc1, args.mu_fc2, args.sd_fc, args.fdr_alpha, args.seed)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    pd.DataFrame({"score": scores}).to_csv(args.out, sep="\t", index=False)
    hist, bin_edges = np.histogram(scores, bins=100)
    centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    pd.DataFrame({"bin_center": centers, "count": hist}).to_csv(args.hist, sep="\t", index=False)
    print(f"Wrote {args.out} and {args.hist}")

    # Optional simple plot
    try:
        import matplotlib
        try:
            matplotlib.use("Agg", force=True)
        except Exception:
            pass
        import matplotlib.pyplot as plt
        os.makedirs(os.path.dirname(args.png), exist_ok=True)
        fig, ax = plt.subplots(figsize=(4, 3), dpi=150)
        ax.hist(scores, bins=100, color="#4C78A8", edgecolor="white", linewidth=0.3)
        # KDE overlay via Gaussian smoothing of histogram (Scott's rule bandwidth)
        try:
            bin_width = (bin_edges[1] - bin_edges[0]) if len(bin_edges) > 1 else 1.0
            n = float(scores.size)
            sigma = float(np.std(scores, ddof=1))
            iqr = float(np.subtract(*np.percentile(scores, [75, 25])) / 1.34)
            if sigma <= 0 and iqr > 0:
                sigma = iqr
            if sigma <= 0:
                sigma = 1.0
            h = 1.06 * sigma * (n ** (-1.0 / 5.0))
            h = max(h, 1e-12)
            sigma_bins = max(h / bin_width, 1e-6)
            L = int(np.ceil(6.0 * sigma_bins))
            L = max(L, 1)
            xk = np.arange(-L, L + 1, dtype=float)
            kernel = np.exp(-0.5 * (xk / sigma_bins) ** 2)
            s = kernel.sum()
            kernel = kernel / s if s > 0 else np.array([1.0])
            smoothed = np.convolve(hist.astype(float), kernel, mode="same")
            density = smoothed / (n * bin_width)
            ax.plot(centers, density, color="#F58518", linewidth=1.2)
        except Exception:
            pass
        ax.set_xlabel("C-score")
        ax.set_ylabel("Count")
        ax.set_title("C-score distribution")
        ax.grid(True, alpha=0.25, linewidth=0.3)
        fig.tight_layout()
        fig.savefig(args.png)
        plt.close(fig)
        print(f"Wrote {args.png}")
    except Exception as e:
        print(f"[cscore] Warning: could not create PNG plot ({e}). Install matplotlib to enable plotting.", file=sys.stderr)

if __name__ == "__main__":
    main()
