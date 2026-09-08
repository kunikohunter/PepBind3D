#!/usr/bin/env python3
"""Paired comparisons across the six benchmark arms on their shared targets.

Every arm predicts the same targets, so comparisons are PAIRED: each target
contributes one difference. Unpaired tests would discard that pairing and
inflate variance with between-target difficulty, which dominates the spread.

Tests
-----
selected / oracle  Wilcoxon signed-rank (two-sided, exact where n allows).
                   Reported with the Hodges-Lehmann estimate -- the median of
                   the Walsh averages of the paired differences, which is the
                   location estimate Wilcoxon actually tests -- and matched-
                   pairs rank-biserial correlation as effect size.
<=2A               McNemar exact test on the discordant pairs. A paired binary
                   outcome needs McNemar, not Wilcoxon.

Multiplicity
------------
15 arm pairs per metric. Holm-Bonferroni within each metric family; Holm is
uniformly more powerful than Bonferroni and still controls FWER without
assuming independence, which these correlated arms plainly violate.

Usage:  paired_tests.py [--metric selected|oracle|sub2|all]
"""
import csv, itertools, sys
from pathlib import Path
import numpy as np
from scipy import stats

BENCH = Path("/data/p_csb_meiler/huntek1/benchmark")
COFOLD = BENCH / "cofold_arm"
COFOLD_ARMS = ["af2_reuse", "boltz1", "boltz2", "chai1", "protenix"]
SUB2_THRESHOLD = 2.0


def load_arms():
    arms = {}
    for a in COFOLD_ARMS:
        p = COFOLD / a / "results" / f"{a}_summary.csv"
        if p.exists():
            arms[a] = {r["pdb_id"]: r for r in csv.DictReader(open(p))}
    ros = BENCH / "rosetta_arm" / "results" / "arm_summary.csv"
    if ros.exists():
        arms["rosetta"] = {r["pdb_id"]: r for r in csv.DictReader(open(ros))}
    if not arms:
        sys.exit("no arm summaries found")
    return arms


def hodges_lehmann(d):
    """Median of Walsh averages (d_i + d_j)/2 over i <= j."""
    n = len(d)
    walsh = [(d[i] + d[j]) / 2.0 for i in range(n) for j in range(i, n)]
    return float(np.median(walsh))


def rank_biserial(d):
    """Matched-pairs rank-biserial: (W+ - W-) / (W+ + W-), zeros dropped."""
    d = d[d != 0]
    if len(d) == 0:
        return 0.0
    r = stats.rankdata(np.abs(d))
    wp, wn = r[d > 0].sum(), r[d < 0].sum()
    return float((wp - wn) / (wp + wn))


def holm(pvals):
    m = len(pvals)
    order = np.argsort(pvals)
    adj = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(order):
        val = (m - rank) * pvals[idx]
        running = max(running, val)
        adj[idx] = min(running, 1.0)
    return adj


def run_continuous(arms, common, metric):
    vals = {a: np.array([float(arms[a][t][metric]) for t in common]) for a in arms}
    rows = []
    for x, y in itertools.combinations(sorted(arms), 2):
        d = vals[x] - vals[y]            # negative => x better (lower RMSD)
        if np.all(d == 0):
            continue
        stat, p = stats.wilcoxon(d, alternative="two-sided", zero_method="wilcox")
        rows.append({
            "a": x, "b": y, "n": len(d),
            "median_diff": float(np.median(d)),
            "hl": hodges_lehmann(d),
            "rbc": rank_biserial(d),
            "p": p,
        })
    for r, padj in zip(rows, holm([r["p"] for r in rows])):
        r["p_holm"] = padj
    return rows


def run_sub2(arms, common):
    ok = {a: np.array([float(arms[a][t]["selected"]) <= SUB2_THRESHOLD for t in common])
          for a in arms}
    rows = []
    for x, y in itertools.combinations(sorted(arms), 2):
        b = int(np.sum(ok[x] & ~ok[y]))   # x succeeds, y fails
        c = int(np.sum(~ok[x] & ok[y]))
        if b + c == 0:
            p = 1.0
        else:
            p = float(stats.binomtest(b, b + c, 0.5).pvalue)
        rows.append({"a": x, "b_arm": y, "b": b, "c": c,
                     "rate_a": float(ok[x].mean()), "rate_b": float(ok[y].mean()), "p": p})
    for r, padj in zip(rows, holm([r["p"] for r in rows])):
        r["p_holm"] = padj
    return rows


def run_ensemble_value(arms, common):
    """Within-arm: is the decoy the model's own confidence picks worse than the
    best decoy in its 25-member ensemble? One-sided Wilcoxon (selected > oracle).

    This is the direct test of whether generating 25 samples adds value over a
    single prediction. It is a WITHIN-arm comparison, so it carries no
    multiplicity correction across arms -- each arm is its own question.
    """
    print("=== ensemble value: selected vs oracle within each arm ===")
    print(f"{'arm':<11}{'selected':>10}{'oracle':>9}{'gap':>8}{'rank_sel':>10}{'p':>12}")
    for a in sorted(arms):
        sel = np.array([float(arms[a][t]["selected"]) for t in common])
        orc = np.array([float(arms[a][t]["oracle"]) for t in common])
        _, p = stats.wilcoxon(sel - orc, alternative="greater")
        rk = [arms[a][t].get("rank_of_selected", "") for t in common]
        rk = [float(x) for x in rk if x not in ("", "None")]
        mr = np.median(rk) if rk else float("nan")
        print(f"{a:<11}{np.median(sel):>10.3f}{np.median(orc):>9.3f}"
              f"{np.median(sel-orc):>8.3f}{mr:>10.1f}{p:>12.2e}")
    print("\ngap = median(selected - oracle); rank_sel = median rank of the "
          "selected decoy (1 = confidence picked the best of 25)\n")


def main():
    metric = "all"
    if len(sys.argv) > 1 and sys.argv[1] == "--metric":
        metric = sys.argv[2]
    arms = load_arms()
    common = sorted(set.intersection(*[set(d) for d in arms.values()]))
    print(f"arms: {', '.join(f'{a} (n={len(d)})' for a, d in sorted(arms.items()))}")
    print(f"shared targets: {len(common)}\n")

    for m in (["selected", "oracle"] if metric in ("all",) else
              ([metric] if metric in ("selected", "oracle") else [])):
        print(f"=== {m} RMSD -- Wilcoxon signed-rank, Holm-corrected over 15 pairs ===")
        print(f"{'comparison':<24}{'med_diff':>10}{'HL':>9}{'rbc':>8}{'p':>12}{'p_holm':>12}  better")
        for r in sorted(run_continuous(arms, common, m), key=lambda r: r["p"]):
            better = r["a"] if r["hl"] < 0 else r["b"]
            sig = "" if r["p_holm"] >= 0.05 else ("  *" if r["p_holm"] >= 0.01 else "  **")
            print(f"{r['a']+' vs '+r['b']:<24}{r['median_diff']:>10.3f}{r['hl']:>9.3f}"
                  f"{r['rbc']:>8.3f}{r['p']:>12.2e}{r['p_holm']:>12.2e}  {better}{sig}")
        print()

    if metric in ("all", "ensemble"):
        run_ensemble_value(arms, common)

    if metric in ("all", "sub2"):
        print(f"=== fraction of targets with selected <= {SUB2_THRESHOLD} A -- McNemar exact ===")
        print(f"{'comparison':<24}{'rate_a':>8}{'rate_b':>8}{'b':>5}{'c':>5}{'p':>12}{'p_holm':>12}  better")
        for r in sorted(run_sub2(arms, common), key=lambda r: r["p"]):
            better = r["a"] if r["b"] > r["c"] else r["b_arm"]
            sig = "" if r["p_holm"] >= 0.05 else ("  *" if r["p_holm"] >= 0.01 else "  **")
            print(f"{r['a']+' vs '+r['b_arm']:<24}{r['rate_a']:>8.3f}{r['rate_b']:>8.3f}"
                  f"{r['b']:>5}{r['c']:>5}{r['p']:>12.2e}{r['p_holm']:>12.2e}  {better}{sig}")


if __name__ == "__main__":
    main()
