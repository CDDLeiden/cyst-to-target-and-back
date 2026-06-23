"""Cell-death proxy for the ADPKD cyst-swelling experiments.

For both the target-validation and compound-exploration screens, the image
analysis reports, per well, the fraction of segmented nuclei that lack a
corresponding actin/cytoskeleton (organoid) signal (the Fraction_dead_cells
column). Apoptotic cells lose their cytoskeleton while their nuclear stain is
still detected, so this fraction is a relative proxy for cell death rather than
an absolute dead-cell count. A baseline of ~20% is present even in untreated
wells, because some segmented organoids are not captured during analysis and
their nuclei are then left without a matching object.

This script compares that proxy for each single-agent treatment (at the active
1 uM concentration) against the vehicle (DMSO) and FSK-stimulated controls and
against the staurosporine toxic control. It shows that the cyst-swelling
reduction of the validated compounds is not accompanied by signs of cell death,
and produces the supplementary figure for the target-validation experiment (the
cleanest toxic-control window). Summary numbers for both experiments are printed
to stdout.

Run from the repository root or from this directory:
    python cytotoxicity_analysis.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT_FIG = (
    ROOT
    / "figures"
    / "nuclei-based-toxicity-proxy.png"
)

# IBM colorblind-safe palette used throughout the paper's figures.
BLUE, MAGENTA, ORANGE, GREY = "#648FFF", "#DC267F", "#FE6100", "#9A9A9A"

DMSO_LIKE = {"DMSO"}

# Fraction of nuclei lacking a corresponding actin signal; used as a cell-death
# proxy. Plotted as a percentage.
PROXY_LABEL = "Nuclei without actin signal (%)"

# Compounds that reduced cyst swelling, by experiment (used only for colouring).
CS_REDUCERS = {
    "validation": {"Z4509024390", "Z211311146", "CPA", "Capadenoson",
                   "Aldosterone", "Esaxerenone"},
    "exploration": {"MIPS521", "Capadenoson", "27070328", "1249141",
                    "Esaxerenone"},
}

# Display order of single-agent treatments, grouped by target.
ORDER = {
    "validation": ["Z4509024390", "Z211311146", "A-804598", "JNJ-47965567",
                   "CPA", "Capadenoson", "DPCPX",
                   "Aldosterone", "Esaxerenone", "Finerenone"],
    "exploration": ["Capadenoson", "DPCPX", "MIPS521", "824745", "1237561",
                    "1249141", "1823372", "22755240", "27070328",
                    "Esaxerenone", "Finerenone", "Apararenone", "Benidipine",
                    "Z90308949", "Z318400112", "Z95680027"],
}

COLS = {
    "validation": ("treatment_type", "Treatment", "Treatment concentration",
                   "Treatment 2", "Treatment 2 concentration",
                   "Fraction_dead_cells",
                   "target_validation/"
                   "ADPKD-TargetValidationScreen_Batch3791_and_Batch3753.csv"),
    "exploration": ("plate.layout.info.Type", "plate.layout.info.Treatment",
                    "plate.layout.info.Treatment concentration",
                    "plate.layout.info.Treatment 2",
                    "plate.layout.info.Treatment 2 concentration",
                    "Fraction_dead_cells",
                    "compound_exploration/"
                    "ADPKD-CpdExplorationScreen_Batch4042_and_Batch4064.csv"),
}


def load(experiment):
    """Return a tidy frame with columns: kind, agent, proxy_pct.

    A single-agent well has exactly one real compound; the partner slot is
    either empty or DMSO. Across batches the active compound is stored either
    in the Treatment or in the Treatment 2 slot, so both are resolved here.
    proxy_pct is the percentage of segmented nuclei lacking a corresponding
    actin signal (cell-death proxy).
    """
    kind, t1, c1, t2, c2, proxy_col, rel = COLS[experiment]
    df = pd.read_csv(DATA / rel)
    df = df.rename(columns={kind: "kind", t1: "t1", c1: "c1", t2: "t2",
                            c2: "c2", proxy_col: "proxy"})
    df["proxy_pct"] = df["proxy"] * 100

    def agent_dose(row):
        a, b = row["t1"], row["t2"]
        a_real = pd.notna(a) and a not in DMSO_LIKE
        b_real = pd.notna(b) and b not in DMSO_LIKE
        if a_real and not b_real:
            return pd.Series({"agent": a, "dose": row["c1"]})
        if b_real and not a_real:
            return pd.Series({"agent": b, "dose": row["c2"]})
        return pd.Series({"agent": None, "dose": None})  # control or combination

    treat = df[df["kind"] == "treatment"].copy()
    treat[["agent", "dose"]] = treat.apply(agent_dose, axis=1)
    treat = treat.dropna(subset=["agent"])
    treat = treat[treat["dose"] == 1.0][["agent", "proxy_pct"]]
    treat["kind"] = "agent"

    ctrls = df[df["kind"].isin(["solvent_ctrl", "stim_only"])][["kind", "proxy_pct"]]
    ctrls = ctrls.assign(agent=ctrls["kind"])
    tox = df[(df["kind"] == "tox_ctrl") & (df["c1"] == 1.0)][["proxy_pct"]]
    tox = tox.assign(kind="tox_ctrl", agent="Staurosporine")
    return pd.concat([ctrls, treat, tox], ignore_index=True)


def conditions(tidy, experiment):
    """Yield (label, color, values) for each plotted condition, in order."""
    yield ("Vehicle\n(DMSO)", BLUE, tidy.loc[tidy["agent"] == "solvent_ctrl", "proxy_pct"])
    yield ("FSK\n(stim.)", BLUE, tidy.loc[tidy["agent"] == "stim_only", "proxy_pct"])
    for cpd in ORDER[experiment]:
        vals = tidy.loc[tidy["agent"] == cpd, "proxy_pct"]
        if vals.empty:
            continue
        color = MAGENTA if cpd in CS_REDUCERS[experiment] else GREY
        yield (cpd, color, vals)
    yield ("Stauro-\nsporine", ORANGE, tidy.loc[tidy["agent"] == "Staurosporine", "proxy_pct"])


def panel(ax, tidy, experiment):
    labels, colors, data = [], [], []
    for label, color, vals in conditions(tidy, experiment):
        labels.append(label)
        colors.append(color)
        data.append(vals.values)

    pos = list(range(len(data)))
    bp = ax.boxplot(data, positions=pos, widths=0.6, showfliers=False,
                    patch_artist=True, medianprops=dict(color="black"))
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
    for p, vals, color in zip(pos, data, colors):
        ax.scatter([p] * len(vals), vals, s=16, color=color,
                   edgecolor="black", linewidth=0.3, zorder=3)

    ax.axhline(float(pd.Series(data[0]).median()), ls="--", lw=0.8,
               color="grey", zorder=0)
    ax.set_xticks(pos)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(PROXY_LABEL)
    ax.set_ylim(0, 100)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor=BLUE, alpha=0.55, label="Control"),
        Patch(facecolor=MAGENTA, alpha=0.55, label="CS-reducing compound"),
        Patch(facecolor=GREY, alpha=0.55, label="Inactive / CS-enhancing compound"),
        Patch(facecolor=ORANGE, alpha=0.55, label="Staurosporine (toxic control)"),
    ], fontsize=7, loc="upper left", ncol=2, frameon=False)


def summary(tidy, experiment):
    print(f"\n=== {experiment.upper()} — median nuclei without actin signal (%) ===")
    for label, _, vals in conditions(tidy, experiment):
        print(f"  {label.replace(chr(10), ' '):<18} n={len(vals):<3} median={vals.median():5.1f}")


def main():
    val, exp = load("validation"), load("exploration")
    summary(val, "validation")
    summary(exp, "exploration")

    fig, ax = plt.subplots(figsize=(6, 4))
    panel(ax, val, "validation")
    fig.tight_layout()
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
    print(f"\nSaved figure to {OUT_FIG}")


if __name__ == "__main__":
    main()
