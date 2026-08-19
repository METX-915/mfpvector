"""
Regenerates the adsorption figure (manuscript Fig 3, file Fig2.pdf) from data/.
Every plotted quantity is read from data/, with one stated exception noted below.

  (a) polymer-PFAS interaction-energy heatmap        <- data/experiment_summary.json
  (b) polymer-level mean IE with spread              <- data/experiment_summary.json
  (c) per-replicate IE distribution, zero-contact     <- data/raw_results.csv
      replicates marked explicitly
  (d) interaction energy vs umbrella-sampling free    <- data/experiment_summary.json
      energy for the two PMF-anchored systems            + data/pmf_result.json (PVDF)
                                                         + ECTFE dG = -28.5 kJ/mol, an author-
                                                         reported result (ESI Table S6); its
                                                         window-level output is not yet archived

Usage:  python code/make_figures.py
Output: figures/Fig2.pdf (vector, for submission) and figures/Fig2.png
"""
import json, os, csv, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures")
os.makedirs(FIG, exist_ok=True)
DATA = os.path.join(HERE, "..", "data")

plt.rcParams.update({
    "font.size": 13, "axes.spines.top": False, "axes.spines.right": False,
    "font.family": "DejaVu Sans", "axes.titleweight": "bold",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
BLUE, GREEN, ORANGE, GREY = "#0072B2", "#009E73", "#D55E00", "#7f7f7f"

summary = json.load(open(os.path.join(DATA, "experiment_summary.json")))
polymers = ["PTFE", "PVDF", "ECTFE", "PE"]
comps = ["PFOA", "PFBS", "GenX", "6:2 FTOH"]
key = {"PFOA": "PFOA", "PFBS": "PFBS", "GenX": "GenX", "6:2 FTOH": "FTOH62"}
IE = {p: [] for p in polymers}
for p in polymers:
    for c in comps:
        m = next(s for s in summary["systems"] if s["polymer"] == p and s["pfas"] == key[c])
        IE[p].append(m["meanIE"])

rep = collections.defaultdict(list)
for r in csv.DictReader(open(os.path.join(DATA, "raw_results.csv"))):
    rep[(r["polymer"], r["pfas"])].append(float(r["totalIE_kJmol"]))

fig, ax = plt.subplots(2, 2, figsize=(13, 10))
fig.subplots_adjust(hspace=0.42, wspace=0.40)

# (a) heatmap -- unchanged
M = np.array([IE[p] for p in polymers])
e = np.arange(-0.5, 4.5)
im = ax[0, 0].pcolormesh(e, e, M, cmap="viridis", shading="flat",
                         linewidth=0, edgecolors="face", rasterized=False)
ax[0, 0].set_aspect("auto"); ax[0, 0].invert_yaxis()
ax[0, 0].set_xticks(range(4)); ax[0, 0].set_xticklabels(comps, rotation=25, ha="right")
ax[0, 0].set_yticks(range(4)); ax[0, 0].set_yticklabels(polymers)
for i in range(4):
    for j in range(4):
        v = M[i, j]
        ax[0, 0].text(j, i, f"{v:.0f}", ha="center", va="center",
                      color="white" if v < -45 else "black", fontweight="bold")
cb = fig.colorbar(im, ax=ax[0, 0], fraction=0.046, pad=0.04); cb.set_label("IE (kJ/mol)")
ax[0, 0].set_title("(a) Polymer–PFAS interaction energy"); ax[0, 0].spines[:].set_visible(False)

# (b) polymer-level mean IE -- unchanged
means = [np.mean(IE[p]) for p in polymers]; sds = [np.std(IE[p]) for p in polymers]
ax[0, 1].bar(range(4), means, yerr=sds, capsize=5, color=[BLUE, BLUE, BLUE, GREY],
             edgecolor="black", linewidth=0.6)
ax[0, 1].axhline(0, color="grey", lw=0.8); ax[0, 1].set_xticks(range(4))
ax[0, 1].set_xticklabels(polymers)
ax[0, 1].set_ylabel("Mean IE (kJ/mol)"); ax[0, 1].set_title("(b) Adsorption strength by polymer")
ax[0, 1].set_ylim(-88, 2)
ax[0, 1].text(0.02, 0.03, "blue = fluoropolymer · grey = baseline",
              transform=ax[0, 1].transAxes, fontsize=11, color="grey")

# (c) per-replicate spread, zero-contact replicates marked
MK = {"PFOA": "o", "PFBS": "s", "GenX": "^", "6:2 FTOH": "D"}
rng = np.random.default_rng(11)
axc = ax[1, 0]
for i, p in enumerate(polymers):
    for k, c in enumerate(comps):
        v = np.array(rep[(p, key[c])])
        x = i + (k - 1.5) * 0.16 + rng.uniform(-0.03, 0.03, len(v))
        nz, z = v != 0.0, v == 0.0
        axc.scatter(x[nz], v[nz], marker=MK[c], s=42, facecolor=BLUE, edgecolor="black",
                    linewidth=0.5, alpha=0.85, zorder=3)
        if z.any():
            xz = x[z].copy()
            if len(xz) > 1:                      # separate co-located zero markers
                xz = xz.mean() + np.linspace(-0.055, 0.055, len(xz))
            axc.scatter(xz, v[z], marker=MK[c], s=95, facecolor="none", edgecolor=ORANGE,
                        linewidth=2.0, zorder=4)
    axc.plot([i - 0.34, i + 0.34], [means[i]] * 2, color="black", lw=2.4, zorder=5)
axc.axhline(0, color=ORANGE, lw=1.2, ls="--")
axc.set_xticks(range(4)); axc.set_xticklabels(polymers)
axc.set_ylabel("Replicate IE (kJ/mol)"); axc.set_xlim(-0.6, 3.6); axc.set_ylim(-95, 12)
axc.set_title("(c) Per-replicate spread (79 trajectories)")
axc.annotate("4 replicates beyond the 1.0 nm\ncutoff (IE = 0; fully desorbed)",
             xy=(2.35, 0), xytext=(1.05, -17), fontsize=12, color=ORANGE,
             arrowprops=dict(arrowstyle="->", color=ORANGE, lw=1.4))
axc.legend(handles=[Line2D([0], [0], marker=MK[c], ls="", mfc=BLUE, mec="black", ms=8, label=c)
                    for c in comps] +
                   [Line2D([0], [0], color="black", lw=2.4, label="polymer mean")],
           loc="lower left", fontsize=11.5, ncol=1, frameon=False, handletextpad=0.4, labelspacing=0.3)

# (d) interaction energy vs umbrella-sampling free energy
sysname = ["PVDF + PFOA", "ECTFE + PFOA"]
ie_anchor = [next(s["meanIE"] for s in summary["systems"] if s["polymer"] == "PVDF" and s["pfas"] == "PFOA"),
             next(s["meanIE"] for s in summary["systems"] if s["polymer"] == "ECTFE" and s["pfas"] == "PFOA")]
pmf = json.load(open(os.path.join(DATA, "pmf_result.json")))
assert pmf["system"] == "PVDF+PFOA"
dg = [pmf["delta_G_ads_kJ_per_mol"], -28.5]      # ECTFE value: ESI Table S6
x = np.arange(2); w = 0.34
axd = ax[1, 1]
b1 = axd.bar(x - w / 2, ie_anchor, w, color=BLUE, edgecolor="black", linewidth=0.6,
             label="Short-range IE")
b2 = axd.bar(x + w / 2, dg, w, color=GREEN, edgecolor="black", linewidth=0.6,
             label="PMF $\\Delta G_{\\mathrm{ads}}$")
for r in list(b1) + list(b2):
    axd.text(r.get_x() + r.get_width() / 2, r.get_height() + 2.4, f"{r.get_height():.1f}",
             ha="center", va="bottom", fontsize=12, fontweight="bold", color="white")
for i in range(2):
    axd.text(i, 5.4, f"|IE| = {abs(ie_anchor[i] / dg[i]):.1f} × |$\\Delta G$|",
             ha="center", va="center", fontsize=11.5, fontweight="bold", color=ORANGE)
axd.axhline(0, color="grey", lw=0.8)
axd.set_xticks(x); axd.set_xticklabels(sysname)
axd.set_ylabel("Energy (kJ/mol)"); axd.set_ylim(-88, 13)
axd.set_title("(d) Interaction energy is not a free energy")
axd.legend(loc="lower center", fontsize=12, frameon=False, ncol=2,
           bbox_to_anchor=(0.5, 0.02), columnspacing=2.2)

for ext, kw in (("pdf", {}), ("png", {"dpi": 400})):
    out = os.path.join(FIG, f"Fig2.{ext}")
    fig.savefig(out, bbox_inches="tight", **kw)
    print("wrote", os.path.normpath(out))
