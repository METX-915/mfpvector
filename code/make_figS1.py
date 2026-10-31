"""Fig S1: PMF profiles (v2, 2026-10-20).
Usage: make_figS1.py figdir POLY=prof.xvg[:bs.xvg[:zlo:zhi]] ...   (order PTFE PVDF ECTFE PE)
prof.xvg = gmx wham -o profile; bs.xvg = gmx wham -bsres (z, mean, sd) for the +/-1 sigma band;
zlo:zhi = valid range after excluding sparse edge bins (report_pmf.py). Profiles are shifted so the
mean of the outer 25 % of the valid range (bulk plateau) is zero, as in report_pmf.py."""
import sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
figdir = sys.argv[1]
col = {"PTFE": "#D55E00", "PVDF": "#0072B2", "ECTFE": "#009E73", "PE": "#7f7f7f"}
def xvg(f): return np.array([[float(x) for x in ln.split()[:3]] for ln in open(f) if ln[0] not in "#@"])
fig, ax = plt.subplots(figsize=(6, 4))
for arg in sys.argv[2:]:
    p, spec = arg.split("=", 1); parts = spec.split(":")
    d = xvg(parts[0]); z, g = d[:, 0], d[:, 1]
    sd = None
    if len(parts) > 1 and parts[1]:
        b = xvg(parts[1]); sd = np.interp(z, b[:, 0], b[:, 2])
    if len(parts) > 3:
        m = (z >= float(parts[2])) & (z <= float(parts[3])); z, g = z[m], g[m]; sd = None if sd is None else sd[m]
    ref = g[z >= z.min() + 0.75 * (z.max() - z.min())].mean(); g = g - ref
    ax.plot(z, g, color=col[p], lw=1.8, label=f"{p}–PFOA")
    if sd is not None: ax.fill_between(z, g - sd, g + sd, color=col[p], alpha=0.18, lw=0)
ax.axhline(0, color="grey", lw=0.8, ls="--")
ax.set_xlabel("PFOA–slab centre-of-mass separation $z$ (nm)"); ax.set_ylabel("PMF (kJ mol$^{-1}$)")
ax.legend(frameon=False)
fig.savefig(f"{figdir}/FigS1.pdf", bbox_inches="tight"); fig.savefig(f"{figdir}/FigS1.png", dpi=200, bbox_inches="tight"); print("ok")
