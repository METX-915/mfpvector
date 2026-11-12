"""Surface partition coefficient K_surf = int [exp(-W/kT) - 1] dz from each PFOA PMF (data/pmf/*_prof.xvg),
and particle-water K_d = K_surf * a_s for smooth spheres a_s = 6/(rho d). Values quoted in SI Table S16.
Usage: python code/ksurf.py data/pmf"""
import sys, numpy as np
d0 = sys.argv[1] if len(sys.argv) > 1 else "data/pmf"; kT = 2.494  # kJ/mol at 300 K
for p in ["PTFE", "PVDF", "ECTFE", "PE"]:
    a = np.array([[float(x) for x in l.split()[:2]] for l in open(f"{d0}/{p}_prof.xvg") if l[0] not in "#@" and l.strip()])
    z, W = a[:, 0], a[:, 1]; W = W - W[-len(W) // 4:].mean()
    K = np.trapezoid(np.exp(-W / kT) - 1, z) * 1e-9  # nm -> m
    print(f"{p}: Wmin {W.min():.1f} kJ/mol, K_surf {K:.2e} m, K_d(20-300 um, rho 2 g/cm3) {K*10*1e3:.2f}-{K*150*1e3:.1f} L/kg")
