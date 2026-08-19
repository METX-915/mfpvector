"""
mfpvector RQ1 — Computational mass-balance bounding:
Can microfluoroplastic (MFP)-derived organofluorine be of the same order
as the measured EOF pool in river sediment?

DATA PROVENANCE (모든 입력 출처 명시; 허구 수치 없음):
[R1] EOF in Upper Ganges sediment = 180-780 ng F/g
     Müller et al. 2026, Environment International, DOI 10.1016/j.envint.2026.110104 (abstract, 실검증)
[R2] non-oxidised target PFAS explain <=1.3% of EOF (same paper) -> unexplained EOF ~ 98.7%+
[R3] total microplastics in river sediment ~ 140-27,259 particles/kg dry
     (literature range: Saigon 140-1200; Rhine 228-3763; Olifants/Sabie 2022-27259 items/kg) — WebSearch 2026-07
[R4] MFP = 2-8% of total microplastic COUNT in dust/sediment
     C&EN-reported LDIR study (2026-06) [FLAGGED: S2 미색인, DOI 수집불가, 보도기반 — 게재 전 CrossRef 확정 필요]
[R5] fluorine wt%: PTFE 76.0%, PVDF 59.4%, ECTFE ~39%, FEP 76.0% (화학량론 계산)
[R6] polymer density: PTFE 2.20, PVDF 1.78 g/cm3 (물성값); MFP-weighted ~2.0
[P1] effective particle size: log-uniform 20-300 um  [PRIOR: 퇴적물 MP 크기분포, 지배적 불확실성]
[P2] shape factor (fragment): 0.3  [PRIOR: 구형 대비 부피보정]

원칙: 절대 attribution 아닌 BOUNDING. 크기분포(d^3)가 지배적 불확실성 -> AL 우선측정 대상.
"""
import numpy as np
rng = np.random.default_rng(20260710)
N = 200_000

# --- reference organofluorine scale (measured) ---
EOF = rng.uniform(180, 780, N)                      # ng F / g dry  [R1]

# --- MFP particle number per g sediment ---
mp_per_kg = 10**rng.uniform(np.log10(140), np.log10(27259), N)  # particles/kg [R3] log-uniform
mp_per_g  = mp_per_kg / 1000.0
mfp_frac  = rng.uniform(0.02, 0.08, N)              # [R4]
mfp_per_g = mp_per_g * mfp_frac

# --- mass per MFP particle ---
d_um   = 10**rng.uniform(np.log10(20), np.log10(300), N)   # um [P1] dominant uncertainty
d_cm   = d_um * 1e-4
shape  = 0.3                                          # [P2]
vol_cm3 = shape * (np.pi/6.0) * d_cm**3               # cm^3
density = rng.uniform(1.78, 2.20, N)                 # g/cm3 [R6]
mass_particle_g = vol_cm3 * density                  # g per particle

# --- fluorine content ---
# aggregate F wt% across the 6 fluoropolymers, weighted uniform over plausible mix
fwt = rng.uniform(0.39, 0.76, N)                     # [R5] range ECTFE..PTFE/FEP

# --- MFP-derived fluorine concentration ---
mfp_mass_per_g = mfp_per_g * mass_particle_g          # g MFP / g sediment
mfp_F_ng_per_g = mfp_mass_per_g * fwt * 1e9           # ng F / g sediment

ratio = mfp_F_ng_per_g / EOF                          # MFP-F relative to measured EOF scale

def pct(x):
    q = np.percentile(x,[2.5,25,50,75,97.5])
    return dict(p2_5=q[0],p25=q[1],median=q[2],p75=q[3],p97_5=q[4])

print("=== RQ1 mass-balance bounding (N=%d draws) ==="%N)
print("MFP-derived fluorine [ng F / g dry sediment]:")
for k,v in pct(mfp_F_ng_per_g).items(): print(f"   {k:>7}: {v:.3g}")
print("\nRatio  MFP-F / measured EOF  [-]:")
for k,v in pct(ratio).items(): print(f"   {k:>7}: {v:.3g}")
print(f"\nP(MFP-F >= 1%% of EOF)   = {np.mean(ratio>=0.01):.3f}")
print(f"P(MFP-F >= 10%% of EOF)  = {np.mean(ratio>=0.10):.3f}")
print(f"P(MFP-F >= 100%% of EOF) = {np.mean(ratio>=1.0):.3f}")

# sensitivity: variance contribution via rank correlation (Spearman) of log-inputs to log-output
from scipy.stats import spearmanr
logout = np.log(mfp_F_ng_per_g)
inputs = {"total_MP(R3)":np.log(mp_per_kg),"MFP_frac(R4)":np.log(mfp_frac),
          "particle_size_d(P1)":np.log(d_um),"density(R6)":np.log(density),"F_wt%(R5)":np.log(fwt)}
print("\nSensitivity (|Spearman rho| of input vs MFP-F):")
for k,v in sorted(inputs.items(), key=lambda kv:-abs(spearmanr(kv[1],logout).statistic)):
    print(f"   {k:>20}: {spearmanr(v,logout).statistic:+.3f}")
