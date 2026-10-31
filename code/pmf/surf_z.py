#!/usr/bin/env python3
"""슬래브의 Gibbs 분할면(고분자 질량밀도가 내부값의 1/2 이 되는 위치)을
POL 질량중심 기준 거리로 계산한다. PMF 의 x 축을 표면 기준으로 바꾸는 데 쓴다.

사용: surf_z.py <gro> <index.ndx>            (한 계)
      surf_z.py --all                        (네 계 + PE_TALL 를 한 번에)
출력 한 줄: 계이름  z_surface(nm)  내부밀도(kg/m^3)  계면폭 10-90%(nm)
"""
import sys, os, glob

MASS = {'C': 12.011, 'F': 18.998, 'H': 1.008, 'O': 15.999, 'N': 14.007,
        'S': 32.06, 'CL': 35.45, 'NA': 22.99}
BIN = 0.05           # nm

def surf(gro, ndx, name=""):
    cur, groups = None, {}
    for ln in open(ndx):
        s = ln.strip()
        if s.startswith('['):
            cur = s.strip('[] ').strip(); groups.setdefault(cur, []); continue
        if cur:
            groups[cur] += [int(x) for x in s.split()]
    pol = groups["POL"]
    L = open(gro).read().splitlines(); n = int(L[1])
    bx, by, bz = [float(x) for x in L[2 + n].split()[:3]]
    A = L[2:2 + n]

    def m(l):
        e = ''.join(c for c in l[10:15].strip() if c.isalpha()).upper()
        return MASS.get(e[:2], MASS.get(e[:1], 12.0))

    ws = [m(A[i - 1]) for i in pol]
    zs = [float(A[i - 1][36:44]) for i in pol]
    # 질량중심(주기경계 고려: 첫 원자를 기준으로 최소상 이미지로 모은 뒤 평균)
    z0 = zs[0]
    zu = [z - bz * round((z - z0) / bz) for z in zs]
    com = sum(w * z for w, z in zip(ws, zu)) / sum(ws)
    d = [z - com for z in zu]

    lo, hi = min(d), max(d)
    nb = int((hi - lo) / BIN) + 2
    h = [0.0] * nb
    for w, dz in zip(ws, d):
        h[int((dz - lo) / BIN)] += w
    vol = bx * by * BIN                                   # nm^3
    rho = [x / vol * 1.66054 * 1000 / 1000 for x in h]    # amu/nm^3 -> kg/m^3
    cen = [lo + (k + 0.5) * BIN for k in range(nb)]

    core = [r for c, r in zip(cen, rho) if abs(c) < 0.5]
    bulk = sum(core) / len(core)

    def cross(frac, side):
        t = frac * bulk
        rng = range(nb - 1, 0, -1) if side > 0 else range(0, nb - 1)
        prev = None
        for k in rng:
            if rho[k] >= t:
                if prev is None:
                    return cen[k]
                a, b = (cen[k], cen[prev]) if side > 0 else (cen[prev], cen[k])
                ra, rb = (rho[k], rho[prev]) if side > 0 else (rho[prev], rho[k])
                if rb == ra:
                    return b
                return a + (b - a) * (t - ra) / (rb - ra)
            prev = k
        return cen[-1] if side > 0 else cen[0]

    up, dn = cross(0.5, +1), cross(0.5, -1)
    zsurf = (up - dn) / 2.0                    # 반두께 = COM 에서 분할면까지
    width = abs(cross(0.10, +1) - cross(0.90, +1))
    print(f"{name or os.path.basename(os.path.dirname(gro)):10s} "
          f"z_surface = {zsurf:5.3f} nm   내부밀도 {bulk:6.1f} kg/m^3   "
          f"계면폭(10-90%) {width:4.2f} nm")
    return zsurf

if len(sys.argv) == 2 and sys.argv[1] == "--all":
    P = os.path.expanduser("~/mfpvector/pmf")
    for s in ["PTFE_PFOA", "PVDF_PFOA", "ECTFE_PFOA", "PE_PFOA", "PE_TALL"]:
        g, x = f"{P}/{s}/npt.gro", f"{P}/{s}/index.ndx"
        if os.path.exists(g) and os.path.exists(x):
            surf(g, x, s)
        else:
            print(f"{s:10s} (npt.gro 없음 — 건너뜀)")
else:
    surf(sys.argv[1], sys.argv[2])
