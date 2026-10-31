#!/usr/bin/env python3
"""PMF 판정 v6 (2026-09-17).
v5 대비 수정:
 - 창 중첩을 파일 순서가 아니라 '좌표 순'으로 계산한다.
   (보간 창을 목록 끝에 덧붙이면 v5 는 3.89 와 2.52 를 인접으로 보고 0% 를 냈다.)
 - 중첩 부족 지점마다 두 창의 좌표와 그 중간값을 찍는다 -> 어디에 창을 넣을지 바로 나온다.
 - 중첩 부족이 있으면 '판정: 사용 가능' 을 내지 않는다.
사용: report_pmf.py <계이름> profile.xvg bsResult.xvg histo.xvg
"""
import sys, json
import numpy as np
def read_xvg(path):
    return np.array([[float(x) for x in l.split()] for l in open(path) if l[0] not in "#@" and l.strip()])
name, fp, fb, fh = sys.argv[1:5]
prof = read_xvg(fp); z, g = prof[:, 0], prof[:, 1]
H = read_xvg(fh); tot = H[:, 1:].sum(1)
good = np.isfinite(g) & (tot >= 0.10 * tot.max())
zt = z.min() + 0.75 * (z.max() - z.min()); tail = good & (z >= zt)
g0 = g - g[tail].mean()
i = int(np.argmin(np.where(good, g0, np.inf))); dG, zmin = float(g0[i]), float(z[i])
err = float("nan")
try:
    bs = read_xvg(fb)
    if bs.shape[1] >= 3: err = float(np.sqrt(bs[tail, 2].mean() ** 2 + bs[i, 2] ** 2))
except Exception: pass
plat_sd = float(g0[tail].std()); slope = float(np.polyfit(z[tail], g0[tail], 1)[0])
zg = z[good]
print(f"{name}")
print(f"  dG_ads = {dG:7.2f} +/- {err:.2f} kJ/mol   (최소점 z = {zmin:.3f} nm)")
print(f"  평탄부: z >= {zt:.2f} nm, 표준편차 {plat_sd:.2f}, 기울기 {slope:+.1f} kJ/mol/nm")
print(f"  유효 구간: {zg.min():.3f} – {zg.max():.3f} nm (표본 희박 구간 {(~good).sum()}개 제외)")
flags = []
if abs(zmin - zg.min()) < 0.05: flags.append("최소점이 안쪽 가장자리 — 안쪽 창 추가 권장")
if abs(slope) > 5: flags.append("평탄부 기울기 큼 — 바깥쪽 창 추가 권장")
ovmin = float("nan")
try:
    Hm = H[:, 1:]; hz = H[:, 0]; n = Hm.shape[1]
    csum = np.maximum(Hm.sum(0), 1e-12)
    cen = (Hm * hz[:, None]).sum(0) / csum          # 창별 히스토그램 무게중심
    order = np.argsort(cen)                          # 좌표 순 정렬 (핵심 수정)
    Hs = Hm[:, order]; cs = cen[order]
    ov, bad = [], []
    for k in range(n - 1):
        a, b = Hs[:, k], Hs[:, k + 1]
        frac = np.minimum(a, b).sum() / min(a.sum(), b.sum())
        ov.append(frac)
        if frac < 0.10: bad.append(k)
    ov = np.array(ov); ovmin = float(ov.min())
    print(f"  인접 창 중첩(좌표순): 평균 {ov.mean()*100:.0f}%, 최소 {ov.min()*100:.0f}%")
    gaps = np.diff(cs)
    print(f"  창 간격: 최소 {gaps.min():.3f} / 최대 {gaps.max():.3f} nm, 창 {n}개")
    if bad:
        print(f"  *** 중첩 부족 {len(bad)} 곳 — 각 구간 중간에 창 추가 필요")
        for k in bad:
            print(f"        w{order[k]:02d}({cs[k]:.3f}) – w{order[k+1]:02d}({cs[k+1]:.3f})"
                  f" : {ov[k]*100:.0f}%   -> 추가 지점 {0.5*(cs[k]+cs[k+1]):.3f} nm")
        flags.append(f"중첩 부족 {len(bad)} 곳 — 창 보간 필요")
    else:
        print("  중첩 양호 (모든 인접 쌍 10% 이상)")
except Exception as e:
    print(f"  히스토그램 확인 실패: {e}")
if not np.isnan(err) and abs(err) > 0.15 * abs(dG):
    flags.append("bootstrap 오차가 dG 의 15% 초과 — 전 창 연장 권장")
for f in flags: print("  *** " + f)
print("  판정: " + ("사용 가능" if not flags else "보완 필요"))
json.dump(dict(system=name, dG=round(dG, 2), err=round(err, 2), zmin=zmin,
               plateau_from=round(zt, 2), plateau_sd=round(plat_sd, 2),
               overlap_min=round(ovmin, 3)), open(f"{name}_pmf.json", "w"), indent=1)
