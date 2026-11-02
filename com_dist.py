#!/usr/bin/env python3
"""두 index 그룹의 질량중심 z 거리(최소영상)를 출력. 사용: com_dist.py conf.gro index.ndx POL PFA
v2 (2026-10-11, 강구 중령 수정 반영): 그룹 질량중심을 주기 경계를 고려한 원형 평균으로 계산한다.
 v1 은 단순 평균이라 슬래브가 상자 z 경계에 걸쳐 있으면 거리가 틀렸다 (PVDF: 2.125 vs 실제 2.983 nm) ->
 SMD 스텝 수 부족. ECTFE·PE·9월 PTFE 는 SMD 가 목표 거리에 도달했으므로 결과에는 영향 없음."""
import sys, math
gro, ndx, g1, g2 = sys.argv[1:5]
groups, cur = {}, None
for ln in open(ndx):
    s = ln.strip()
    if s.startswith('['): cur = s.strip('[] ').strip(); groups.setdefault(cur, []); continue
    if cur: groups[cur] += [int(x) for x in s.split()]
L = open(gro).read().splitlines(); n = int(L[1]); bz = float(L[2 + n].split()[2])
MASS = {'C': 12.011, 'F': 18.998, 'H': 1.008, 'O': 15.999, 'N': 14.007, 'S': 32.06, 'CL': 35.45, 'NA': 22.99}
def mass(name):
    e = ''.join(c for c in name if c.isalpha()).upper()
    return MASS.get(e[:2], MASS.get(e[:1], 12.0))
def comz(ids):
    # 원형 평균 (GROMACS pull 의 cosine-weighted COM 과 같은 취지): 상자 경계에 걸친 그룹도 바르게 처리
    m = cs = sn = 0.0
    for i in ids:
        l = L[1 + i]; w = mass(l[10:15].strip()); th = 2 * math.pi * float(l[36:44]) / bz
        m += w; cs += w * math.cos(th); sn += w * math.sin(th)
    return (math.atan2(sn / m, cs / m) / (2 * math.pi)) * bz
d = comz(groups[g2]) - comz(groups[g1]); d -= bz * round(d / bz)
print(f"{abs(d):.3f}")
