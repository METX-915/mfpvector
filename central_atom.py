#!/usr/bin/env python3
"""index.ndx 의 한 그룹에서 공간적으로 가장 중심에 가까운 원자의 번호를 출력.

GROMACS 는 당김 기준 그룹이 상자 절반보다 크면 pbcatom 을 직접 지정하라고
요구한다. 자동 선택은 "원자 번호의 가운데"라서 공간 위치와 무관하다.

사용: central_atom.py <conf.gro> <index.ndx> <그룹명>
"""
import sys

gro, ndx, grp = sys.argv[1], sys.argv[2], sys.argv[3]

# --- index.ndx 에서 그룹 원자 번호(1-based) 읽기
idx, cur = [], None
for ln in open(ndx):
    s = ln.strip()
    if s.startswith('['):
        cur = s.strip('[] ').strip()
        continue
    if cur == grp:
        idx += [int(x) for x in s.split()]
if not idx:
    sys.exit(f"index.ndx 에 [ {grp} ] 그룹이 없습니다")

# --- .gro 좌표
L = open(gro).read().splitlines()
n = int(L[1].strip())
pos = {}
want = set(idx)
for k in range(n):
    if (k + 1) in want:
        l = L[2 + k]
        pos[k + 1] = (float(l[20:28]), float(l[28:36]), float(l[36:44]))
if len(pos) != len(idx):
    sys.exit(f"좌표를 다 찾지 못했습니다: {len(pos)}/{len(idx)}")

cx = sum(p[0] for p in pos.values()) / len(pos)
cy = sum(p[1] for p in pos.values()) / len(pos)
cz = sum(p[2] for p in pos.values()) / len(pos)
best = min(pos.items(), key=lambda kv: (kv[1][0]-cx)**2 + (kv[1][1]-cy)**2 + (kv[1][2]-cz)**2)
i, p = best
d = ((p[0]-cx)**2 + (p[1]-cy)**2 + (p[2]-cz)**2) ** 0.5
print(i)
print(f"# {grp}: {len(idx)} atoms, centre ({cx:.2f},{cy:.2f},{cz:.2f}) nm, "
      f"pbcatom {i} at ({p[0]:.2f},{p[1]:.2f},{p[2]:.2f}), 거리 {d:.2f} nm",
      file=sys.stderr)
