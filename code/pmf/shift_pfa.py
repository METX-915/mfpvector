#!/usr/bin/env python3
"""PFAS 를 z 로 평행이동해 슬래브 COM 에서 목표 거리에 두고, 겹치는 물 분자를 제거한다.  (v2)

v2 수정 사항
 1) [molecules] 에 SOL 블록이 여러 개인 경우(genion 이 이온을 끼워 넣으면 쪼개진다)
    블록마다 제거 개수를 따로 세어 각각 갱신한다.
    v1 은 첫 블록만 줄여서 top 과 gro 의 원자 순서가 어긋났다.
 2) PFAS 를 원자별이 아니라 분자 통째로 감싼다. 원자별 % bz 는 상자 경계에서 분자를 찢는다.
 3) 끝에 자체 검증: 원자 수, 제거 개수, 갱신된 [molecules] 를 출력한다.

사용: shift_pfa.py in.gro index.ndx system.top TARGET out.gro out.top out.ndx
 - 방향은 현재 부호를 유지(위쪽에 있으면 위쪽에서 목표 거리).
"""
import sys

gro, ndx, top, target, ogro, otop, ondx = sys.argv[1:8]
target = float(target)
SOLNAMES = ("SOL", "WAT", "HOH")

# ---- index ----
groups, cur = {}, None
for ln in open(ndx):
    s = ln.strip()
    if s.startswith('['):
        cur = s.strip('[] ').strip(); groups.setdefault(cur, []); continue
    if cur:
        groups[cur] += [int(x) for x in s.split()]
pol, pfa = groups["POL"], groups["PFA"]

# ---- gro ----
L = open(gro).read().splitlines()
n = int(L[1])
box = L[2 + n]
bx, by, bz = [float(x) for x in box.split()[:3]]
atoms = L[2:2 + n]

MASS = {'C': 12.011, 'F': 18.998, 'H': 1.008, 'O': 15.999, 'N': 14.007,
        'S': 32.06, 'CL': 35.45, 'NA': 22.99}
def mass(l):
    e = ''.join(c for c in l[10:15].strip() if c.isalpha()).upper()
    return MASS.get(e[:2], MASS.get(e[:1], 12.0))
def zc(i): return float(atoms[i - 1][36:44])
def comz(ids):
    m = s = 0.0
    for i in ids:
        w = mass(atoms[i - 1]); m += w; s += w * zc(i)
    return s / m

d = comz(pfa) - comz(pol)
d -= bz * round(d / bz)
sgn = 1 if d >= 0 else -1
shift = sgn * target - d
# 분자를 찢지 않도록 COM 기준으로 통째로 감싼다
newcom = comz(pfa) + shift
shift -= bz * (newcom // bz)
print(f"    현재 부호있는 거리 {d:+.3f} nm -> 목표 {sgn*target:+.3f}, 이동 {shift:+.3f} nm")

pset = set(pfa)
newatoms = []
for i, l in enumerate(atoms, 1):
    if i in pset:
        zz = float(l[36:44]) + shift
        l = l[:36] + f"{zz:8.3f}" + l[44:]
    newatoms.append(l)

# ---- 물 분자 목록(gro 순서) ----
mols = []          # (첫 원자 index0, 끝 index0+1)
i = 0
while i < n:
    res = newatoms[i][5:10].strip()
    if res in SOLNAMES:
        rn = newatoms[i][:5]; j = i
        while j < n and newatoms[j][:5] == rn and newatoms[j][5:10].strip() == res:
            j += 1
        mols.append((i, j)); i = j
    else:
        i += 1

# ---- 겹치는 물 판정 (PFA 원자와 0.28 nm 이내) ----
pxyz = [(float(newatoms[i-1][20:28]), float(newatoms[i-1][28:36]), float(newatoms[i-1][36:44])) for i in pfa]
def near(x, y, zz):
    for px, py, pz in pxyz:
        dx = x - px; dx -= bx * round(dx / bx)
        dy = y - py; dy -= by * round(dy / by)
        dz = zz - pz; dz -= bz * round(dz / bz)
        if dx * dx + dy * dy + dz * dz < 0.28 ** 2:
            return True
    return False

drop = [False] * len(mols)
for k, (a, b) in enumerate(mols):
    for t in range(a, b):
        l = newatoms[t]
        if near(float(l[20:28]), float(l[28:36]), float(l[36:44])):
            drop[k] = True; break
removed = sum(drop)

keep = [True] * n
first_removed = n + 1
for k, (a, b) in enumerate(mols):
    if drop[k]:
        first_removed = min(first_removed, a + 1)
        for t in range(a, b):
            keep[t] = False

# POL/PFA 색인을 그대로 쓰려면 제거된 원자가 모두 그 뒤에 있어야 한다
guard = max(max(pol), max(pfa))
if removed and first_removed <= guard:
    sys.exit(f"!! 제거된 물({first_removed}번)이 POL/PFA 색인({guard}) 앞에 있습니다 — 색인을 다시 만들어야 합니다")

out = [newatoms[k] for k in range(n) if keep[k]]
fixed = [l[:15] + f"{k % 100000:5d}" + l[20:] for k, l in enumerate(out, 1)]
open(ogro, "w").write(L[0] + "\n" + f"{len(fixed):5d}\n" + "\n".join(fixed) + "\n" + box + "\n")

# ---- topology: SOL 블록마다 따로 갱신 ----
T = open(top).read().splitlines()
newT, inmol, idx, blocks = [], False, 0, []
for l in T:
    st = l.split(';')[0].strip()
    if st.startswith('['):
        inmol = st.strip('[] ').lower() == 'molecules'
        newT.append(l); continue
    p = st.split()
    if inmol and len(p) == 2 and p[0] in SOLNAMES and p[1].lstrip('-').isdigit():
        c = int(p[1])
        if idx + c > len(mols):
            sys.exit(f"!! [molecules] 의 물 개수({idx+c})가 gro 의 물 분자 수({len(mols)})보다 많습니다")
        rm = sum(1 for k in range(idx, idx + c) if drop[k])
        idx += c
        blocks.append((p[0], c, rm))
        newT.append(f"{p[0]}{' ' * max(1, 17 - len(p[0]))}{c - rm}")
    else:
        newT.append(l)
if idx != len(mols):
    sys.exit(f"!! [molecules] 물 합계 {idx} != gro 물 분자 수 {len(mols)}")
open(otop, "w").write("\n".join(newT) + "\n")

# ---- index: System 갱신, POL/PFA 그대로 ----
def fmt(ids):
    return "".join(f"{i}\n" if (k + 1) % 15 == 0 else f"{i} " for k, i in enumerate(ids)) + "\n"
with open(ondx, "w") as f:
    f.write("[ System ]\n" + fmt(range(1, len(fixed) + 1)) + "[ POL ]\n" + fmt(pol) + "[ PFA ]\n" + fmt(pfa))

# ---- 자체 검증 ----
print(f"    물 {removed} 분자 제거, 원자 {n} -> {len(fixed)}")
for name, c, rm in blocks:
    print(f"      [molecules] {name} {c} -> {c - rm}   (제거 {rm})")
assert len(fixed) == n - 3 * removed, f"원자 수 불일치: {len(fixed)} != {n - 3*removed}"
