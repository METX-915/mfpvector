#!/usr/bin/env python3
"""창 보간: 지정한 거리에 가장 가까운 '아직 안 쓴' SMD 프레임을 골라 window_list 끝에 덧붙인다.
목록 끝에 붙이므로 기존 창(n00..)의 이름과 결과는 그대로 유지된다.
사용: (계 디렉터리에서) infill.py <win디렉터리> <smd디렉터리> <목표1> <목표2> ...
  예: cd ~/mfpvector/pmf/PE_PFOA && python3 infill.py win2 smd2 1.609 1.819 2.053
"""
import sys, os
win, smd = sys.argv[1], sys.argv[2]
targets = [float(a) for a in sys.argv[3:]]
if not targets: sys.exit("!! 목표 거리를 하나 이상 주십시오")
lst = os.path.join(win, "window_list.txt")
lines = [l.split() for l in open(lst) if l.strip()]
rows = [(float(l.split()[0]), float(l.split()[1])) for l in open(os.path.join(smd, "pullx.xvg"))
        if l[0] not in "#@" and len(l.split()) > 1]
times = [r[0] for r in rows]; xs = [r[1] for r in rows]
nfr = len([f for f in os.listdir(smd) if f.startswith("frame") and f.endswith(".gro")])
if nfr < 10: sys.exit(f"!! SMD 프레임이 {nfr} 개뿐입니다")
dtf = times[-1] / (nfr - 1)
cand, j = [], 1
for k in range(nfr):
    tk = k * dtf
    while j < len(times) - 1 and times[j] < tk: j += 1
    t0, t1, x0, x1 = times[j-1], times[j], xs[j-1], xs[j]
    cand.append((k, x0 if t1 <= t0 else x0 + (x1 - x0) * (tk - t0) / (t1 - t0)))
pre = len(f"{smd}/frame")
used = {int(l[0][pre:-4]) for l in lines}
added = []
with open(lst, "a") as f:
    for tg in targets:
        pool = [c for c in cand if c[0] not in used]
        if not pool: sys.exit("!! 남은 프레임이 없습니다")
        k, xx = min(pool, key=lambda c: abs(c[1] - tg))
        if abs(xx - tg) > 0.05:
            print(f"  경고: 목표 {tg:.3f} 에 가장 가까운 프레임이 {xx:.3f} (차이 {abs(xx-tg):.3f} nm)")
        used.add(k); added.append((k, xx))
        f.write(f"{smd}/frame{k}.gro {xx:.3f}\n")
        print(f"  목표 {tg:.3f} -> frame{k} ({xx:.3f} nm)")
print(f"  {len(added)} 개 추가, 총 {len(lines)+len(added)} 창")
