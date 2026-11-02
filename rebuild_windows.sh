#!/usr/bin/env bash
# 평형된 셀(equil.gro) 하나에서 umbrella 창 전체를 새로 만든다.  v2 (2026-09-15)
#  - 창 선택: pullx 를 프레임 시각에 선형보간한다 (v1 은 "시각이 프레임간격의 정수배" 를
#    부동소수 등호로 검사해 후보가 2개만 남았다. 창이 전부 양끝에 몰려 20시간을 버렸다.)
#  - 창 목록을 쓴 직후 퇴화 검사: 서로 다른 좌표가 창 수의 80% 미만이거나 실제 범위가
#    요청 범위의 80% 미만이면 즉시 중단한다.
#  - 측면 고정: compressibility = 0 4.5e-5. 계면이 둘인 슬래브에서 측면 압력을 1 bar 로
#    묶으면 표면장력을 0 으로 강제해 상자가 끝없이 수축한다(PE 는 4.07 -> 2.62 까지 갔다).
#  - 기본 -ntomp 4 (벤치마크: 71.1 vs 63.0 ns/day)
# 사용: bash rebuild_windows.sh <SYS> <ZSTOP> <ZMAX> [STEP] [NS]      K=1000 등 환경변수 가능
set -eo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate md
SYS=$1; ZSTOP=$2; ZMAX=$3; STEP=${4:-0.13}; NS=${5:-3}; K=${K:-1000}; RATE=${RATE:-0.002}
KIT="$(cd "$(dirname "$0")" && pwd)"; PMF="$HOME/mfpvector/pmf"; NT="${NT:-4}"
GPU=(-nb gpu -pme gpu -bonded gpu -ntmpi 1 -ntomp "$NT" -pin on)
cd "$PMF/$SYS"
BASE="${BASE:-equil.gro}"
[ -f "$BASE" ] || { echo "!! $BASE 가 없습니다. 먼저 equil_cell.sh 를 돌리십시오"; exit 1; }
BZ=$(tail -1 "$BASE" | awk '{print $3}')
awk -v z="$ZMAX" -v b="$BZ" 'BEGIN{if(z>0.44*b){printf "!! ZMAX %.2f > 0.44*box_z(%.2f)\n",z,0.44*b; exit 1}}'
echo "== $SYS 창 재구성: $ZMAX -> $ZSTOP nm, 간격 $STEP, 창당 $NS ns, k=$K, ntomp=$NT  (상자 z $BZ)"
mkdir -p rev2 smd2 win2
FIXXY='s/^compressibility .*/compressibility          = 0 4.5e-5/'

# 1) 벌크 출발 구조
if [ ! -f rev2/start.gro ]; then
  python3 "$KIT/shift_pfa.py" "$BASE" index.ndx system.top "$ZMAX" rev2/shifted.gro system_rev2.top index_rev2.ndx
  printf "integrator=steep\nnsteps=2000\nemtol=500\nemstep=0.01\nnstlist=10\ncutoff-scheme=Verlet\ncoulombtype=PME\nrcoulomb=1.0\nrvdw=1.0\npbc=xyz\n" > rev2/em.mdp
  gmx grompp -f rev2/em.mdp -c rev2/shifted.gro -p system_rev2.top -n index_rev2.ndx -o rev2/em.tpr -maxwarn 1
  gmx mdrun -deffnm rev2/em -nb gpu -ntmpi 1 -ntomp "$NT"
  cp rev2/em.gro rev2/start.gro
fi

# 2) 안쪽으로 SMD
X0=$(python3 "$KIT/com_dist.py" rev2/start.gro index_rev2.ndx POL PFA)
NSTEPS=$(awk -v x="$X0" -v s="$ZSTOP" -v r="$RATE" 'BEGIN{printf "%d",(x-s+0.05)/r/0.002}')
echo "    출발 거리 $X0 nm -> 목표 $ZSTOP nm, 속도 $RATE nm/ps ($NSTEPS steps = $(awk -v n=$NSTEPS 'BEGIN{printf "%.0f",n*0.002}') ps)"
if [ ! -f smd2/smd.gro ]; then
  sed -e "s/^nsteps .*/nsteps                   = $NSTEPS/" \
      -e "s/^pull-coord1-rate .*/pull-coord1-rate         = -$RATE/" \
      -e "s/^gen-vel .*/gen-vel                  = yes\ngen-temp                 = 300\ngen-seed                 = 4242/" \
      -e "$FIXXY" pull_smd.mdp > smd2/smd.mdp
  gmx grompp -f smd2/smd.mdp -c rev2/start.gro -p system_rev2.top -n index_rev2.ndx -o smd2/smd.tpr -maxwarn 1
  gmx mdrun -deffnm smd2/smd -px smd2/pullx.xvg -pf smd2/pullf.xvg "${GPU[@]}"
  ( cd smd2 && gmx trjconv -s smd.tpr -f smd.xtc -o frame.gro -sep <<< "0" )
fi

# 3) 창 고르기 (+ 퇴화 검사)
if [ ! -f win2/window_list.txt ]; then
python3 - "$ZSTOP" "$STEP" "$ZMAX" <<'PY'
import sys, os
zstop, step, zmax = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
rows = [(float(l.split()[0]), float(l.split()[1])) for l in open("smd2/pullx.xvg")
        if l[0] not in "#@" and len(l.split()) > 1]
times = [r[0] for r in rows]; xs = [r[1] for r in rows]
nfr = len([f for f in os.listdir("smd2") if f.startswith("frame") and f.endswith(".gro")])
if nfr < 10: sys.exit(f"!! SMD 프레임이 {nfr} 개뿐입니다")
dtf = times[-1] / (nfr - 1)
# 각 프레임 시각에서 pullx 를 선형보간 -> (프레임번호, 거리)
cand, j = [], 1
for k in range(nfr):
    tk = k * dtf
    while j < len(times) - 1 and times[j] < tk: j += 1
    t0, t1, x0, x1 = times[j-1], times[j], xs[j-1], xs[j]
    x = x0 if t1 <= t0 else x0 + (x1 - x0) * (tk - t0) / (t1 - t0)
    cand.append((k, x))
targets, x = [], zstop
while x <= max(c[1] for c in cand) + 1e-6: targets.append(x); x += step
out, used = [], set()
for tg in targets:
    pool = [c for c in cand if c[0] not in used] or cand
    k, xx = min(pool, key=lambda c: abs(c[1] - tg))
    used.add(k); out.append((f"smd2/frame{k}.gro", xx))
out.sort(key=lambda w: w[1])
coords = [w[1] for w in out]
uniq = len({round(c, 3) for c in coords})
reach = min(c[1] for c in cand)
print(f"    창 {len(out)} 개: {min(coords):.2f} -> {max(coords):.2f} nm  "
      f"(서로 다른 좌표 {uniq}, SMD 도달 {reach:.3f} nm / 목표 {zstop:.2f})")
if uniq < 0.8 * len(out):
    sys.exit(f"!! 창이 퇴화했습니다: {len(out)} 개 중 서로 다른 좌표가 {uniq} 개뿐입니다")
if reach > zstop + 0.30:
    sys.exit(f"!! SMD 가 {reach:.3f} nm 에서 멈췄습니다 (목표 {zstop:.2f}). "
             f"당김 속도를 낮추십시오: RATE=0.001 로 다시 실행")
gaps = [coords[i+1] - coords[i] for i in range(len(coords) - 1)]
if max(gaps) > 2.5 * step:
    sys.exit(f"!! 창 간격이 최대 {max(gaps):.3f} nm 로 벌어졌습니다 (요청 {step})")
with open("win2/window_list.txt", "w") as f:
    for g, c in out: f.write(f"{g} {c:.3f}\n")
print(f"    간격 최소 {min(gaps):.3f} / 최대 {max(gaps):.3f} nm — 정상")
PY
fi
NW=$(grep -c . win2/window_list.txt)
echo "    총 $NW 창, 예상 소요 약 $(awk -v n=$NW -v s=$NS 'BEGIN{printf "%.1f", n*(s*20.3+1)/60}') 시간"

# 4) umbrella
sed -e "s/^nsteps .*/nsteps                   = $((NS*500000))/" \
    -e "s/^pull-coord1-k .*/pull-coord1-k            = $K/" \
    -e "$FIXXY" pull_umbrella.mdp > win2/umb.mdp
i=0
while read -r frame coord; do
  w=$(printf "n%02d" "$i"); i=$((i+1)); d="win2/$w"; mkdir -p "$d"
  [ -f "$d/umb.gro" ] && { echo "   $w 완료"; continue; }
  [ -f "$d/start.gro" ] || cp "$frame" "$d/start.gro"
  [ -f "$d/umb.tpr" ] || gmx grompp -f win2/umb.mdp -c "$d/start.gro" -p system_rev2.top -n index_rev2.ndx -o "$d/umb.tpr" -maxwarn 1
  CPI=(); [ -f "$d/umb.cpt" ] && CPI=(-cpi "$d/umb.cpt")
  echo "   $w  ($coord nm)  $(date '+%H:%M:%S')"
  gmx mdrun -deffnm "$d/umb" "${CPI[@]}" -px "$d/pullx.xvg" -pf "$d/pullf.xvg" "${GPU[@]}"
done < win2/window_list.txt

# 5) WHAM — win2 만 (NO_REPORT=1 이면 생략: 중간 점검용일 뿐 최종값에는 쓰이지 않음)
if [ -n "$NO_REPORT" ]; then echo "REBUILD_DONE $SYS (중간 WHAM 생략)"; exit 0; fi
mkdir -p wham2; : > wham2/tpr-files.dat; : > wham2/pullf-files.dat
for d in win2/n??; do [ -f "$d/umb.gro" ] || continue; echo "$d/umb.tpr" >> wham2/tpr-files.dat; echo "$d/pullf.xvg" >> wham2/pullf-files.dat; done
gmx wham -it wham2/tpr-files.dat -if wham2/pullf-files.dat -b 200 -o wham2/profile.xvg -hist wham2/histo.xvg \
         -unit kJ -bsres wham2/bsResult.xvg -nBootstrap 200 -bs-method b-hist > wham2/wham.log 2>&1
python3 "$KIT/report_pmf.py" "${SYS}_NEW" wham2/profile.xvg wham2/bsResult.xvg wham2/histo.xvg | tee "$PMF/${SYS}_new_report.txt"
echo "REBUILD_DONE $SYS"
