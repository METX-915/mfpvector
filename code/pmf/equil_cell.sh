#!/usr/bin/env bash
# 셀 평형.  v2 (2026-09-15)
#  - 측면 고정: compressibility = 0 4.5e-5
#    계면이 둘인 슬래브에서 x,y 를 1 bar 로 묶으면 표면장력을 0 으로 강제하는 셈이라
#    상자가 끝없이 측면 수축한다. PE 는 4.07 -> 2.62 nm 까지 갔고 수렴하지 않았다.
#    융체는 전단응력을 못 버티므로 x,y 를 고정하고 z 만 커플링하는 것이 옳다.
#  - BASE 로 출발 구조 지정 가능 (기본: 가장 이완된 정방향 창, 없으면 npt.gro)
#  - OUT 으로 출력 이름 지정 가능 (기본 equil -> equil/ 과 equil.gro)
# 사용: BASE=npt.gro bash equil_cell.sh <SYS> [ns]
set -eo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate md
SYS=$1; NS=${2:-8}; OUT=${OUT:-equil}; NT="${NT:-4}"
KIT="$(cd "$(dirname "$0")" && pwd)"; PMF="$HOME/mfpvector/pmf"
GPU=(-nb gpu -pme gpu -bonded gpu -ntmpi 1 -ntomp "$NT" -pin on)
cd "$PMF/$SYS"; mkdir -p "$OUT"

if [ ! -f "$OUT/start.gro" ]; then
  if [ -n "$BASE" ]; then cp "$BASE" "$OUT/start.gro"; echo "  출발 구조: $BASE"
  else
    L=""; for d in windows/w??; do [ -f "$d/umb.gro" ] && L="$d"; done
    if [ -n "$L" ]; then cp "$L/umb.gro" "$OUT/start.gro"; echo "  출발 구조: $L/umb.gro"
    else cp npt.gro "$OUT/start.gro"; echo "  출발 구조: npt.gro"; fi
  fi
fi
echo "  출발 상자: $(tail -1 "$OUT/start.gro")   (x,y 는 여기서 고정됩니다)"

if [ ! -f "$OUT/eq.gro" ]; then
  sed -e "s/^nsteps .*/nsteps                   = $((NS*500000))/" \
      -e "s/^compressibility .*/compressibility          = 0 4.5e-5/" -e "/^pull/d" pull_umbrella.mdp > "$OUT/eq.mdp"
  gmx grompp -f "$OUT/eq.mdp" -c "$OUT/start.gro" -p system.top -n index.ndx -o "$OUT/eq.tpr" -maxwarn 1
  echo "  $NS ns NPT (측면 고정) 시작  $(date '+%H:%M:%S')"
  gmx mdrun -deffnm "$OUT/eq" "${GPU[@]}"
fi

printf "Box-X\nBox-Y\nBox-Z\nDensity\n\n" | gmx energy -f "$OUT/eq.edr" -o "$OUT/box.xvg" > /dev/null 2>&1
echo; echo "  === 상자 추이 (1 ns 평균) ==="
awk '$1!~/^[#@]/{b=int($1/1000); n[b]++; x[b]+=$2; z[b]+=$4; d[b]+=$5}
     END{for(k=0;k<=b;k++) if(n[k]>50) printf "   %3d-%3d ns   x %.3f  z %6.3f   계 밀도 %6.1f\n",k,k+1,x[k]/n[k],z[k]/n[k],d[k]/n[k]}' "$OUT/box.xvg"
cp "$OUT/eq.gro" "$OUT.gro"
echo; echo "  -> $OUT.gro 저장:  $(tail -1 "$OUT.gro")"
python3 "$KIT/surf_z.py" "$OUT.gro" index.ndx 2>/dev/null || true
echo "  z 가 마지막 3 구간에서 0.1 nm 안쪽이면 평형입니다 (x 는 고정이므로 움직이지 않습니다)"
