#!/usr/bin/env bash
# finish_pmf.sh (2026-10-19) — v2.1 결과(10-16)의 '보완 필요' 를 마무리한다. rerun_pmf.sh 가 끝난 뒤에 실행.
#  (1) 4-8 ns 보고서에 남은 '추가 지점' 마다 창을 추가한다. SMD 프레임을 쓰지 않고, 4-8 ns 평균 거리가 가장 가까운
#      기존 창의 마지막 구조에서 출발해 기준 위치를 목표값에 고정한다(pull-coord1-start=no, init=목표; 안쪽 창과 같은 방법).
#      -> v2.1 의 보간이 같은 지점을 반복해서 넣은 문제(SMD 프레임이 없는 구간)를 피한다. 이름 j01, j02 ...
#      최대 FILLROUND(2) 회 반복.
#  (2) 블록 일치 규칙(사전 확정): |dG(4-6) - dG(6-8)| > 2*sqrt(s46^2+s68^2) 이면 전 창을 +4 ns 연장(12 ns)하고
#      4-8 / 8-12 / 4-12 ns 로 다시 분석한다(ECTFE 를 13 ns 로 늘린 것과 같은 조치). 최종값은 4-12 ns.
#  (3) 결과 묶음 ~/pmf_final_fin_MMDD.tgz
# 사용:  tmux new -d -s finish "bash ~/rerun_kit/finish_pmf.sh"      진행: tail -f ~/work/pmf/finish.log
# 환경변수: SYSTEMS="PTFE PVDF"  NT=4  FILLROUND=2  EXTRA_PTFE="1.41 1.58"  EXTRA_PVDF="1.70 1.96 2.47" (보고서 지점에 더해 강제 추가)
set -o pipefail
KITDIR="$(cd "$(dirname "$0")" && pwd)"
W=$HOME/mfpvector; P=$W/pmf; K=$HOME/work/pmf; OUT=${OUT:-$W/pmf_final}
LOG=$K/finish.log; SYSTEMS=${SYSTEMS:-"PTFE PVDF"}; NT=${NT:-4}; FILLROUND=${FILLROUND:-2}
GPU=(-nb gpu -pme gpu -bonded gpu -ntmpi 1 -ntomp "$NT" -pin on)
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate md || { echo "conda 환경 md 없음"; exit 1; }
mkdir -p "$OUT"; cp -n "$KITDIR"/pmf_scripts/* "$K"/ 2>/dev/null
log() { echo "$*  [$(date '+%m-%d %H:%M')]" >> "$LOG"; }
die() { log "!! $*"; exit 1; }
lastt() { if [ -f "$1/pullx.xvg" ]; then grep -v '^[#@]' "$1/pullx.xvg" | tail -1 | awk '{printf "%d",$1}'; else echo 0; fi; }
meanx() { grep -v '^[#@]' "$1/pullx.xvg" | awk -v a=$2 -v b=$3 '$1>=a && $1<=b {s+=$2;n++} END{if(n) printf "%.3f", s/n; else print "nan"}'; }
LASTNS=8     # 현재 창 길이 (ns); 연장하면 12
blocks() {   # $1 $2 $3 = 블록 정의 "a b" x3 (ps)
  cd "$P/$S"; mkdir -p wham2; : > wham2/tpr-all.dat; : > wham2/pullf-all.dat
  for d in win2/n?? win2/i?? win2/j??; do [ -f "$d/umb.gro" ] || continue
    echo "$d/umb.tpr" >> wham2/tpr-all.dat; echo "$d/pullf.xvg" >> wham2/pullf-all.dat; done
  local WT=$(( $(nproc) / 3 )); [ $WT -lt 1 ] && WT=1
  log "$S WHAM 3 구간 시작 ($1 / $2 / $3 ps, $(grep -c . wham2/tpr-all.dat) 창)"
  for B in "$1" "$2" "$3"; do set -- $B
    OMP_NUM_THREADS=$WT gmx wham -it wham2/tpr-all.dat -if wham2/pullf-all.dat -b $1 -e $2 \
      -o wham2/prof_$1_$2.xvg -hist wham2/hist_$1_$2.xvg -unit kJ \
      -bsres wham2/bs_$1_$2.xvg -nBootstrap 200 -bs-method b-hist > wham2/wham_$1_$2.log 2>&1 < /dev/null &
  done
  wait
  : > "$K/report_${T}_v4_blocks.txt"
  for B in "$BL1" "$BL2" "$BL3"; do set -- $B
    [ -s wham2/bs_$1_$2.xvg ] || die "$S WHAM $1-$2 결과 없음 (wham2/wham_$1_$2.log)"
    { echo "== $1-$2 ps"; python3 $K/report_pmf.py ${T}_$1_$2 wham2/prof_$1_$2.xvg wham2/bs_$1_$2.xvg wham2/hist_$1_$2.xvg; } >> "$K/report_${T}_v4_blocks.txt" 2>&1
  done
  log "$S WHAM 완료"
  set -- $BL3; for k in prof bs hist; do cp wham2/${k}_$1_$2.xvg "$OUT/${T}_$k.xvg"; done
  cp "$K/report_${T}_v4_blocks.txt" "$OUT/"
}
secF() { awk -v h="== ${BL3// /-} ps" '$0==h{f=1;next} /^== /{f=0} f' "$K/report_${T}_v4_blocks.txt"; }
dgOf() { awk -v h="== ${1// /-} ps" '$0==h{f=1;next} /^== /{f=0} f && /dG_ads/ {print; exit}' "$K/report_${T}_v4_blocks.txt" | sed -E 's/.*=[ ]*(-?[0-9.]+)[ ]*\+\/-[ ]*([0-9.]+).*/\1 \2/'; }
addwin() {   # $1 = 목표 거리 -> 가장 가까운 기존 창의 마지막 구조에서 init 고정 창 생성·8 ns (또는 LASTNS)
  cd "$P/$S"; local tg=$1 best="" bd=9
  for d in win2/n?? win2/i?? win2/j??; do [ -f "$d/umb.gro" ] || continue
    m=$(meanx $d 4000 $((LASTNS*1000))); [ "$m" = nan ] && continue
    awk -v a=$m -v b=$tg -v c=$bd 'BEGIN{exit !((a-b)^2 < c^2)}' && { bd=$(awk -v a=$m -v b=$tg 'BEGIN{d=a-b; if(d<0)d=-d; print d}'); best=$d; }
  done
  [ -n "$best" ] || die "$S 출발 창을 찾지 못함"
  local j=1 d; while [ -d "win2/$(printf j%02d $j)" ] && [ "$(lastt win2/$(printf j%02d $j))" -ge $((LASTNS*1000-10)) ]; do j=$((j+1)); done
  d=win2/$(printf j%02d $j); mkdir -p $d
  if [ ! -f $d/umb.tpr ]; then
    sed -e "s/^nsteps .*/nsteps                   = $((LASTNS*500000))/" \
        -e "s/^pull-coord1-start .*/pull-coord1-start        = no\npull-coord1-init         = $tg/" win2/umb.mdp > $d/umb.mdp
    gmx grompp -f $d/umb.mdp -c $best/umb.gro -p system_rev2.top -n index_rev2.ndx -o $d/umb.tpr -maxwarn 1 > $d/grompp.log 2>&1 < /dev/null \
      || die "$S ${d#win2/} grompp 실패 ($P/$S/$d/grompp.log)"
    echo "$tg $best" > $d/target.txt
  fi
  CPI=(); [ -f $d/umb.cpt ] && CPI=(-cpi $d/umb.cpt)
  log "$S ${d#win2/} 목표 $tg nm (출발 ${best#win2/}, 평균 $bd 차이) ${LASTNS} ns 시작"
  gmx mdrun -deffnm $d/umb "${CPI[@]}" -px $d/pullx.xvg -pf $d/pullf.xvg "${GPU[@]}" > $d/mdrun.log 2>&1 < /dev/null \
    || die "$S ${d#win2/} mdrun 실패 ($P/$S/$d/mdrun.log)"
  log "   ${d#win2/} 4-${LASTNS} ns 평균 거리 $(meanx $d 4000 $((LASTNS*1000))) nm"
}
log "== finish 시작 SYSTEMS=$SYSTEMS"
for T in $SYSTEMS; do
  S=${T}_PFOA; cd "$P/$S" || die "$P/$S 없음"
  [ -f win2/umb.mdp ] || die "$S win2/umb.mdp 없음 (rerun_pmf.sh 가 끝난 계만 가능)"
  # 현재 창 길이 확인
  LASTNS=8; [ "$(lastt win2/n00)" -ge 11990 ] && LASTNS=12
  BL1="4000 6000"; BL2="6000 8000"; BL3="4000 8000"
  [ $LASTNS -eq 12 ] && { BL1="4000 8000"; BL2="8000 12000"; BL3="4000 12000"; }
  [ -f "$K/report_${T}_v4_blocks.txt" ] || { cp "$K/report_${T}_v3_blocks.txt" "$K/report_${T}_v4_blocks.txt" 2>/dev/null || blocks "$BL1" "$BL2" "$BL3"; }
  # (1) 남은 중첩 부족 + 강제 추가 지점
  EXTRA_VAR="EXTRA_$T"; EXTRA="${!EXTRA_VAR}"
  for r in $(seq 1 $FILLROUND); do
    TG="$(secF | grep -a "추가 지점" | sed -E 's/.*추가 지점 ([0-9.]+).*/\1/' | tr '\n' ' ')"
    [ $r -eq 1 ] && TG="$TG $EXTRA"
    TG=$(echo $TG | tr ' ' '\n' | awk 'NF{printf "%.2f\n",$1}' | sort -u | tr '\n' ' ')
    [ -z "$(echo $TG)" ] && { log "$S 라운드 $r: 중첩 부족 없음"; break; }
    [ -f win2/.fin_r$r ] && { log "$S 라운드 $r 이미 수행 ($(cat win2/.fin_r$r))"; blocks "$BL1" "$BL2" "$BL3"; continue; }
    log "$S 라운드 $r 추가 창: $TG"
    for tg in $TG; do addwin $tg; done
    echo "$TG" > win2/.fin_r$r
    blocks "$BL1" "$BL2" "$BL3"
  done
  # (2) 블록 일치 규칙 -> 필요 시 전 창 +4 ns
  if [ $LASTNS -eq 8 ]; then
    read d46 e46 <<< "$(dgOf "$BL1")"; read d68 e68 <<< "$(dgOf "$BL2")"
    if awk -v a=$d46 -v b=$d68 -v x=$e46 -v y=$e68 'BEGIN{d=a-b; if(d<0)d=-d; exit !(d > 2*sqrt(x*x+y*y))}'; then
      log "$S 블록 차이 |$d46 - $d68| > 2σ -> 전 창 +4 ns 연장 (12 ns), 4-12 ns 로 재분석"
      for d in win2/n?? win2/i?? win2/j??; do [ -f "$d/umb.cpt" ] || continue
        [ "$(lastt $d)" -ge 11990 ] && continue
        gmx convert-tpr -s $d/umb.tpr -until 12000 -o $d/umb.tpr > /dev/null 2>&1 < /dev/null
        log "$S ${d#win2/} -> 12 ns"
        gmx mdrun -deffnm $d/umb -cpi $d/umb.cpt -px $d/pullx.xvg -pf $d/pullf.xvg "${GPU[@]}" > $d/ext12.log 2>&1 < /dev/null \
          || die "$S ${d#win2/} 연장 실패"
      done
      LASTNS=12; BL1="4000 8000"; BL2="8000 12000"; BL3="4000 12000"
      blocks "$BL1" "$BL2" "$BL3"
      # 연장 후 남은 중첩 부족이 있으면 한 번 더 (12 ns 창으로)
      TG="$(secF | grep -a "추가 지점" | sed -E 's/.*추가 지점 ([0-9.]+).*/\1/' | tr '\n' ' ')"
      if [ -n "$(echo $TG)" ] && [ ! -f win2/.fin_r12 ]; then
        log "$S 12 ns 추가 창: $TG"; for tg in $TG; do addwin $tg; done; echo "$TG" > win2/.fin_r12; blocks "$BL1" "$BL2" "$BL3"
      fi
    else log "$S 블록 차이 |$d46 - $d68| ≤ 2σ — 연장 불필요"; fi
  fi
  { echo "$S 최종 (${BL3// /-} ps):"; secF | grep -a "dG_ads\|중첩\|가장자리\|판정"; } >> "$LOG"
  { echo "$S 블록:"; for B in "$BL1" "$BL2" "$BL3"; do echo "  ${B// /-} ps: $(dgOf "$B")"; done; } >> "$LOG"
  log "$S finish 완료 (창 길이 $LASTNS ns)"
done
cp "$LOG" "$OUT/finish.log"; date '+%F %H:%M:%S' > "$OUT/DONE_finish.txt"
STAMP=$(date +%m%d)
tar czf "$HOME/pmf_final_fin_$STAMP.tgz" -C "$W" pmf_final && log "결과 묶음: ~/pmf_final_fin_$STAMP.tgz (남한열에게 전달)"
tar czf "$HOME/pmf_raw_fin_$STAMP.tgz" -C "$W" pmf systems slabs && log "원자료 묶음: ~/pmf_raw_fin_$STAMP.tgz (보관)"
log "== finish 전체 완료"
