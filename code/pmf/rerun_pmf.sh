#!/usr/bin/env bash
# rerun_pmf.sh v2.2 (2026-10-19) — v2.1 버그 수정: blocks() 에 mkdir -p wham2, WHAM 결과 없으면 중단 (강구 중령 보고)
# rerun_pmf.sh v2 (2026-10-10) — PTFE·PVDF + PFOA 흡착 PMF 를 계 구축부터 최종 WHAM 까지 한 번에.
# v2: 계 구축은 원고에 쓰인 원본 코드(batch_run.py + patch_build.py, compact.sh 1단계)를 그대로 사용.
#     v1 은 패치 전 구버전 batch_run.build(40 사슬, packmol 직접)를 써서 원고의 계와 달랐음.
# ECTFE·PE 최종값(JHM\final)과 같은 프로토콜:
#   8 ns 셀 평형(측면 고정, z 만 C-rescale) -> 창 범위 z_surf-0.32 ~ min(0.44*box_z-0.04, z_surf+1.85)
#   -> 안쪽 SMD 0.002 nm/ps -> 0.13 nm 간격 창, k=1000 -> 창당 8 ns(3+5) -> 앞 4 ns 폐기, 4-8 ns WHAM(b-hist 200)
#   -> 중첩 <10% 구간 보간 -> 최소점이 안쪽 가장자리면 안쪽 창 2개(-0.12, -0.24 nm)
# 사용:  tmux new -d -s rerun "bash ~/rerun_kit/rerun_pmf.sh"     진행: tail -f ~/work/pmf/rerun.log
#        끊겨도 같은 명령으로 다시 실행하면 끝난 단계는 건너뛰고 이어서 진행한다.
# 환경변수: SYSTEMS="PTFE PVDF" (기본)  NT=4 (OpenMP 스레드)  MAXROUND=3 (보간 반복 한도)  OUT=~/mfpvector/pmf_final
# 끝나면 ~/pmf_final_MMDD.tgz (결과, 전달용) 와 ~/pmf_raw_MMDD.tgz (원자료, 보관용) 를 만든다.
set -o pipefail
KITDIR="$(cd "$(dirname "$0")" && pwd)"
W=$HOME/mfpvector; P=$W/pmf; K=$HOME/work/pmf
OUT=${OUT:-$W/pmf_final}; LOG=$K/rerun.log
SYSTEMS=${SYSTEMS:-"PTFE PVDF"}; NT=${NT:-4}; MAXROUND=${MAXROUND:-3}
GPU=(-nb gpu -pme gpu -bonded gpu -ntmpi 1 -ntomp "$NT" -pin on)
mkdir -p "$W" "$P" "$K" "$OUT"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
log() { echo "$*  [$(date '+%m-%d %H:%M')]" >> "$LOG"; }
die() { log "!! $*"; exit 1; }
cp -n "$KITDIR"/pmf_scripts/* "$K"/ 2>/dev/null

# ---------------------------------------------------------------- 0. 사전 점검 (--check: 화면에 바로 결과 표시 후 종료)
preflight() {
  local bad=0
  ok()  { printf "  [OK]   %s\n" "$*"; }
  ng()  { printf "  [실패] %s\n" "$*"; bad=1; }
  echo "== 사전 점검"
  command -v nvidia-smi >/dev/null && nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | sed 's/^/  GPU: /' || ng "nvidia-smi 없음 (NVIDIA 드라이버 확인)"
  conda activate md 2>/dev/null && ok "conda 환경 md" || ng "conda 환경 md 없음 → bash ~/rerun_kit/1_setup.sh"
  if command -v gmx >/dev/null; then
    gmx --version 2>/dev/null | grep -q "GPU support:.*CUDA" && ok "$(gmx --version 2>/dev/null | grep 'GROMACS version' | xargs) / CUDA" || ng "GROMACS 에 CUDA 지원 없음"
  else ng "gmx 없음"; fi
  command -v packmol >/dev/null && ok "packmol" || ng "packmol 없음"
  conda activate top 2>/dev/null && ok "conda 환경 top" || ng "conda 환경 top 없음"
  command -v acpype >/dev/null && ok "acpype" || ng "acpype 없음"
  python -c "import rdkit" 2>/dev/null && ok "rdkit" || ng "rdkit 없음"
  command -v tmux >/dev/null && ok "tmux" || ng "tmux 없음 → sudo apt install -y tmux"
  local free=$(df -Pk "$HOME" | awk 'NR==2{printf "%d",$4/1048576}')
  [ "$free" -ge 30 ] && ok "여유 디스크 ${free} GB" || ng "여유 디스크 ${free} GB (30 GB 이상 필요)"
  [ $bad -eq 0 ] && echo "== 사전 점검 통과 — 재계산을 시작합니다" || echo "== 사전 점검 실패 — 위 [실패] 항목을 해결한 뒤 다시 실행하십시오"
  return $bad
}
if [ "$1" = "--check" ]; then preflight; exit $?; fi
log "== rerun 시작  SYSTEMS=$SYSTEMS  NT=$NT"
preflight >> "$LOG" 2>&1 || die "사전 점검 실패 (로그의 [실패] 항목 확인)"

# ---------------------------------------------------------------- 1. 토폴로지 (GAFF2 / AM1-BCC, 30 백본 탄소)
# 원고·SI 와 같은 사슬: PTFE (CF2CF2)15 = 92 원자, PVDF (CH2CF2)15 = 92 원자, PFOA = 25 원자
declare -A NAT=([PTFE]=92 [PVDF]=92 [PFOA]=25)
declare -A RHO=([PTFE]=1.883 [PVDF]=1.482)      # SI Table S1 건조 슬래브 밀도 (g/cm3)
declare -A ZS0=([PTFE]=1.717 [PVDF]=1.782)      # 지난 실행 z_surf (참고 비교용)
for n in $SYSTEMS; do [ -n "${RHO[$n]}" ] || die "지원하지 않는 계: $n (PTFE, PVDF 만)"; done
SLAB="$KITDIR/build_kit/slab_tools.py"
REUSE=${REUSE_DIR:-$HOME/mfpvector_c12}
conda activate top || die "conda 환경 top 없음 (1_setup.sh 먼저)"
cd "$W"
for f in make_structures.py compact_conf.py polymers pfas build; do [ -e "$f" ] || cp -r "$KITDIR/build_kit/$f" .; done
[ -f polymers/.made ] || { python make_structures.py --backbone 30 > "$K/make_structures.log" 2>&1 && touch polymers/.made; } \
  || die "make_structures 실패 ($K/make_structures.log)"
for n in $SYSTEMS PFOA; do
  if [ "$(python3 "$SLAB" natoms "$n.acpype/${n}_GMX.itp" 2>/dev/null)" = "${NAT[$n]}" ]; then continue; fi
  # 이전에 만든 30 탄소 토폴로지가 있으면 재사용 (acpype PTFE 는 수 시간 걸림)
  found=""
  if [ -d "$REUSE" ]; then
    for f in $(find "$REUSE" -maxdepth 4 -name "${n}_GMX.itp" 2>/dev/null); do
      [ "$(python3 "$SLAB" natoms "$f")" = "${NAT[$n]}" ] && { found=$(dirname "$f"); break; }
    done
  fi
  if [ -n "$found" ]; then
    rm -rf "$n.acpype"; cp -r "$found" "$n.acpype"; log "토폴로지 재사용 $n <- $found (${NAT[$n]} 원자)"
  else
    [ -f polymers/.compact ] && die "$n 토폴로지가 없는데 사슬 구조가 이미 compact 로 바뀜 — ~/mfpvector/polymers 를 지우고 다시 실행"
    if [ "$n" = PFOA ]; then m=pfas/PFOA.pdb; q=-1; else m=polymers/$n.pdb; q=0; fi
    log "acpype $n 시작 (PTFE 는 수 시간)"
    rm -rf "$n.acpype"
    acpype -i "$m" -b "$n" -n "$q" -a gaff2 -c bcc -s 86400 > "$K/acpype_$n.log" 2>&1 || die "acpype $n 실패 ($K/acpype_$n.log)"
  fi
  [ "$(python3 "$SLAB" natoms "$n.acpype/${n}_GMX.itp")" = "${NAT[$n]}" ] || die "$n 토폴로지 원자 수가 ${NAT[$n]} 이 아님"
done
# 원본 compact.sh 1단계: ETKDG 300 개 중 bounding box 가 가장 작은 구조로 좌표만 교체 (토폴로지는 그대로)
if [ ! -f polymers/.compact ]; then
  python compact_conf.py $SYSTEMS >> "$LOG" 2>&1 && touch polymers/.compact || die "compact 구조 생성 실패"
fi

# ---------------------------------------------------------------- 2. 슬래브 + 스크린 계 (원본 코드 그대로)
# 원 저장소 batch_run.py + patch_build.py(2026-09, 원고 SI 1절의 2단계 구축) 를 그대로 사용한다.
#   packmol(1.5 A, PTFE 4x4x15.2 / PVDF 4x4x14.8 nm) -> 진공 EM
#   -> 1 ns 300 K 200 bar (tau-p 1) -> 2 ns 500 K -> 3 ns 300 K -> 1 ns 300 K (1 bar, tau-p 2, 단계마다 새 속도)
#   -> 바닥 0.2 nm 위로 이동, 4x4x8 nm -> PFOA insert-molecules (표면 +0.7 nm) -> 용매화·중화 -> EM -> NPT 200 ps
conda activate md || die "conda 환경 md 없음"
gmx --version 2>/dev/null | grep -E "GROMACS version|GPU support" >> "$LOG"
mkdir -p "$W/code"
[ -f "$W/code/batch_run.py" ] || cp "$KITDIR/build_kit/batch_run.py" "$W/code/"
cp -n "$KITDIR/build_kit/patch_build.py" "$KITDIR/build_kit/build_one.py" "$W/" 2>/dev/null
grep -q SLAB_BOXZ "$W/code/batch_run.py" || { (cd "$W" && python3 patch_build.py >> "$LOG" 2>&1) && grep -q SLAB_BOXZ "$W/code/batch_run.py"; } \
  || die "batch_run.py 패치 실패"
for poly in $SYSTEMS; do
  if [ ! -f "$W/systems/${poly}_PFOA/npt.gro" ]; then
    log "슬래브·스크린 계 구축 $poly+PFOA (원본 batch_run.build, 슬래브 7 ns 포함)"
    ( cd "$W" && python3 build_one.py "$poly" > "$K/build_$poly.log" 2>&1 ) || die "$poly 구축 실패 ($K/build_$poly.log, $W/slabs/$poly)"
  fi
  # 판정: 건조 슬래브 내부 밀도(밀도 분포 평탄부)를 SI 값과 비교. 원본 density.txt(마지막 프레임 상자 기준)도 기록
  D=$W/slabs/$poly
  if [ ! -f "$D/check.txt" ]; then
    NP=$(sed -n 2p "$D/slab.gro" | tr -d ' ')
    python3 "$SLAB" index "$D/slab.gro" $NP 0 "$D/pol.ndx" > /dev/null
    RI=$(python3 "$K/surf_z.py" "$D/slab.gro" "$D/pol.ndx" | grep -o "내부밀도 *[0-9.]*" | awk '{print $2}')
    DEV=$(awk -v a="$RI" -v b=${RHO[$poly]} 'BEGIN{printf "%+.1f", 100*(a/1000-b)/b}')
    echo "$poly dry slab interior $RI kg/m3 (SI ${RHO[$poly]} g/cm3, $DEV %); 원본 기록(마지막 프레임): $(cat "$D/density.txt" 2>/dev/null)" > "$D/check.txt"
    log "$(cat "$D/check.txt")"
    awk -v x=$DEV 'BEGIN{exit !(x<=5 && x>=-5)}' || { rm -f "$D/check.txt"; die "$poly 건조 슬래브 내부 밀도가 SI 값과 5% 넘게 다름 — 중단 (이 로그를 남한열에게 전달)"; }
  fi
  log "스크린 계 $poly+PFOA: 상자 $(tail -1 "$W/systems/${poly}_PFOA/npt.gro")"
done

# ---------------------------------------------------------------- 3. PMF
lastt() { if [ -f "$1/pullx.xvg" ]; then grep -v '^[#@]' "$1/pullx.xvg" | tail -1 | awk '{printf "%d",$1}'; else echo 0; fi; }
complete8() {
  ( cd "$K" && NO_REPORT=1 NT=$NT bash rebuild_windows.sh "$S" $ZSTOP $ZMAX 0.13 3 >> "rb_$S.log" 2>&1 ) || { log "$S 창 구성/3 ns 실패 ($K/rb_$S.log)"; return 1; }
  cd "$P/$S"; local i=0
  while read -r fr c; do d=$(printf "win2/n%02d" $i); i=$((i+1))
    t=$(lastt $d); [ "$t" -ge 7990 ] && continue
    gmx convert-tpr -s $d/umb.tpr -until 8000 -o $d/umb.tpr > /dev/null 2>&1 < /dev/null
    log "$S ${d#win2/} ($c nm)  $t -> 8000 ps"
    gmx mdrun -deffnm $d/umb -cpi $d/umb.cpt -px $d/pullx.xvg -pf $d/pullf.xvg "${GPU[@]}" > $d/ext.log 2>&1 < /dev/null \
      || { log "$S ${d#win2/} MDRUN 실패"; return 1; }
  done < win2/window_list.txt
}
blocks() {      # 3 구간 WHAM 을 동시에 (CPU 코어를 나눠 씀) — 분석 조건은 그대로
  cd "$P/$S"; mkdir -p wham2; : > wham2/tpr-all.dat; : > wham2/pullf-all.dat
  for d in win2/n?? win2/i??; do [ -f "$d/umb.gro" ] || continue
    echo "$d/umb.tpr" >> wham2/tpr-all.dat; echo "$d/pullf.xvg" >> wham2/pullf-all.dat; done
  local WT=$(( $(nproc) / 3 )); [ $WT -lt 1 ] && WT=1
  log "$S WHAM 3 구간 시작 (구간당 $WT 스레드)"
  for B in "4000 6000" "6000 8000" "4000 8000"; do set -- $B
    OMP_NUM_THREADS=$WT gmx wham -it wham2/tpr-all.dat -if wham2/pullf-all.dat -b $1 -e $2 \
      -o wham2/prof_$1_$2.xvg -hist wham2/hist_$1_$2.xvg -unit kJ \
      -bsres wham2/bs_$1_$2.xvg -nBootstrap 200 -bs-method b-hist > wham2/wham_$1_$2.log 2>&1 < /dev/null &
  done
  wait
  for B in "4000 6000" "6000 8000" "4000 8000"; do set -- $B
    [ -s wham2/bs_$1_$2.xvg ] || die "$S WHAM $1-$2 결과 없음 (wham2/wham_$1_$2.log 확인)"; done
  for B in "4000 6000" "6000 8000" "4000 8000"; do set -- $B
    echo "== $1-$2 ps"; python3 $K/report_pmf.py ${T}_$1_$2 wham2/prof_$1_$2.xvg wham2/bs_$1_$2.xvg wham2/hist_$1_$2.xvg
  done > "$K/report_${T}_v3_blocks.txt" 2>&1
  log "$S WHAM 3 구간 완료"
  for k in prof bs hist; do cp wham2/${k}_4000_8000.xvg "$OUT/${T}_$k.xvg"; done
  cp "$K/report_${T}_v3_blocks.txt" "$OUT/"
}
sec48() { awk '/^== 4000-8000/{f=1} f' "$K/report_${T}_v3_blocks.txt"; }

for T in $SYSTEMS; do
  S=${T}_PFOA; D=$P/$S; mkdir -p "$D"; cd "$D"
  for f in system.top index.ndx npt.gro atomtypes.itp water.itp ions.itp; do cp -n "$W/systems/$S/$f" . 2>/dev/null; done
  cp -n "$W/systems/$S/"*_clean.itp . 2>/dev/null
  for m in pull_umbrella pull_smd; do [ -f $m.mdp ] || cp "$K/$m.mdp" $m.mdp; done
  # 3a. 셀 평형 8 ns (측면 고정)
  if [ ! -f equil.gro ]; then
    log "$S 셀 평형 8 ns 시작  (시작 상자 $(tail -1 npt.gro))"
    ( cd "$K" && OUT=equil BASE=npt.gro NT=$NT bash equil_cell.sh $S 8 > eq_$S.log 2>&1 ) || die "$S 평형 실패 ($K/eq_$S.log)"
  fi
  # 3b. pbcatom (슬래브 공간 중심 원자)
  if ! grep -q '^pull-group1-pbcatom' pull_umbrella.mdp; then
    PB=$(python3 "$K/central_atom.py" equil.gro index.ndx POL 2>/dev/null)
    for m in pull_umbrella pull_smd; do echo "pull-group1-pbcatom      = $PB" >> $m.mdp; done
    log "$S pbcatom = $PB"
  fi
  # 3c. 창 범위
  if [ ! -f range.txt ]; then
    ZS=$(python3 "$K/surf_z.py" equil.gro index.ndx | awk '{for(i=1;i<=NF;i++) if($i=="z_surface") print $(i+2)}')
    BZ=$(tail -1 equil.gro | awk '{print $3}')
    python3 -c "zs=$ZS;bz=$BZ;print(f'{zs-0.32:.2f}', f'{min(0.44*bz-0.04, zs+1.85):.2f}', zs)" > range.txt
    log "$S 평형 상자 $(tail -1 equil.gro)  z_surf=$ZS  창 $(cut -d' ' -f1,2 range.txt)"
  fi
  read ZSTOP ZMAX ZS < range.txt
  # 3d. 창 전부 8 ns
  complete8 || { log "$S 재시도 (60 s 뒤, 체크포인트에서 이어서)"; sleep 60; complete8; } || { log "!! $S 중단 — 이 로그를 남한열에게 전달"; continue; }
  log "$S 전 창 8 ns 완료"
  L=$(ls -d win2/n?? | tail -1)
  { echo "$T z_surf(equil)=$ZS"; echo "last window box: $(tail -1 $L/umb.gro)"; python3 "$K/surf_z.py" $L/umb.gro index_rev2.ndx; } > "$OUT/${T}_cell.txt"
  RHOC=$(grep -o "내부밀도 *[0-9.]*" "$OUT/${T}_cell.txt" | awk '{print $2}')
  RHOD=$(awk -v r=${RHO[$T]} 'BEGIN{print r*1000}')
  awk -v r="$RHOC" -v d="$RHOD" 'BEGIN{exit !(r>0.95*d && r<1.05*d)}' \
    && log "$S 셀 점검 정상: 내부 밀도 $RHOC kg/m3 (SI 건조 슬래브 $RHOD), z_surf $ZS nm (지난 실행 ${ZS0[$T]})" \
    || log "$S 셀 점검 주의: 내부 밀도 $RHOC kg/m3, SI 건조 슬래브 $RHOD 과 5% 이상 차이 — 계산은 계속, 결과 해석 시 확인"
  cp "$W/slabs/$T/check.txt" "$OUT/${T}_dryslab.txt" 2>/dev/null
  blocks
  # 3e. 중첩 부족 보간 — 남은 구간이 없어질 때까지 최대 MAXROUND 회 자동 반복
  for r in $(seq 1 $MAXROUND); do
    cd "$D"
    TG=$(sec48 | grep -a "추가 지점" | sed -E 's/.*추가 지점 ([0-9.]+).*/\1/' | tr '\n' ' ')
    [ -z "$TG" ] && { log "$S 중첩 부족 없음 (라운드 $r 시점)"; break; }
    if [ ! -f win2/.infill_r$r ]; then
      cp win2/window_list.txt win2/window_list.pre_r$r
      python3 "$K/infill.py" win2 smd2 $TG >> "$LOG" 2>&1
      echo "$TG" > win2/.infill_r$r            # 재실행 시 같은 라운드 중복 보간 방지
      log "$S 보간 라운드 $r: $TG"
    fi
    complete8 || { log "$S 재시도 (60 s 뒤)"; sleep 60; complete8; } || { log "!! $S 중단 — 이 로그를 남한열에게 전달"; continue 2; }
    blocks
  done
  # 3f. 안쪽 가장자리 최소점 -> 안쪽 창 2개
  cd "$D"
  if sec48 | grep -aq "안쪽 가장자리"; then
    read ZI WI <<< "$(awk '{print $2, NR-1}' win2/window_list.txt | sort -g | head -1)"
    PREV=$(printf "win2/n%02d" $WI)/umb.gro; j=1
    log "$S 최소점 안쪽 가장자리 -> 가장 안쪽 창 $ZI nm 기준 안쪽 창 추가"
    for dz in 0.12 0.24; do
      TGT=$(awk -v z=$ZI -v d=$dz 'BEGIN{printf "%.3f", z-d}'); d=$(printf "win2/i%02d" $j); j=$((j+1)); mkdir -p $d
      if [ "$(lastt $d)" -lt 7990 ]; then
        if [ ! -f $d/umb.tpr ]; then
          sed -e "s/^nsteps .*/nsteps                   = 4000000/" \
              -e "s/^pull-coord1-start .*/pull-coord1-start        = no\npull-coord1-init         = $TGT/" win2/umb.mdp > $d/umb.mdp
          gmx grompp -f $d/umb.mdp -c $PREV -p system_rev2.top -n index_rev2.ndx -o $d/umb.tpr -maxwarn 1 > $d/grompp.log 2>&1 < /dev/null \
            || { log "$S ${d#win2/} grompp 실패"; break; }
        fi
        CPI=(); [ -f $d/umb.cpt ] && CPI=(-cpi $d/umb.cpt)
        log "$S ${d#win2/} init $TGT nm 8 ns"
        gmx mdrun -deffnm $d/umb "${CPI[@]}" -px $d/pullx.xvg -pf $d/pullf.xvg "${GPU[@]}" > $d/mdrun.log 2>&1 < /dev/null \
          || { log "$S ${d#win2/} mdrun 실패"; break; }
      fi
      log "   ${d#win2/} 4-8 ns 평균 거리 $(grep -v '^[#@]' $d/pullx.xvg | awk '$1>=4000{s+=$2;n++}END{if(n)printf "%.3f",s/n}') nm"
      PREV=$d/umb.gro
    done
    blocks
  else log "$S 최소점 가장자리 아님"; fi
  { echo "$S 최종 (4-8 ns):"; sec48 | grep -a "dG_ads\|중첩\|가장자리\|판정"; } >> "$LOG"
  if sec48 | grep -aq "판정: 사용 가능"; then log "$S 최종 판정: 사용 가능"
  else log "$S 최종 판정: 보완 필요 — 자동 보완 한도 도달. 더 돌리지 말고 이 로그를 남한열에게 전달"; fi
  log "$S 완료"
done
# ---------------------------------------------------------------- 4. 결과 묶기
cp "$LOG" "$OUT/"; date '+%F %H:%M:%S' > "$OUT/DONE.txt"
STAMP=$(date +%m%d)
tar czf "$HOME/pmf_final_$STAMP.tgz" -C "$W" pmf_final && log "결과 묶음: ~/pmf_final_$STAMP.tgz  (남한열에게 전달)"
tar czf "$HOME/pmf_raw_$STAMP.tgz" -C "$W" pmf systems && log "원자료 묶음: ~/pmf_raw_$STAMP.tgz  (외장 저장장치에 보관, 삭제 금지)"
log "== 전체 완료"
