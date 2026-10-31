#!/usr/bin/env bash
# ===========================================================================
#  특정 창만 이어서 더 돌린다 (재계산 아님, 뒤에 붙임)
#  사용법:  bash extend_pmf.sh PVDF_PFOA 5 w03 w04 w05
#           bash extend_pmf.sh PVDF_PFOA 5            (모든 창 5 ns 추가)
#  그 뒤 wham 을 다시 돌리려면:  bash run_pmf.sh PVDF_PFOA
# ===========================================================================
set -eo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate md
SYS="${1:?사용법: extend_pmf.sh <계> <추가 ns> [창...]}"
NS="${2:?추가할 ns}"
shift 2
NT="${NT:-8}"
GPU=(-nb gpu -pme gpu -bonded gpu -ntmpi 1 -ntomp "$NT" -pin on)
cd "$HOME/mfpvector/pmf/$SYS"
WINS=("$@"); [ ${#WINS[@]} -eq 0 ] && WINS=($(ls -d ${WINDIR:-windows}/* | xargs -n1 basename))
for w in "${WINS[@]}"; do
  d="${WINDIR:-windows}/$w"
  [ -f "$d/umb.cpt" ] || { echo "$w: 체크포인트 없음 — 건너뜀"; continue; }
  echo "$w  +${NS} ns  $(date '+%H:%M:%S')"
  gmx convert-tpr -s "$d/umb.tpr" -extend "$((NS*1000))" -o "$d/umb.tpr"
  gmx mdrun -deffnm "$d/umb" -cpi "$d/umb.cpt" -px "$d/pullx.xvg" -pf "$d/pullf.xvg" "${GPU[@]}"
done
echo "완료. wham 재실행:  bash run_pmf.sh $SYS"
