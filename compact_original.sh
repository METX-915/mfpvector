#!/usr/bin/env bash
# ===========================================================================
#  (1) 고분자 4종의 "가장 뭉친" 형태를 골라 polymers/*.pdb 를 교체
#      - SMILES / 원자 순서 / 결합이 동일하므로 기존 .acpype 토폴로지 그대로 사용
#      - 좌표만 바뀜
#  (2) PTFE·ECTFE 는 상자 높이를 올려가며 packmol 통과 지점을 탐색
#  실행:  cd ~/work && bash compact.sh
# ===========================================================================
set -eo pipefail
source "$HOME/miniforge3/etc/profile.d/conda.sh"

cd "$HOME/work"
mkdir -p backup_conf
cp polymers/*.pdb backup_conf/ 2>/dev/null || true

# ---------------------------------------------------------------- 1. 구조
conda activate top
python - <<'PY'
import os, sys, itertools, math
sys.path.insert(0, os.path.expanduser("~/work"))
import make_structures as MS
from rdkit import Chem
from rdkit.Chem import AllChem

def bbox_vol(mol, cid):
    c = mol.GetConformer(cid)
    p = [tuple(c.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    e = [max(x[k] for x in p) - min(x[k] for x in p) for k in range(3)]
    return e[0]*e[1]*e[2]/1000.0, e          # nm^3, extents in A

print(f"{'poly':7s} {'원자':>4s} {'기존 bbox':>10s} {'새 bbox':>9s} {'감소':>7s}")
print("-"*46)
for name, spec in MS.POLYMERS.items():
    old_n = sum(1 for l in open(f"polymers/{name}.pdb")
                if l.startswith(("ATOM","HETATM")))
    ox = [(float(l[30:38]),float(l[38:46]),float(l[46:54]))
          for l in open(f"polymers/{name}.pdb") if l.startswith(("ATOM","HETATM"))]
    oe = [max(q[k] for q in ox)-min(q[k] for q in ox) for k in range(3)]
    ov = oe[0]*oe[1]*oe[2]/1000.0

    mol, nu, smi = MS.build_chain(name, spec, 30)
    mol = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3(); ps.randomSeed = 20260710; ps.useRandomCoords = True
    ps.pruneRmsThresh = -1.0
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=300, params=ps)
    if not cids:
        sys.exit(f"{name}: 컨포머 생성 실패")
    AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000, numThreads=0)
    best = min(cids, key=lambda c: bbox_vol(mol, c)[0])
    nv, ne = bbox_vol(mol, best)

    if mol.GetNumAtoms() != old_n:
        sys.exit(f"{name}: 원자수 불일치 {mol.GetNumAtoms()} != {old_n} — 중단")

    keep = Chem.Mol(mol); keep.RemoveAllConformers()
    keep.AddConformer(mol.GetConformer(best), assignId=True)
    MS.write_pdb(keep, f"polymers/{name}.pdb", name)
    print(f"{name:7s} {old_n:4d} {ov:9.2f}  {nv:8.2f}  {(1-nv/ov)*100:6.0f}%")
PY

# ---------------------------------------------------------------- 2. packmol
conda activate md
: > compact.out
scan () {                       # $1 이름  $2 사슬수  $3.. 시도할 zt 목록
  local p=$1 n=$2; shift 2
  for z in "$@"; do
    printf 'tolerance 1.5\nfiletype pdb\noutput /tmp/cp.pdb\nstructure polymers/%s.pdb\n  number %s\n  inside box 2. 2. 2. 38. 38. %s\nend structure\n' "$p" "$n" "$z" > /tmp/cp.inp
    local s=$SECONDS
    if packmol < /tmp/cp.inp > /tmp/cp.log 2>&1; then
      echo "$p N=$n zt=$z OK $((SECONDS-s))s h0=$(python3 -c "print(($z+2)/10)")" | tee -a compact.out
      return
    fi
    echo "$p N=$n zt=$z FAIL $((SECONDS-s))s" | tee -a compact.out
  done
  echo "$p 모든 높이 실패" | tee -a compact.out
}
scan PE    65 175.0
scan PVDF  54 146.0
scan ECTFE 42 122.0 160.0 200.0 250.0
scan PTFE  42 115.0 150.0 190.0 240.0
echo DONE | tee -a compact.out
