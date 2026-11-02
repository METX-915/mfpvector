#!/usr/bin/env python3
"""compact.sh 1단계(원본, 2026-09)와 같은 코드: 고분자별 ETKDG 300 개(seed 20260710) -> MMFF94 -> bounding-box 부피가
가장 작은 구조로 polymers/<name>.pdb 좌표만 교체. 원자 순서·결합 동일 -> 기존 acpype 토폴로지 그대로 사용.
사용: (~/mfpvector 에서) python compact_conf.py PTFE PVDF"""
import os, sys
sys.path.insert(0, os.getcwd())
import make_structures as MS
from rdkit import Chem
from rdkit.Chem import AllChem

def bbox_vol(mol, cid):
    c = mol.GetConformer(cid)
    p = [tuple(c.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    e = [max(x[k] for x in p) - min(x[k] for x in p) for k in range(3)]
    return e[0]*e[1]*e[2]/1000.0, e

for name in sys.argv[1:]:
    spec = MS.POLYMERS[name]
    old_n = sum(1 for l in open(f"polymers/{name}.pdb") if l.startswith(("ATOM", "HETATM")))
    ox = [(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in open(f"polymers/{name}.pdb") if l.startswith(("ATOM", "HETATM"))]
    oe = [max(q[k] for q in ox) - min(q[k] for q in ox) for k in range(3)]
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
    print(f"{name}: compact 구조 (bbox {ov:.2f} -> {nv:.2f} nm3, {old_n} 원자)")
