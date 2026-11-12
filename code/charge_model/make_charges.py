"""Alternative partial-charge sets (Gasteiger-Marsili, MMFF94) for each molecule, in GROMACS itp atom order.
Atom order of <M>.acpype/<M>_GMX.itp follows <M>.acpype/<M>.pdb; bond orders are assigned from the build SMILES."""
import sys, json
from rdkit import Chem
from rdkit.Chem import AllChem, rdPartialCharges
K = sys.argv[1]; out = sys.argv[2]
UNIT = {"PE": ("CC", 15), "PTFE": ("C(F)(F)C(F)(F)", 15), "PVDF": ("CC(F)(F)", 15), "ECTFE": ("CCC(F)(Cl)C(F)(F)", 8)}
PFAS = {"PFOA": "[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
        "PFBS": "[O-]S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
        "GenX": "[O-]C(=O)C(F)(C(F)(F)F)OC(F)(F)C(F)(F)C(F)(F)F",
        "FTOH62": "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"}
MASS = {12: "C", 19: "F", 1: "H", 35: "Cl", 16: "O", 32: "S"}
def itp_atoms(p):
    sec = open(p).read().split("[ atoms ]")[1].split("[")[0]
    rows = [l.split(";")[0].split() for l in sec.splitlines()]
    rows = [r for r in rows if len(r) >= 8]
    return [MASS[round(float(r[7]))] for r in rows], [float(r[6]) for r in rows]
res = {"gasteiger": {}, "mmff94": {}, "am1bcc": {}}
for m in ["PTFE", "PVDF", "ECTFE", "PE", "PFOA", "PFBS", "FTOH62", "GenX"]:
    try:
        el, q0 = itp_atoms(f"{K}/{m}.acpype/{m}_GMX.itp")
    except FileNotFoundError:
        print(m, "no topology, skipped"); continue
    smi = ("[H]" + UNIT[m][0] * UNIT[m][1] + "[H]") if m in UNIT else PFAS[m]
    tmpl = Chem.AddHs(Chem.MolFromSmiles(smi))
    pdb = Chem.MolFromPDBFile(f"{K}/{m}.acpype/{m}.pdb", removeHs=False, sanitize=False, proximityBonding=True)
    assert [a.GetSymbol().capitalize() for a in pdb.GetAtoms()] == el, (m, "element order mismatch")
    mol = AllChem.AssignBondOrdersFromTemplate(Chem.RemoveHs(tmpl, sanitize=True), pdb) if False else None
    # robust route: map template onto pdb graph by substructure match on heavy+H graph
    pdbm = Chem.RWMol(pdb)
    for b in pdbm.GetBonds(): b.SetBondType(Chem.BondType.SINGLE)
    for a in pdbm.GetAtoms(): a.SetNoImplicit(True); a.SetFormalCharge(0)
    q = Chem.MolFromSmarts(Chem.MolToSmarts(tmpl).replace("=", "~").replace("-", "~"))
    qq = Chem.RWMol(tmpl)
    params = Chem.AdjustQueryParameters(); 
    match = pdbm.GetMol().GetSubstructMatch(Chem.MolFromSmarts(Chem.MolToSmarts(Chem.RWMol(tmpl)).replace("=", "~")), useChirality=False) if False else None
    # map by generic graph isomorphism: compare element+connectivity
    from rdkit.Chem import rdFMCS
    gen_t = Chem.RWMol(tmpl)
    for b in gen_t.GetBonds(): b.SetBondType(Chem.BondType.SINGLE)
    for a in gen_t.GetAtoms(): a.SetFormalCharge(0); a.SetNoImplicit(True)
    gp = pdbm.GetMol(); gt = gen_t.GetMol()
    gt.UpdatePropertyCache(False); gp.UpdatePropertyCache(False)
    match = gp.GetSubstructMatch(gt)
    assert len(match) == gt.GetNumAtoms() == gp.GetNumAtoms(), (m, "graph match failed", len(match))
    # transfer bond orders and formal charges from template to pdb-ordered molecule
    newm = Chem.RWMol(gp)
    for b in tmpl.GetBonds():
        i, j = match[b.GetBeginAtomIdx()], match[b.GetEndAtomIdx()]
        newm.GetBondBetweenAtoms(i, j).SetBondType(b.GetBondType())
    for a in tmpl.GetAtoms():
        newm.GetAtomWithIdx(match[a.GetIdx()]).SetFormalCharge(a.GetFormalCharge())
    for a in newm.GetAtoms(): a.SetNoImplicit(True)
    newm = newm.GetMol(); Chem.SanitizeMol(newm)
    rdPartialCharges.ComputeGasteigerCharges(newm)
    g = [float(a.GetProp("_GasteigerCharge")) for a in newm.GetAtoms()]
    props = AllChem.MMFFGetMoleculeProperties(newm)
    mm = [props.GetMMFFPartialCharge(i) for i in range(newm.GetNumAtoms())]
    tot = round(sum(q0))
    # remove rounding residue so each molecule carries exactly its formal charge
    def fix(c): d = (tot - sum(c)) / len(c); return [round(x + d, 6) for x in c]
    res["gasteiger"][m] = fix(g); res["mmff94"][m] = fix(mm); res["am1bcc"][m] = q0
    print(f"{m}: n={len(el)} total {tot}; AM1-BCC F/H means "
          + ", ".join(f"{k} F={sum(c for c,e in zip(res[k][m],el) if e=='F')/max(1,el.count('F')):+.3f} H={sum(c for c,e in zip(res[k][m],el) if e=='H')/max(1,el.count('H')):+.3f}" for k in ("am1bcc","gasteiger","mmff94")))
json.dump(res, open(out, "w"))
