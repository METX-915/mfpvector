#!/usr/bin/env python3
"""
Regenerate the polymer oligomer and PFAS structures that batch_run.py expects
but that were never archived.

Design decision that has to be recorded in the Methods, because the original
chain length is not recoverable from the archive:

    every oligomer is built with 30 backbone carbon atoms.

That keeps the four chains geometrically comparable (same backbone length,
same packmol packing behaviour) while letting the chemistry differ, and it
sits in the middle of the 14-21 monomer range implied by packing 40 chains of
each polymer at its experimental bulk density into the slab region defined in
pack.inp. The resulting slab density is not assumed -- check_density.py
measures it after NPT and compares it with the experimental value, which is
the acceptance test for this rebuild.

Chains are hydrogen/fluorine terminated (no artificial methyl caps beyond the
natural chain ends).

Outputs (matching the paths batch_run.py hardcodes):
    polymers/{PTFE,PVDF,ECTFE,PE}.pdb
    pfas/{PFOA,PFBS,GenX,FTOH62}.pdb

Usage: python3 make_structures.py [--backbone 30] [--outdir .]
"""
import argparse
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

# --------------------------------------------------------------------------
# Polymers. Each entry: (repeat-unit SMILES fragment, backbone C per unit,
# head cap, tail cap).  Chains are written head-to-tail.
# --------------------------------------------------------------------------
POLYMERS = {
    # -(CH2-CH2)-
    "PE":    {"unit": "CC",             "bb": 2, "head": "[H]", "tail": "[H]"},
    # -(CF2-CF2)-
    "PTFE":  {"unit": "C(F)(F)C(F)(F)", "bb": 2, "head": "[H]", "tail": "[H]"},
    # -(CH2-CF2)-  head-to-tail PVDF
    "PVDF":  {"unit": "CC(F)(F)",       "bb": 2, "head": "[H]", "tail": "[H]"},
    # alternating ethylene / chlorotrifluoroethylene: -(CH2-CH2-CFCl-CF2)-
    "ECTFE": {"unit": "CCC(F)(Cl)C(F)(F)", "bb": 4, "head": "[H]", "tail": "[H]"},
}

# --------------------------------------------------------------------------
# PFAS sorbates. Protonation states follow batch_run.py: PFOA/PFBS/GenX are
# modelled as the anions (net charge -1), 6:2 FTOH as the neutral alcohol.
# --------------------------------------------------------------------------
PFAS = {
    # perfluorooctanoate, C7F15COO-
    "PFOA":   {"smiles": "[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
               "charge": -1},
    # perfluorobutanesulfonate, C4F9SO3-
    "PFBS":   {"smiles": "[O-]S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
               "charge": -1},
    # HFPO-DA (GenX) anion, CF3CF2CF2-O-CF(CF3)COO-
    "GenX":   {"smiles": "[O-]C(=O)C(F)(C(F)(F)F)OC(F)(F)C(F)(F)C(F)(F)F",
               "charge": -1},
    # 6:2 fluorotelomer alcohol, C6F13CH2CH2OH
    "FTOH62": {"smiles": "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
               "charge": 0},
}


def build_chain(name, spec, n_backbone):
    """Assemble an oligomer SMILES with the requested backbone carbon count."""
    per_unit = spec["bb"]
    # ECTFE's repeat unit spans 4 backbone carbons, so the target is not always
    # divisible. Round to the nearest achievable chain and report what was built
    # rather than silently using a different length.
    n_units = max(1, round(n_backbone / per_unit))
    smi = spec["head"] + spec["unit"] * n_units + spec["tail"]
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise SystemExit(f"{name}: RDKit could not parse the assembled SMILES:\n  {smi}")
    return mol, n_units, smi


def embed(mol, name, seed=20260710):
    """3D embed + MMFF (fall back to UFF) minimise."""
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = True          # long chains fail without this
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise SystemExit(f"{name}: 3D embedding failed")
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=4000)
        ff = "MMFF94"
    except Exception:
        AllChem.UFFOptimizeMolecule(mol, maxIters=4000)
        ff = "UFF"
    return mol, ff


def formula_and_mass(mol):
    from rdkit.Chem import rdMolDescriptors, Descriptors
    return rdMolDescriptors.CalcMolFormula(mol), Descriptors.MolWt(mol)


def write_pdb(mol, path, resname):
    for atom in mol.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName(resname.ljust(3)[:3])
        info.SetResidueNumber(1)
        info.SetName(f" {atom.GetSymbol():<3}"[:4])
        info.SetIsHeteroAtom(False)
        atom.SetMonomerInfo(info)
    Chem.MolToPDBFile(mol, path, flavor=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", type=int, default=30,
                    help="backbone carbon atoms per polymer chain (default 30)")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    os.makedirs(os.path.join(args.outdir, "polymers"), exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "pfas"), exist_ok=True)

    print(f"target backbone carbons per chain: {args.backbone}")
    print()
    print(f"{'name':8s} {'units':>6s} {'bbC':>5s} {'formula':>22s} {'MW':>9s} {'atoms':>6s} {'x40':>7s}  ff")
    print("-" * 84)

    rows = []
    for name, spec in POLYMERS.items():
        mol, n_units, smi = build_chain(name, spec, args.backbone)
        mol, ff = embed(mol, name)
        formula, mw = formula_and_mass(mol)
        n_at = mol.GetNumAtoms()
        write_pdb(mol, os.path.join(args.outdir, "polymers", f"{name}.pdb"), name)
        bb = n_units * spec["bb"]
        print(f"{name:8s} {n_units:6d} {bb:5d} {formula:>22s} {mw:9.2f} {n_at:6d} {n_at*40:7d}  {ff}")
        rows.append((name, n_units, formula, mw, n_at))

    print()
    print(f"{'name':8s} {'charge':>6s} {'formula':>22s} {'MW':>9s} {'atoms':>6s}  ff")
    print("-" * 68)
    for name, spec in PFAS.items():
        mol = Chem.MolFromSmiles(spec["smiles"])
        if mol is None:
            raise SystemExit(f"{name}: bad SMILES")
        q = Chem.GetFormalCharge(mol)
        if q != spec["charge"]:
            raise SystemExit(f"{name}: formal charge {q} != expected {spec['charge']}")
        mol, ff = embed(mol, name)
        formula, mw = formula_and_mass(mol)
        write_pdb(mol, os.path.join(args.outdir, "pfas", f"{name}.pdb"), name)
        print(f"{name:8s} {q:6d} {formula:>22s} {mw:9.2f} {mol.GetNumAtoms():6d}  {ff}")

    # ---- slab bookkeeping: what packmol will be asked to fit ----
    print()
    print("Slab check (pack.inp region 2-38 A x 2-38 A x 2-28 A = 33.7 nm^3,")
    print("editconf box 4x4x8 nm; NPT then sets the true density):")
    EXP = {"PTFE": 2.20, "PVDF": 1.78, "ECTFE": 1.68, "PE": 0.94}
    NA = 6.02214076e23
    for name, n_units, formula, mw, n_at in rows:
        mass_g = 40 * mw / NA
        for label, vol_nm3 in (("packmol region", 33.7), ("slab in 4x4x2.6 nm", 41.6)):
            rho = mass_g / (vol_nm3 * 1e-21)
            if label == "slab in 4x4x2.6 nm":
                print(f"  {name:6s} 40 chains -> {rho:5.2f} g/cm3 in {label}"
                      f"   (experimental {EXP[name]:.2f})")


if __name__ == "__main__":
    sys.exit(main())
