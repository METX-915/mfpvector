#!/usr/bin/env python3
"""
mfpvector batch MD driver — RQ2 relative PFAS adsorption on (fluoro)polymer slabs.
20 systems (5 polymers x 4 PFAS, pristine) x 5 seeds x 6 ns GPU production.
Metric: polymer-PFAS short-range interaction energy (LJ+Coul, via CPU -rerun),
        + min-distance/contact fraction. Absolute values NOT claimed -> relative ranking.
Run in background:  nohup python scripts/batch_run.py > results/batch.out 2>&1 &
"""
import os, subprocess, sys, time, csv, glob

os.chdir(os.environ.get("MFP_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
R="results"; os.makedirs(R, exist_ok=True)
POLYS=os.environ.get("MFP_POLYS","PTFE,PVDF").split(",")
_ALLQ={"PFOA":-1,"PFBS":-1,"GenX":-1,"FTOH62":0}
PFASS={k:_ALLQ[k] for k in os.environ.get("MFP_PFAS","PFBS,FTOH62").split(",")}
SEEDS=[1,2,3,4,5]
NPOLY=40; PROD_NS=6; BXY=40.0
CSV=os.path.join(R,os.environ.get("MFP_CSV","raw_results_realization2.csv"))
LOG=open(os.path.join(R,"progress.log"),"a")

def log(m):
    s=f"[{time.strftime('%m-%d %H:%M:%S')}] {m}"
    print(s,flush=True); LOG.write(s+"\n"); LOG.flush()

def run(cmd, inp=None, shell=False, quiet=True):
    r=subprocess.run(cmd, input=inp, capture_output=True, text=True, shell=shell)
    if r.returncode!=0 and shell and isinstance(cmd,str) and cmd.startswith("packmol"):
        # packmol may exit non-zero ("not all molecules placed") while still writing a usable output;
        # the vacuum EM that follows resolves residual overlaps (same tolerance as the October rebuild)
        inpf=cmd.split("<")[1].strip(); out=[l.split()[1] for l in open(inpf) if l.startswith("output")][0]
        if os.path.exists(out) and os.path.getsize(out)>0:
            log(f"packmol exit {r.returncode}; output present, continuing"); return r
    if r.returncode!=0:
        raise RuntimeError(f"CMD FAIL: {cmd}\nSTDOUT {r.stdout[-1200:]}\nSTDERR {r.stderr[-1200:]}")
    return r

def natoms_pdb(p):
    return sum(1 for l in open(p) if l.startswith(("ATOM","HETATM")))

# ---- topology helpers (proven in build_system.py) ----
def extract_at(t):
    out=[]; inb=False
    for ln in t.splitlines():
        if ln.strip().startswith('['):
            inb=ln.strip().lower().startswith('[ atomtypes'); continue
        if inb and ln.strip() and not ln.strip().startswith(';'): out.append(ln.rstrip())
    return out
def strip_at(t): return t[t.lower().find('[ moleculetype ]'):]
def rename_mol(t,newname):
    lines=t.splitlines(); out=[]; done=False; hdr=False
    for ln in lines:
        if ln.strip().lower().startswith('[ moleculetype'): hdr=True; out.append(ln); continue
        if hdr and not done and ln.strip() and not ln.strip().startswith(';'):
            out.append(f"{newname}    3"); done=True; continue
        out.append(ln)
    return "\n".join(out)+"\n"

WATER=open("build/water.itp").read(); IONS=open("build/ions.itp").read()
EXTRA=[" OW OW 0.00000 0.00000 A 3.15061e-01 6.36386e-01 ; TIP3P O",
       " HW HW 0.00000 0.00000 A 0.00000e+00 0.00000e+00 ; TIP3P H",
       " Na Na 0.00000 0.00000 A 2.43928e-01 3.65846e-01 ; JC08 Na",
       " Cl Cl 0.00000 0.00000 A 4.47766e-01 1.48913e-01 ; JC08 Cl"]

COMMON=("cutoff-scheme=Verlet\ncoulombtype=PME\nrcoulomb=1.0\nrvdw=1.0\npbc=xyz\n"
        "constraints=h-bonds\nconstraint-algorithm=lincs\ntcoupl=V-rescale\n"
        "tc-grps=System\ntau-t=0.1\nref-t=300\n")
NT=str(min(8,os.cpu_count() or 8))
GPU=["-nb","gpu","-pme","gpu","-bonded","gpu","-update","gpu","-ntmpi","1","-ntomp",NT,"-pin","on"]

# ============================================================================
#  2-stage slab build  (replaces the original single-shot packmol build)
#
#  stage 1  packmol at low density -> vacuum EM -> semiisotropic NPT
#           the melt collapses to the force field's own equilibrium density.
#           built ONCE per polymer and cached in slabs/<poly>/, because the
#           dry slab does not depend on which PFAS is adsorbed.
#           the resulting z and density are logged -- this is the physical
#           acceptance test the original pipeline never performed.
#  stage 2  slab -> box 4x4x8 -> insert PFAS above the surface -> solvate
#           -> genion -> EM -> NPT.  identical to the original from solvate on.
# ============================================================================

NPOLY_BY  = {"PE":65, "PVDF":54, "ECTFE":42, "PTFE":42, "PP":40}
HZ0_BY    = {"PE":17.7,"PVDF":14.8,"ECTFE":16.2,"PTFE":15.2,"PP":2.60}
# reference is the AMORPHOUS density, not the semicrystalline bulk value.
# an MD melt has no crystallites, so bulk PTFE 2.20 / PVDF 1.78 are the wrong
# comparison.  Callister: PTFE rho_c 2.28, rho_a 1.98.  PVDF rho_a 1.68.
# ECTFE amorphous is not well documented -- estimated from bulk 1.68 at ~55%
# crystallinity, so treat it as indicative only.
EXP_RHO   = {"PE":0.855,"PVDF":1.68,"ECTFE":1.55,"PTFE":1.98,"PP":0.85}
PACK_TOL  = "1.5"
# compression schedule: (tag, ns, T/K, ref-p, tau-p)
#   cmp1  collapse the void fast under overpressure
#   cmp2  anneal at 500 K -- chains are immobile at 300 K, so densification
#         stalls; this is what makes the melt actually equilibrate
#   cmp3  cool to 300 K and settle
#   cmp4  1 ns probe: if z still drifts, the slab is not converged
CMP_STAGES = (("cmp1", 1.0, 300, "1.0 200.0", "1.0"),
              ("cmp2", 2.0, 500, "1.0 1.0",   "2.0"),
              ("cmp3", 3.0, 300, "1.0 1.0",   "2.0"),
              ("cmp4", 1.0, 300, "1.0 1.0",   "2.0"))
ZMIN_OK   = 2.2    # 압축 후 상자 z 하한 (PBC 최소 = 2 x rvdw = 2.0 nm)
SLAB_BOXZ = 8.0


def itp_mass(path):
    """Sum the mass column of [ atoms ] in a GROMACS .itp (g/mol)."""
    m = 0.0; inat = False
    for ln in open(path):
        s = ln.strip()
        if s.startswith("["):
            inat = s.replace(" ", "").startswith("[atoms]")
            continue
        if not inat or not s or s.startswith(";"):
            continue
        p = s.split()
        if len(p) > 7:
            try: m += float(p[7])
            except ValueError: pass
    return m


def gro_read(path):
    L = open(path).read().splitlines()
    n = int(L[1].strip())
    return L[0], n, L[2:2+n], L[2+n]


def gro_box(path):
    return [float(x) for x in gro_read(path)[3].split()[:3]]


def gro_maxz(path):
    return max(float(l[36:44]) for l in gro_read(path)[2])


def gro_shift(src, dst, bx, by, bz, zmin=0.2):
    """Translate in z so the lowest atom sits at zmin, and set a new box."""
    title, n, at, _ = gro_read(src)
    dz = zmin - min(float(l[36:44]) for l in at)
    out = [l[:36] + "%8.3f" % (float(l[36:44]) + dz) + l[44:] for l in at]
    with open(dst, "w") as f:
        f.write(title + "\n%5d\n" % n)
        f.write("\n".join(out) + "\n")
        f.write("%10.5f%10.5f%10.5f\n" % (bx, by, bz))


def build(poly, pfas, wd):
    os.makedirs(wd, exist_ok=True)
    NPOLY = NPOLY_BY[poly]
    nc = PFASS[pfas]
    poly_itp = open(f"{poly}.acpype/{poly}_GMX.itp").read()
    pfas_itp = open(f"{pfas}.acpype/{pfas}_GMX.itp").read()
    at = {}
    for t in (poly_itp, pfas_itp):
        for ln in extract_at(t): at[ln.split()[0]] = ln
    with open(f"{wd}/atomtypes.itp", "w") as f:
        f.write("[ atomtypes ]\n;name bond_type mass charge ptype sigma epsilon\n")
        for v in at.values(): f.write(v + "\n")
        for e in EXTRA: f.write(e + "\n")
    open(f"{wd}/{poly}_clean.itp", "w").write(rename_mol(strip_at(poly_itp), poly))
    open(f"{wd}/{pfas}_clean.itp", "w").write(rename_mol(strip_at(pfas_itp), pfas))
    open(f"{wd}/water.itp", "w").write(WATER)
    open(f"{wd}/ions.itp", "w").write(IONS)

    INC = ('[ defaults ]\n1 2 yes 0.5 0.8333333333\n\n'
           f'#include "atomtypes.itp"\n#include "{poly}_clean.itp"\n'
           f'#include "{pfas}_clean.itp"\n'
           f'#include "water.itp"\n#include "ions.itp"\n\n')
    EM = ("integrator=steep\nnsteps=5000\nemtol=1000\nemstep=0.01\nnstlist=10\n"
          + COMMON.split('tcoupl')[0])

    # ------------------------------------------------ stage 1: dry slab (cached)
    sd = f"slabs/{poly}"
    if not os.path.exists(f"{sd}/slab.gro"):
        os.makedirs(sd, exist_ok=True)
        for fn in ("atomtypes.itp", f"{poly}_clean.itp", f"{pfas}_clean.itp",
                   "water.itp", "ions.itp"):
            open(f"{sd}/{fn}", "w").write(open(f"{wd}/{fn}").read())
        hz0 = HZ0_BY[poly]
        zt  = hz0 * 10.0 - 2.0
        open(f"{sd}/pack.inp", "w").write(
            f"tolerance {PACK_TOL}\nfiletype pdb\noutput {sd}/dry.pdb\n"
            f"structure polymers/{poly}.pdb\n  number {NPOLY}\n"
            f"  inside box 2. 2. 2. {BXY-2} {BXY-2} {zt:.1f}\nend structure\n")
        log(f"SLAB {poly}: packing {NPOLY} chains into 3.6 x 3.6 x {hz0-0.4:.2f} nm")
        run(f"packmol < {sd}/pack.inp", shell=True)
        run(["gmx", "editconf", "-f", f"{sd}/dry.pdb", "-o", f"{sd}/dry.gro",
             "-box", "4.0", "4.0", f"{hz0:.3f}"])
        open(f"{sd}/dry.top", "w").write(
            INC + f"[ system ]\ndry {poly} slab\n\n[ molecules ]\n{poly} {NPOLY}\n")
        open(f"{sd}/em.mdp", "w").write(EM)
        run(["gmx", "grompp", "-f", f"{sd}/em.mdp", "-c", f"{sd}/dry.gro",
             "-p", f"{sd}/dry.top", "-o", f"{sd}/em.tpr", "-maxwarn", "5"])
        run(["gmx", "mdrun", "-deffnm", f"{sd}/em",
             "-nb", "gpu", "-ntmpi", "1", "-ntomp", NT])
        # -- compression is TWO stages.  at 1 bar the barostat would need ~40 ns
        #    to collapse a 6x void (C-rescale rate = compressibility x dP / tau-p),
        #    and it would finish "successfully" only half collapsed.  200 bar does
        #    it in ~200 ps; the second stage then finds the 1 bar equilibrium.
        src = f"{sd}/em.gro"
        for tag, ns, T, refp, taup in CMP_STAGES:
            open(f"{sd}/{tag}.mdp", "w").write(
                f"integrator=md\ndt=0.002\nnsteps={int(ns*1000/0.002)}\n"
                f"nstlist=20\nnstxout-compressed=25000\n"
                + COMMON.replace("ref-t=300", f"ref-t={T}") +
                f"pcoupl=C-rescale\npcoupltype=semiisotropic\ntau-p={taup}\n"
                f"ref-p={refp}\ncompressibility=0 4.5e-5\n"
                f"gen-vel=yes\ngen-temp={T}\n")
            run(["gmx", "grompp", "-f", f"{sd}/{tag}.mdp", "-c", src,
                 "-p", f"{sd}/dry.top", "-o", f"{sd}/{tag}.tpr", "-maxwarn", "5"])
            run(["gmx", "mdrun", "-deffnm", f"{sd}/{tag}"] + GPU)
            src = f"{sd}/{tag}.gro"
            log(f"SLAB {poly}: {tag} ({T} K, {ns} ns) -> z = "
                f"{gro_box(f'{sd}/{tag}.gro')[2]:.2f} nm")

        z3   = gro_box(f"{sd}/cmp3.gro")[2]
        zeq  = gro_box(f"{sd}/cmp4.gro")[2]
        drift = abs(zeq - z3) / z3 * 100.0
        log(f"SLAB {poly}: last 1 ns drift = {drift:.1f}% "
            + ("(수렴)" if drift <= 2.0 else "*** 미수렴 — 어닐링을 늘려야 함 ***"))
        if zeq < ZMIN_OK:
            raise RuntimeError(f"{poly}: slab collapsed to z = {zeq:.2f} nm, below the "
                               f"{ZMIN_OK} nm PBC floor (2 x rvdw). check the force field.")
        mass = NPOLY * itp_mass(f"{sd}/{poly}_clean.itp")
        rho  = mass / 6.02214076e23 / (4.0 * 4.0 * zeq * 1e-21)
        log(f"SLAB {poly}: z {hz0:.2f} -> {zeq:.2f} nm, rho = {rho:.2f} g/cm3 "
            f"(amorphous ref {EXP_RHO[poly]:.2f}, {rho/EXP_RHO[poly]*100:.0f}%)")
        open(f"{sd}/density.txt", "w").write(
            f"{poly}\t{NPOLY}\t{hz0:.3f}\t{zeq:.3f}\t{rho:.4f}\t{EXP_RHO[poly]:.2f}\n")

        run(["gmx", "trjconv", "-f", f"{sd}/cmp4.gro", "-s", f"{sd}/cmp4.tpr",
             "-pbc", "mol", "-o", f"{sd}/whole.gro"], inp="0\n")
        gro_shift(f"{sd}/whole.gro", f"{sd}/slab.gro", 4.0, 4.0, SLAB_BOXZ)

    # ------------------------------------------- stage 2: slab + PFAS + solvent
    ztop = gro_maxz(f"{sd}/slab.gro")
    open(f"{wd}/pos.dat", "w").write(f"2.0 2.0 {ztop + 0.7:.3f}\n")
    run(["gmx", "insert-molecules", "-f", f"{sd}/slab.gro",
         "-ci", f"pfas/{pfas}.pdb", "-ip", f"{wd}/pos.dat",
         "-dr", "0.5", "0.5", "0.5", "-try", "500", "-o", f"{wd}/box.gro"])

    apoly = natoms_pdb(f"polymers/{poly}.pdb") * NPOLY
    apfas = natoms_pdb(f"pfas/{pfas}.pdb")
    got = gro_read(f"{wd}/box.gro")[1]
    if got != apoly + apfas:
        raise RuntimeError(f"insert-molecules produced {got} atoms, "
                           f"expected {apoly + apfas} -- PFAS was not placed")

    open(f"{wd}/system.top", "w").write(
        INC + f"[ system ]\n{poly}+{pfas}\n\n"
              f"[ molecules ]\n{poly} {NPOLY}\n{pfas} 1\n")
    open(f"{wd}/em.mdp", "w").write(EM)
    open(f"{wd}/ions.mdp", "w").write(EM)
    run(["gmx", "solvate", "-cp", f"{wd}/box.gro", "-cs", "spc216.gro",
         "-o", f"{wd}/solv.gro", "-p", f"{wd}/system.top"])
    run(["gmx", "grompp", "-f", f"{wd}/ions.mdp", "-c", f"{wd}/solv.gro",
         "-p", f"{wd}/system.top", "-o", f"{wd}/ions.tpr", "-maxwarn", "5"])
    if nc != 0:
        run(["gmx", "genion", "-s", f"{wd}/ions.tpr", "-o", f"{wd}/ions.gro",
             "-p", f"{wd}/system.top", "-pname", "NA", "-nname", "CL",
             "-neutral"], inp="SOL\n")
        cstart = f"{wd}/ions.gro"
    else:
        cstart = f"{wd}/solv.gro"
    run(["gmx", "grompp", "-f", f"{wd}/em.mdp", "-c", cstart,
         "-p", f"{wd}/system.top", "-o", f"{wd}/em.tpr", "-maxwarn", "5"])
    run(["gmx", "mdrun", "-deffnm", f"{wd}/em",
         "-nb", "gpu", "-ntmpi", "1", "-ntomp", NT])

    def fmt(a, b):
        return "".join(f"{i}\n" if (k+1) % 15 == 0 else f"{i} "
                       for k, i in enumerate(range(a, b+1)))
    # run_seed()'s -rerun grompp passes "-n index.ndx".  GROMACS then resolves
    # group names ONLY against this file, so the mdp's tc-grps=System has to be
    # in here too -- otherwise every production run completes and is thrown away
    # at the rerun step.
    ntot = gro_read(f"{wd}/em.gro")[1]
    with open(f"{wd}/index.ndx", "w") as f:
        f.write("[ System ]\n" + fmt(1, ntot) + "\n")
        f.write("[ POL ]\n" + fmt(1, apoly) + "\n")
        f.write("[ PFA ]\n" + fmt(apoly+1, apoly+apfas) + "\n")

    npt = ("integrator=md\ndt=0.002\nnsteps=100000\nnstlist=20\n" + COMMON +
           "pcoupl=C-rescale\npcoupltype=isotropic\ntau-p=2.0\nref-p=1.0\n"
           "compressibility=4.5e-5\ngen-vel=yes\ngen-temp=300\n")
    open(f"{wd}/npt.mdp", "w").write(npt)
    run(["gmx", "grompp", "-f", f"{wd}/npt.mdp", "-c", f"{wd}/em.gro",
         "-p", f"{wd}/system.top", "-o", f"{wd}/npt.tpr", "-maxwarn", "5"])
    run(["gmx", "mdrun", "-deffnm", f"{wd}/npt"] + GPU)
    return apoly, apfas
def run_seed(poly, pfas, wd, seed):
    nsteps=int(PROD_NS*1000/0.002)
    prod=(f"integrator=md\ndt=0.002\nnsteps={nsteps}\nnstlist=40\nnstenergy=5000\n"
          f"nstxout-compressed=5000\n"+COMMON+
          "pcoupl=C-rescale\npcoupltype=isotropic\ntau-p=2.0\nref-p=1.0\ncompressibility=4.5e-5\n"
          f"gen-vel=yes\ngen-temp=300\ngen-seed={seed*1000+7}\n")
    tag=f"{wd}/s{seed}"
    open(f"{tag}_prod.mdp","w").write(prod)
    run(["gmx","grompp","-f",f"{tag}_prod.mdp","-c",f"{wd}/npt.gro","-p",f"{wd}/system.top",
         "-o",f"{tag}_prod.tpr","-maxwarn","5"])
    run(["gmx","mdrun","-deffnm",f"{tag}_prod"]+GPU)
    # rerun with energygrps on CPU for interaction energy
    re=(f"integrator=md\ndt=0.002\nnsteps=0\nnstlist=40\n"+COMMON+
        "energygrps=POL PFA\n")
    open(f"{tag}_re.mdp","w").write(re)
    run(["gmx","grompp","-f",f"{tag}_re.mdp","-c",f"{wd}/npt.gro","-p",f"{wd}/system.top",
         "-n",f"{wd}/index.ndx","-o",f"{tag}_re.tpr","-maxwarn","5"])
    run(["gmx","mdrun","-s",f"{tag}_re.tpr","-rerun",f"{tag}_prod.xtc","-deffnm",f"{tag}_re","-nb","cpu","-ntmpi","1","-ntomp",NT])
    # extract energies (last 4 ns => skip first 2 ns)
    e=run(["gmx","energy","-f",f"{tag}_re.edr","-o",f"{tag}_ie.xvg","-b","2000"],
          inp="Coul-SR:POL-PFA\nLJ-SR:POL-PFA\n\n")
    coul=lj=float("nan")
    for ln in e.stdout.splitlines():
        if ln.startswith("Coul-SR:POL-PFA"): coul=float(ln.split()[1])
        if ln.startswith("LJ-SR:POL-PFA"):  lj=float(ln.split()[1])
    return coul, lj

def main():
    new=not os.path.exists(CSV)
    f=open(CSV,"a",newline=""); w=csv.writer(f)
    if new: w.writerow(["polymer","pfas","seed","coulSR_kJmol","ljSR_kJmol","totalIE_kJmol"]); f.flush()
    done=set()
    if not new:
        for r in csv.DictReader(open(CSV)): done.add((r["polymer"],r["pfas"],r["seed"]))
    total=len(POLYS)*len(PFASS)*len(SEEDS); cnt=0
    log(f"BATCH START: {total} runs ({len(POLYS)}x{len(PFASS)}x{len(SEEDS)})")
    for poly in POLYS:
        for pfas in PFASS:
            wd=f"systems/{poly}_{pfas}"
            need=[s for s in SEEDS if (poly,pfas,str(s)) not in done]
            if not need:
                log(f"skip {poly}+{pfas} (done)"); cnt+=len(SEEDS); continue
            try:
                log(f"BUILD {poly}+{pfas}")
                build(poly,pfas,wd)
            except Exception as ex:
                log(f"BUILD FAIL {poly}+{pfas}: {ex}"); continue
            for s in need:
                cnt=(POLYS.index(poly)*len(PFASS)+list(PFASS).index(pfas))*len(SEEDS)+SEEDS.index(s)+1
                try:
                    t0=time.time()
                    coul,lj=run_seed(poly,pfas,wd,s)
                    tot=coul+lj
                    w.writerow([poly,pfas,s,f"{coul:.3f}",f"{lj:.3f}",f"{tot:.3f}"]); f.flush()
                    log(f"[{cnt}/{total}] {poly}+{pfas} s{s}: IE={tot:.1f} kJ/mol ({time.time()-t0:.0f}s)")
                except Exception as ex:
                    log(f"RUN FAIL {poly}+{pfas} s{s}: {ex}")
            # cleanup big files to save disk
            for pat in ("*_prod.trr","*_re.trr","#*#"):
                for fn in glob.glob(f"{wd}/{pat}"):
                    try: os.remove(fn)
                    except OSError: pass
    log("BATCH DONE")
    f.close()

if __name__=="__main__":
    main()
