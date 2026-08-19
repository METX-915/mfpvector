#!/usr/bin/env python3
"""
mfpvector batch MD driver — RQ2 relative PFAS adsorption on (fluoro)polymer slabs.
20 systems (5 polymers x 4 PFAS, pristine) x 5 seeds x 6 ns GPU production.
Metric: polymer-PFAS short-range interaction energy (LJ+Coul, via CPU -rerun),
        + min-distance/contact fraction. Absolute values NOT claimed -> relative ranking.
Run in background:  nohup python scripts/batch_run.py > results/batch.out 2>&1 &
"""
import os, subprocess, sys, time, csv, glob

os.chdir(os.path.expanduser("~/mfpvector"))
R="results"; os.makedirs(R, exist_ok=True)
POLYS=["PTFE","PVDF","ECTFE","PE","PP"]
PFASS={"PFOA":-1,"PFBS":-1,"GenX":-1,"FTOH62":0}
SEEDS=[1,2,3,4,5]
NPOLY=40; PROD_NS=6; BXY=40.0
CSV=os.path.join(R,"raw_results.csv")
LOG=open(os.path.join(R,"progress.log"),"a")

def log(m):
    s=f"[{time.strftime('%m-%d %H:%M:%S')}] {m}"
    print(s,flush=True); LOG.write(s+"\n"); LOG.flush()

def run(cmd, inp=None, shell=False, quiet=True):
    r=subprocess.run(cmd, input=inp, capture_output=True, text=True, shell=shell)
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
GPU=["-nb","gpu","-pme","gpu","-bonded","gpu","-update","gpu","-ntmpi","1","-ntomp","8","-pin","on"]

def build(poly, pfas, wd):
    os.makedirs(wd, exist_ok=True)
    nc=PFASS[pfas]
    poly_itp=open(f"{poly}.acpype/{poly}_GMX.itp").read()
    pfas_itp=open(f"{pfas}.acpype/{pfas}_GMX.itp").read()
    at={}
    for t in (poly_itp,pfas_itp):
        for ln in extract_at(t): at[ln.split()[0]]=ln
    with open(f"{wd}/atomtypes.itp","w") as f:
        f.write("[ atomtypes ]\n;name bond_type mass charge ptype sigma epsilon\n")
        for v in at.values(): f.write(v+"\n")
        for e in EXTRA: f.write(e+"\n")
    open(f"{wd}/{poly}_clean.itp","w").write(rename_mol(strip_at(poly_itp),poly))
    open(f"{wd}/{pfas}_clean.itp","w").write(rename_mol(strip_at(pfas_itp),pfas))
    open(f"{wd}/water.itp","w").write(WATER); open(f"{wd}/ions.itp","w").write(IONS)
    # packmol: PFAS placed just above slab surface (z 28-42) to sample adsorption
    pack=f"""tolerance 2.0
filetype pdb
output {wd}/dry.pdb
structure polymers/{poly}.pdb
  number {NPOLY}
  inside box 2. 2. 2. {BXY-2} {BXY-2} 28.
end structure
structure pfas/{pfas}.pdb
  number 1
  inside box 8. 8. 29. {BXY-8} {BXY-8} 40.
end structure
"""
    open(f"{wd}/pack.inp","w").write(pack)
    run(f"packmol < {wd}/pack.inp", shell=True)
    run(["gmx","editconf","-f",f"{wd}/dry.pdb","-o",f"{wd}/box.gro","-box","4.0","4.0","8.0"])
    top=(f'[ defaults ]\n1 2 yes 0.5 0.8333333333\n\n'
         f'#include "atomtypes.itp"\n#include "{poly}_clean.itp"\n#include "{pfas}_clean.itp"\n'
         f'#include "water.itp"\n#include "ions.itp"\n\n[ system ]\n{poly}+{pfas}\n\n'
         f'[ molecules ]\n{poly} {NPOLY}\n{pfas} 1\n')
    open(f"{wd}/system.top","w").write(top)
    emmdp="integrator=steep\nnsteps=3000\nemtol=1000\nemstep=0.01\nnstlist=10\n"+COMMON.split('tcoupl')[0]
    open(f"{wd}/em.mdp","w").write(emmdp); open(f"{wd}/ions.mdp","w").write(emmdp)
    run(["gmx","solvate","-cp",f"{wd}/box.gro","-cs","spc216.gro","-o",f"{wd}/solv.gro","-p",f"{wd}/system.top"])
    run(["gmx","grompp","-f",f"{wd}/ions.mdp","-c",f"{wd}/solv.gro","-p",f"{wd}/system.top","-o",f"{wd}/ions.tpr","-maxwarn","5"])
    if nc!=0:
        run(["gmx","genion","-s",f"{wd}/ions.tpr","-o",f"{wd}/ions.gro","-p",f"{wd}/system.top",
             "-pname","NA","-nname","CL","-neutral"], inp="SOL\n")
        cstart=f"{wd}/ions.gro"
    else:
        cstart=f"{wd}/solv.gro"
    run(["gmx","grompp","-f",f"{wd}/em.mdp","-c",cstart,"-p",f"{wd}/system.top","-o",f"{wd}/em.tpr","-maxwarn","5"])
    run(["gmx","mdrun","-deffnm",f"{wd}/em","-nb","gpu","-ntmpi","1","-ntomp","8"])
    # index by atom-number ranges (poly first, pfas next)
    apoly=natoms_pdb(f"polymers/{poly}.pdb")*NPOLY
    apfas=natoms_pdb(f"pfas/{pfas}.pdb")
    def fmt(a,b):
        return "".join(f"{i}\n" if (k+1)%15==0 else f"{i} " for k,i in enumerate(range(a,b+1)))
    with open(f"{wd}/index.ndx","w") as f:
        f.write("[ POL ]\n"+fmt(1,apoly)+"\n[ PFA ]\n"+fmt(apoly+1,apoly+apfas)+"\n")
    # NPT once (position from em)
    npt=("integrator=md\ndt=0.002\nnsteps=100000\nnstlist=20\n"+COMMON+
         "pcoupl=C-rescale\npcoupltype=isotropic\ntau-p=2.0\nref-p=1.0\ncompressibility=4.5e-5\ngen-vel=yes\ngen-temp=300\n")
    open(f"{wd}/npt.mdp","w").write(npt)
    run(["gmx","grompp","-f",f"{wd}/npt.mdp","-c",f"{wd}/em.gro","-p",f"{wd}/system.top","-o",f"{wd}/npt.tpr","-maxwarn","5"])
    run(["gmx","mdrun","-deffnm",f"{wd}/npt"]+GPU)
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
    run(["gmx","mdrun","-s",f"{tag}_re.tpr","-rerun",f"{tag}_prod.xtc","-deffnm",f"{tag}_re","-nb","cpu","-ntmpi","1","-ntomp","8"])
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
