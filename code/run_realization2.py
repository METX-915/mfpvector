#!/usr/bin/env python3
"""Second slab realisation of the PFOA interaction-energy screen (PTFE, PVDF).

Protocol is identical to code/batch_run.py::run_seed of the mfpvector repository
(release v1.0.x): from the archived 200 ps NPT cell of the rebuilt slab (npt.gro),
five production replicates of 6 ns (2 fs, V-rescale 300 K, C-rescale 1 bar isotropic,
PME, 1.0 nm cutoffs, LINCS h-bonds) differing only in gen-seed = seed*1000+7, then a
CPU -rerun with energygrps=POL PFA and gmx energy over the final 4 ns (-b 2000).

Usage (inside the conda env that has gmx, e.g. `conda activate md`):
    python run_realization2.py            # both systems, seeds 1-5
    python run_realization2.py PTFE_PFOA  # one system
Resumable: finished seeds (present in results_realization2.csv) are skipped.
Outputs: results_realization2.csv, <SYS>/s<k>_ie.xvg, <SYS>/s<k>_re.edr, <SYS>/s<k>_prod.log
"""
import os, sys, csv, subprocess, time, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
SYSTEMS = sys.argv[1:] or ["PTFE_PFOA", "PVDF_PFOA"]
SEEDS = [1, 2, 3, 4, 5]
PROD_NS = 6
NT = str(min(8, os.cpu_count() or 8))
CSV = os.path.join(HERE, "results_realization2.csv")
LOG = os.path.join(HERE, "run_realization2.log")

COMMON = ("cutoff-scheme=Verlet\ncoulombtype=PME\nrcoulomb=1.0\nrvdw=1.0\npbc=xyz\n"
          "constraints=h-bonds\nconstraint-algorithm=lincs\ntcoupl=V-rescale\n"
          "tc-grps=System\ntau-t=0.1\nref-t=300\n")
GPU = ["-nb", "gpu", "-pme", "gpu", "-bonded", "gpu", "-update", "gpu",
       "-ntmpi", "1", "-ntomp", NT, "-pin", "on"]


def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    open(LOG, "a").write(line + "\n")


def run(cmd, inp=None, cwd=None):
    r = subprocess.run(cmd, input=inp, capture_output=True, text=True, cwd=cwd)
    if r.returncode != 0:
        open(LOG, "a").write(r.stdout + "\n" + r.stderr + "\n")
        raise RuntimeError(f"command failed: {' '.join(cmd)}\n{r.stderr[-2000:]}")
    return r


def run_seed(wd, seed):
    nsteps = int(PROD_NS * 1000 / 0.002)
    tag = f"s{seed}"
    prod = (f"integrator=md\ndt=0.002\nnsteps={nsteps}\nnstlist=40\nnstenergy=5000\n"
            f"nstxout-compressed=5000\n" + COMMON +
            "pcoupl=C-rescale\npcoupltype=isotropic\ntau-p=2.0\nref-p=1.0\ncompressibility=4.5e-5\n"
            f"gen-vel=yes\ngen-temp=300\ngen-seed={seed*1000+7}\n")
    open(os.path.join(wd, f"{tag}_prod.mdp"), "w").write(prod)
    run(["gmx", "grompp", "-f", f"{tag}_prod.mdp", "-c", "npt.gro", "-p", "system.top",
         "-o", f"{tag}_prod.tpr", "-maxwarn", "5"], cwd=wd)
    try:
        run(["gmx", "mdrun", "-deffnm", f"{tag}_prod"] + GPU, cwd=wd)
    except RuntimeError as e:
        log(f"GPU mdrun failed for {wd}/{tag}, retrying with -nb gpu only: {str(e)[-300:]}")
        run(["gmx", "mdrun", "-deffnm", f"{tag}_prod", "-nb", "gpu", "-ntmpi", "1", "-ntomp", NT], cwd=wd)
    re_mdp = (f"integrator=md\ndt=0.002\nnsteps=0\nnstlist=40\n" + COMMON + "energygrps=POL PFA\n")
    open(os.path.join(wd, f"{tag}_re.mdp"), "w").write(re_mdp)
    run(["gmx", "grompp", "-f", f"{tag}_re.mdp", "-c", "npt.gro", "-p", "system.top",
         "-n", "index.ndx", "-o", f"{tag}_re.tpr", "-maxwarn", "5"], cwd=wd)
    run(["gmx", "mdrun", "-s", f"{tag}_re.tpr", "-rerun", f"{tag}_prod.xtc", "-deffnm", f"{tag}_re",
         "-nb", "cpu", "-ntmpi", "1", "-ntomp", NT], cwd=wd)
    e = run(["gmx", "energy", "-f", f"{tag}_re.edr", "-o", f"{tag}_ie.xvg", "-b", "2000"],
            inp="Coul-SR:POL-PFA\nLJ-SR:POL-PFA\n\n", cwd=wd)
    coul = lj = float("nan")
    for ln in e.stdout.splitlines():
        if ln.startswith("Coul-SR:POL-PFA"):
            coul = float(ln.split()[1])
        if ln.startswith("LJ-SR:POL-PFA"):
            lj = float(ln.split()[1])
    # keep the kit small: drop the trajectory and checkpoints once energies are extracted
    for ext in ("_prod.xtc", "_prod.trr", "_prod.cpt", "_prod_prev.cpt", "_prod.edr", "_re.trr", "_re.gro"):
        p = os.path.join(wd, tag + ext)
        if os.path.exists(p):
            os.remove(p)
    return coul, lj


def main():
    if shutil.which("gmx") is None:
        sys.exit("gmx not found on PATH. Activate the GROMACS conda environment first (e.g. `conda activate md`).")
    done = set()
    new = not os.path.exists(CSV)
    if not new:
        for r in csv.DictReader(open(CSV)):
            done.add((r["system"], r["seed"]))
    f = open(CSV, "a", newline="")
    w = csv.writer(f)
    if new:
        w.writerow(["system", "polymer", "pfas", "realization", "seed", "coulSR_kJmol", "ljSR_kJmol", "totalIE_kJmol"])
        f.flush()
    log(f"START realization-2 screen: {SYSTEMS} x seeds {SEEDS}, PROD_NS={PROD_NS}, ntomp={NT}")
    for s in SYSTEMS:
        wd = os.path.join(HERE, s)
        if not os.path.isfile(os.path.join(wd, "npt.gro")):
            log(f"skip {s}: npt.gro missing")
            continue
        poly, pfas = s.split("_")
        for seed in SEEDS:
            if (s, str(seed)) in done:
                log(f"skip {s} seed {seed} (done)")
                continue
            t0 = time.time()
            log(f"RUN {s} seed {seed}")
            try:
                coul, lj = run_seed(wd, seed)
            except Exception as ex:
                log(f"FAIL {s} seed {seed}: {ex}")
                continue
            w.writerow([s, poly, pfas, 2, seed, f"{coul:.3f}", f"{lj:.3f}", f"{coul+lj:.3f}"])
            f.flush()
            log(f"DONE {s} seed {seed}: Coul {coul:.2f} LJ {lj:.2f} total {coul+lj:.2f} kJ/mol ({(time.time()-t0)/60:.1f} min)")
    log("ALL DONE")


if __name__ == "__main__":
    main()
