#!/usr/bin/env python3
"""Charge-model sensitivity of the short-range Coulomb term (SI Section on charge models).

For every screen cell under systems/<P>_<X>/ that still has production trajectories (s<k>_prod.xtc), the polymer and
compound partial charges are replaced by an alternative charge model (Gasteiger-Marsili or MMFF94, from charges_alt.json,
in the atom order of the GROMACS topologies) and the trajectory is re-evaluated with the same -rerun settings as the screen
(energygrps POL PFA, final 4 ns). Only the charges change; the configurations are those sampled with AM1-BCC charges, so the
result is a single-point (perturbative) test of how the sign and size of the Coulomb term depend on the charge model.

Usage (in ~/real2AB, conda env md):  python charge_sensitivity.py      (resumable; output results/charge_sensitivity.csv)
"""
import os, re, csv, json, glob, subprocess, time

HERE = os.path.dirname(os.path.abspath(__file__))
CH = json.load(open(os.path.join(HERE, "charges_alt.json")))
OUT = os.path.join(HERE, "results", "charge_sensitivity.csv")
NT = str(min(8, os.cpu_count() or 8))
MODELS = ["gasteiger", "mmff94"]


def run(cmd, inp=None, cwd=None):
    r = subprocess.run(cmd, input=inp, capture_output=True, text=True, cwd=cwd)
    if r.returncode != 0:
        raise RuntimeError(" ".join(cmd) + "\n" + r.stderr[-1500:])
    return r


def recharge(src, dst, q):
    """Copy a *_clean.itp replacing the charge column of [ atoms ] by q (in order)."""
    out, inat, k = [], False, 0
    for ln in open(src):
        s = ln.strip()
        if s.startswith("["):
            inat = s.replace(" ", "").lower().startswith("[atoms]")
            out.append(ln); continue
        body = ln.split(";")[0]
        f = body.split()
        if inat and len(f) >= 8:
            f[6] = f"{q[k]:.6f}"; k += 1
            out.append("  ".join(f) + "\n"); continue
        out.append(ln)
    assert k == len(q), (src, k, len(q))
    open(dst, "w").write("".join(out))


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    done = set()
    if os.path.exists(OUT):
        for r in csv.DictReader(open(OUT)):
            done.add((r["system"], r["seed"], r["model"]))
    new = not os.path.exists(OUT)
    fo = open(OUT, "a", newline=""); w = csv.writer(fo)
    if new:
        w.writerow(["system", "polymer", "pfas", "seed", "model", "coulSR_kJmol", "ljSR_kJmol"]); fo.flush()
    for wd in sorted(glob.glob(os.path.join(HERE, "systems", "*_*"))):
        sysname = os.path.basename(wd); poly, pfas = sysname.split("_", 1)
        if poly not in CH["gasteiger"] or pfas not in CH["gasteiger"]:
            print("skip (no charges)", sysname); continue
        xtcs = sorted(glob.glob(os.path.join(wd, "s*_prod.xtc")))
        if not xtcs:
            print("skip (no trajectories)", sysname); continue
        top = open(os.path.join(wd, "system.top")).read()
        for model in MODELS:
            recharge(os.path.join(wd, f"{poly}_clean.itp"), os.path.join(wd, f"{poly}_{model}.itp"), CH[model][poly])
            recharge(os.path.join(wd, f"{pfas}_clean.itp"), os.path.join(wd, f"{pfas}_{model}.itp"), CH[model][pfas])
            t2 = top.replace(f'"{poly}_clean.itp"', f'"{poly}_{model}.itp"').replace(f'"{pfas}_clean.itp"', f'"{pfas}_{model}.itp"')
            open(os.path.join(wd, f"system_{model}.top"), "w").write(t2)
            for x in xtcs:
                seed = re.search(r"s(\d+)_prod\.xtc$", x).group(1)
                if (sysname, seed, model) in done:
                    continue
                t0 = time.time(); tag = f"s{seed}_re_{model}"
                run(["gmx", "grompp", "-f", f"s{seed}_re.mdp", "-c", "npt.gro", "-p", f"system_{model}.top",
                     "-n", "index.ndx", "-o", f"{tag}.tpr", "-maxwarn", "5"], cwd=wd)
                run(["gmx", "mdrun", "-s", f"{tag}.tpr", "-rerun", f"s{seed}_prod.xtc", "-deffnm", tag,
                     "-nb", "cpu", "-ntmpi", "1", "-ntomp", NT], cwd=wd)
                e = run(["gmx", "energy", "-f", f"{tag}.edr", "-o", f"{tag}_ie.xvg", "-b", "2000"],
                        inp="Coul-SR:POL-PFA\nLJ-SR:POL-PFA\n\n", cwd=wd)
                coul = lj = float("nan")
                for ln in e.stdout.splitlines():
                    if ln.startswith("Coul-SR:POL-PFA"): coul = float(ln.split()[1])
                    if ln.startswith("LJ-SR:POL-PFA"): lj = float(ln.split()[1])
                for ext in (".trr", ".gro", ".log", ".edr"):
                    p = os.path.join(wd, tag + ext)
                    if ext in (".trr", ".gro") and os.path.exists(p): os.remove(p)
                w.writerow([sysname, poly, pfas, seed, model, f"{coul:.3f}", f"{lj:.3f}"]); fo.flush()
                print(f"[{time.strftime('%m-%d %H:%M:%S')}] {sysname} s{seed} {model}: Coul {coul:.2f} LJ {lj:.2f} ({time.time()-t0:.0f}s)", flush=True)
    print("CHARGE SENSITIVITY DONE", flush=True)


if __name__ == "__main__":
    main()
