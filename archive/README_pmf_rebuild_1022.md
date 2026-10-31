# mfpvector PTFE/PVDF + PFOA umbrella-sampling PMF — input & settings archive

Created 2026-10-22 on the recompute PC (RTX 3090, WSL2 Ubuntu, GROMACS 2026.3).
Final run: kit v2.1 (2026-10-14 ~ 10-16) + finish step v2.2 (2026-10-19).

## Contents
- `kit/`            : recompute kit scripts (rerun_pmf.sh, finish_pmf.sh, pmf_scripts/, build_kit/) — v2.2
- `work_pmf/`       : script copies actually executed + all stage logs (rerun.log, finish.log, build/eq/rb logs) and block reports
- `mfpvector/`      : run tree (topologies *.acpype, polymers/pfas structures, patched code/batch_run.py,
                      slabs/ (dry-slab build), systems/ (screened systems), pmf/ (equilibration, SMD, umbrella windows, WHAM), pmf_final/)
- `ENVIRONMENT.txt` : OS / GPU / GROMACS build / conda package lists
- `MANIFEST.tsv`    : every file with size (bytes)

Per window (`pmf/<sys>_PFOA/win2/{n,i,j}NN/`): umb.mdp (where changed), umb.tpr (full run input), umb.gro (final), logs.
Pull coordinate/force time series are in the companion archive `pmf_pull_1022.tgz`.

## Excluded (kept on the recompute PC / raw archive pmf_raw_fin_1019.tgz)
trajectories (*.xtc, *.trr), checkpoints (*.cpt), energy files (*.edr), GROMACS backups (#*#),
SMD frame dumps (smd2/frame*.gro), recipient cover letters (kit *.md).

## Note on modifications vs. original kit
see work_pmf/rerun.log (NOTE lines) and kit/rerun_pmf.sh header:
acpype time limit 24 h, [ System ] index group, com_dist.py PBC-aware COM, wham2 mkdir fix.
