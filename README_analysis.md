# Analysis pipeline (JHM submission, 2026-10)

All numbers in the manuscript and Supplementary material are generated from `data/raw_results.csv`
(80 rows: polymer, compound, seed, Coulomb, Lennard-Jones, total short-range interaction energy, kJ/mol)
and from the archived WHAM outputs in `data/pmf/`. Run from the repository root:

    python code/rq1_massbalance.py                                   # Monte Carlo mass-balance bound (seed 20260710)
    python code/analyze_screen.py data/raw_results.csv out           # stats -> out/experiment_summary.json, results_table.csv
    python code/make_fig3.py out out/fig                             # evidence tier (out/tier.json) + Fig 3
    python code/make_fig2.py out data/raw_results.csv out/fig data/pmf/pmf_summary.json   # Fig 2
    python code/make_ga.py out out/fig                               # graphical abstract
    python code/make_numbers_tex.py out data/pmf/pmf_summary.json    # numbers.tex macros
    python code/make_si_tables.py out data/raw_results.csv data/pmf/pmf_summary.json      # SI tables
    python code/make_figS1.py out/fig PTFE=data/pmf/PTFE_prof.xvg:data/pmf/PTFE_bs.xvg:1.147:3.674 \
        PVDF=data/pmf/PVDF_prof.xvg:data/pmf/PVDF_bs.xvg:1.512:3.642 \
        ECTFE=data/pmf/ECTFE_prof.xvg:data/pmf/ECTFE_bs.xvg:1.653:3.966 \
        PE=data/pmf/PE_prof.xvg:data/pmf/PE_bs.xvg:1.408:3.677                             # Fig S1
    python code/emit_md_inputs.py                                    # .mdp / packmol inputs of the production runs

`data/pmf/pmf_summary.json` is assembled from the `report_*` files in `data/pmf/` (values: dG_ads, error, windows,
analysed segment); `code/pmf_from_report.py` parses a concatenation of those reports.
`code/batch_run.py` is the two-stage slab-construction and production driver; `code/pmf/rerun_pmf.sh` and
`finish_pmf.sh` with `pull_smd.mdp` / `pull_umbrella.mdp` are the umbrella-sampling drivers.
Statistics: compound-blocked paired t (df = 3), Friedman, replicate-level Welch t / Mann-Whitney, bootstrap CIs (10^4).
Superseded scripts removed in v1.0.1: analyze_rq2.py, make_figures.py, make_fig4.py, data/mfp_vi.json.
