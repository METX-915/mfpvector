# Analysis pipeline (JHM revision, 2026-09)

All numbers in the manuscript and Supplementary material are generated from `data/raw_results.csv`
(80 rows: polymer, compound, seed, Coulomb, Lennard-Jones, total short-range interaction energy, kJ/mol)
and, for the free-energy anchors, from `pmf/*/wham/profile.xvg`.

    python code/analyze_screen.py  data/raw_results.csv out          # stats -> out/experiment_summary.json, results_table.csv
    python code/make_fig3.py       out out/fig                       # evidence-tier matrix (writes out/tier.json, Fig3)
    python code/make_fig2.py       out data/raw_results.csv out/fig [out/pmf.json]   # adsorption figure
    python code/make_ga.py         out out/fig                       # graphical abstract
    python code/make_numbers_tex.py out [out/pmf.json]               # numbers.tex macros used by manuscript.tex / SI.tex
    python code/make_si_tables.py  out data/raw_results.csv [out/pmf.json]   # SI tables
    python code/pmf_from_report.py pmf/reports.txt out/pmf.json      # parse report_pmf.py output for the four PFOA systems
    python code/make_figS1.py      out/fig pmf/PTFE_PFOA/wham/profile.xvg pmf/PVDF_PFOA/wham/profile.xvg pmf/ECTFE_PFOA/wham/profile.xvg pmf/PE_PFOA/wham/profile.xvg

`code/batch_run.py` is the slab-construction and production driver (two-stage packmol -> compress/anneal -> solvate -> 5 x 6 ns);
`pmf/run_pmf.sh` with `pull_smd.mdp` / `pull_umbrella.mdp` is the umbrella-sampling driver.
Statistics: compound-blocked paired t (df = 3), Friedman, replicate-level Welch t / Mann-Whitney, bootstrap CIs (10^4).
