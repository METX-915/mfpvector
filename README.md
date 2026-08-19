# mfpvector

Data and analysis code for the manuscript:

> **Microfluoroplastics as a Non-Extractable Organofluorine Pool: A Measurement Framework for Environmental Fluorine Mass Balance**
> HanYeol Nam, Yoon Jeong Jang, Ku Kang
> CBRN Defense Research Institute, Republic of Korea Ministry of National Defense
> Submitted to *Environmental Science: Processes & Impacts* (RSC)

Canonical URL: <https://github.com/METX-915/mfpvector>

This repository archives the analysis scripts, per-replicate outputs, and figure
data needed to reproduce the numbers in the manuscript and its Electronic
Supplementary Information (ESI). Scope and known gaps are stated explicitly below.

---

## What this study does

Aggregate fluorine metrics (TF, TOF, EOF, AOF) miss insoluble fluoropolymer
particles, so microfluoroplastic (MFP) fluorine belongs to the *non-extractable*
organofluorine pool rather than to the unexplained extractable gap. The repository
supports three quantitative components:

1. **Mass-balance bounding** — a Monte Carlo scenario (2×10⁵ draws, fixed seed) placing
   MFP-derived sediment fluorine relative to a measured EOF reference, with global
   sensitivity analysis.
2. **Adsorption interaction-energy screen** — GROMACS MD of four PFAS-related
   compounds on PTFE / PVDF / ECTFE and a PE baseline (79 trajectories), plus two
   umbrella-sampling PMF free-energy anchors for PFOA.
3. **Evidence-tier prioritization** — an ordinal PPARα receptor tier × normalized
   adsorption, reported as a decision matrix (not a calibrated risk index).

---

## Repository layout

```
mfpvector/
├── code/
│   ├── rq1_massbalance.py      # Monte Carlo mass-balance bound (seed 20260710) + Spearman sensitivity
│   ├── analyze_rq2.py          # adsorption interaction-energy analysis (means, CI, effect sizes)
│   ├── make_figures.py         # regenerates the adsorption 4-panel figure (manuscript Fig 3) from data/
│   ├── make_fig4.py            # regenerates the evidence-tier figure (manuscript Fig 4) from data/
│   ├── batch_run.py            # MD driver that produced data/raw_results.csv
│   └── emit_md_inputs.py       # re-emits every .mdp / packmol input from batch_run.py
├── data/
│   ├── raw_results.csv         # per-replicate short-range interaction energies (79 trajectories)
│   ├── results_table.csv       # per-cell mean IE with bootstrap 95% CI
│   ├── experiment_summary.json # per-cell means/CI/SD + effect sizes vs PE
│   ├── mfp_vi.json             # normalized adsorption (Â) and receptor evidence tier (Ĥ) → tier matrix
│   ├── rq1_summary.json        # mass-balance percentiles + sensitivity (Spearman ρ)
│   ├── rq3_docking.json        # PPARα (PDB 1K7L) smina docking scores (consistency check only)
│   ├── pmf_result.json         # PVDF+PFOA umbrella-sampling PMF anchor (ΔG = −13.8 kJ/mol)
│   └── tfeof_meta.json         # DOI-verified TF/TOF–EOF/AOF literature set
├── figures/                    # created by make_figures.py (not tracked)
├── md_inputs/                  # created by emit_md_inputs.py (not tracked)
└── out/                        # created by analyze_rq2.py (not tracked)
```

Figure ↔ file mapping (as compiled in the manuscript): **Fig 1** framework schematic ·
**Fig 2** mass-balance bounding · **Fig 3** adsorption screen (produced by
`make_figures.py`) · **Fig 4** evidence-tier prioritization.

---

## Reproducing the results

```bash
pip install -r requirements.txt

# Mass-balance bound + sensitivity (prints percentiles; matches data/rq1_summary.json)
python code/rq1_massbalance.py

# Adsorption means, bootstrap 95% CI, effect sizes
# writes out/results_table.csv, byte-identical to data/results_table.csv
python code/analyze_rq2.py

# Regenerate the adsorption figure -> figures/Fig2.pdf (vector) and figures/Fig2.png
python code/make_figures.py
python code/make_fig4.py

# Re-emit the GROMACS .mdp and packmol inputs used for the production runs.
# Parameters are read out of code/batch_run.py with `ast`, not retyped, so the
# emitted files are the ones the runs actually used.
python code/emit_md_inputs.py
```

All scripts resolve their paths relative to the repository root and run from a fresh
clone with no configuration. Random draws use fixed seeds (`20260710` for the Monte
Carlo, `7` for the bootstrap), so the numeric outputs are bit-reproducible.

---

## Data provenance and integrity notes

- Interaction energies are used only as **ordinal ranking proxies**; they are not free
  energies. Absolute PFOA adsorption free energies come from umbrella-sampling PMFs.
- **Every panel of the adsorption figure is generated from the files in `data/`**, with one
  stated exception: the ECTFE+PFOA PMF free energy (ΔG = −28.5 kJ/mol) plotted in panel (d)
  is an author-reported result (ESI Table S6) whose window-level output is not yet archived
  (see Known gaps). The PVDF+PFOA anchor is in `pmf_result.json`.
- **Zero-interaction replicates.** Short-range Coulomb and Lennard-Jones terms are
  truncated at a 1.0 nm cutoff, so a replicate in which the compound leaves the surface
  and stays beyond that cutoff for the whole averaging window returns exactly
  `0.000`. Four of the 79 replicates in `raw_results.csv` are of this type
  (`ECTFE/PFBS/4`, `PE/PFOA/4`, `PE/GenX/1`, `PE/GenX/5`) and are **retained** in the
  cell means, so those means are ensemble averages over bound and spontaneously
  desorbed states. ESI Table S5 reports the full analysis with these four replicates
  removed; no conclusion in the manuscript depends on which treatment is used.
- **PE–6:2 FTOH has four replicates, not five**; the seed-5 trajectory is absent from
  the archived output, which is why the screen contains 79 rather than 80 trajectories.
- `analyze_rq2.py` lists **PP** among the polymers because it was in the original run
  plan, but PP produced no completed trajectories and appears nowhere in the data or
  the manuscript.
- The only statistically significant polymer-level contrast in the analysis as reported
  is **PVDF > PE** (p = 0.004); all other polymer contrasts are descriptive.
- The MFP-fraction prior (2–8 % of microplastic count, `rng.uniform(0.02, 0.08)` in
  `rq1_massbalance.py`) derives from a single LDIR reanalysis and is treated
  as a **bounding prior**, not a universal abundance.

### Not archived (available from the corresponding author on request)


- **ECTFE+PFOA PMF window data (ΔG = −28.5 kJ/mol) is missing**; `pmf_result.json`
  holds the PVDF+PFOA anchor only.
- The acpype/GAFF2 topologies (`*_GMX.itp`), the polymer and PFAS `.pdb` structures, and
  `index.ndx` are **not archived**. They cannot be regenerated here: AM1-BCC charges depend
  on the exact input geometry, so re-deriving them would not give the files the runs used.
- The umbrella-sampling window definitions and the `gmx wham` outputs (`profile.xvg`,
  `bsResult.xvg`) for both PMF anchors are **not archived**.
- The short-timescale isotonic-saline stability check reported in ESI Section S5 / Table S7
  has no archived inputs or output.

The manuscript's Data Availability statement lists these items as available from the
corresponding author on request, not as archived here.

---

## Release status

Public release of this repository has been cleared; no institutional security review is
outstanding. This repository is the archive referenced by the manuscript's Data
Availability statement. A DOI can additionally be minted by enabling the Zenodo–GitHub
integration and publishing a tagged release.

## License

Code: MIT (see `LICENSE`). Data and figures: CC-BY-4.0.

## Citing this repository

If you use these data or scripts, please cite the manuscript above. A versioned
snapshot with a DOI can be created by enabling the Zenodo–GitHub integration and
publishing a release.

## Contact

Corresponding author: Ku Kang — bisu9082@gmail.com
Repository maintained by: HanYeol Nam (@METX-915) — namhy915@gmail.com
