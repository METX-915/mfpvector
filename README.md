# mfpvector

Data and analysis code for the manuscript

> **Fluorine content does not set simulated PFAS sorption to fluoropolymer microplastics: free-energy anchors and a bound on an unquantified organofluorine pool**
> HanYeol Nam, Yoon Jeong Jang, Jin Yoo, Doo-Hee Lee, Ku Kang — CBRN Defense Research Institute, Republic of Korea Ministry of National Defense
> Submitted to *Journal of Hazardous Materials* (2026-10)

Canonical URL: <https://github.com/METX-915/mfpvector>. Archived at Zenodo — concept DOI (all versions)
[10.5281/zenodo.23838772](https://doi.org/10.5281/zenodo.23838772); v1.0.0 = 10.5281/zenodo.23838773 (2026-10-31), v1.0.1 = see Zenodo.

## What the repository supports

1. **Adsorption interaction-energy screen** — GROMACS MD (GAFF2/AM1-BCC, TIP3P, 300 K) of PFOA, PFBS, HFPO-DA and
   6:2 FTOH on amorphous PTFE, PVDF, ECTFE and polyethylene slabs: 4 × 4 × 5 = 80 trajectories, with Coulomb /
   Lennard-Jones decomposition (`data/raw_results.csv`). All 16 cells were repeated on a second, independently
   built slab per polymer (`data/realization2/`), and the Coulomb term was re-evaluated with Gasteiger-Marsili and MMFF94
   charges (`data/charge_model/`).
2. **Umbrella-sampling free-energy anchors** — PMFs of PFOA on all four polymers with one protocol (`data/pmf/`):
   PTFE −30.7 ± 3.5, PVDF −26.6 ± 3.2, ECTFE −32.2 ± 2.3, PE −31.4 ± 3.9 kJ mol⁻¹.
3. **Fluorine mass-balance bound** — Monte Carlo scenario (2 × 10⁵ draws, seed 20260710) placing MFP-derived sediment
   fluorine against a measured EOF reference, with Spearman sensitivity (`code/rq1_massbalance.py`).
4. **Evidence-tier prioritization** — normalized adsorption × ordinal PPARα docking tier (`data/tier.json`,
   `data/rq3_docking.json`); Supplementary Fig. S2.

## Layout

```
code/
  batch_run.py            slab construction (two-stage packmol -> compress/anneal) + screen production driver
  build/                  original builder pieces: patch_build.py, compact_conf.py (ETKDG, seed 20260710),
                          make_structures.py, build_one.py, polymers/*.pdb, pfas/*.pdb
  emit_md_inputs.py       re-emits every .mdp / packmol input from batch_run.py
  analyze_screen.py       screen statistics -> experiment_summary.json, results_table.csv
  make_fig1.py, make_fig2.py, make_fig4.py (Fig 3, mass balance), make_fig3.py (SI Fig S2, tier), make_ga.py,
                          make_figS1.py, make_numbers_tex.py, make_si_tables.py
  sasa.py, ksurf.py       surface composition of the slabs (SI Table S16); surface partition coefficients from the PMFs
  run_realization2.py, realization2_setAB/, realization2_setCD/   second slab realisation drivers (SI Section S11)
  charge_model/           alternative charges (RDKit) and -rerun re-evaluation; SI Section S12 statistics and tables
  pmf_from_report.py      parses report_pmf.py output -> pmf summary json
  rq1_massbalance.py      Monte Carlo mass-balance bound
  pmf/                    umbrella-sampling protocol: rerun_pmf.sh, finish_pmf.sh, equil_cell.sh, shift_pfa.py,
                          pull_smd.mdp, pull_umbrella.mdp, rebuild_windows.sh, extend_pmf.sh, infill.py,
                          report_pmf.py, surf_z.py, com_dist.py, central_atom.py
  README_analysis.md      command sequence that regenerates every number in the manuscript and SI
topologies/               GAFF2/AM1-BCC topologies of the rebuilt PTFE, PVDF and PFOA (acpype: *_GMX.itp/.gro/.top, input pdb, mol2, frcmod)
slabs/                    PTFE and PVDF dry slabs: packmol input, EM and four compression-stage .mdp, dry.gro/slab.gro, density, check
md_inputs/                PTFE_PFOA and PVDF_PFOA screen cells: system.top, index.ndx, npt.gro, em/ions/npt .mdp, clean itp files
archive/                  complete settings archive of the October 2026 PTFE/PVDF PMF recomputation (recompute PC, RTX 3090,
                          GROMACS 2026.3): three zips (build+logs, PTFE windows, PVDF windows) holding every .mdp/.gro/.top/.ndx,
                          cell equilibration, SMD, all umbrella windows (umb.mdp, umb.gro), WHAM in/out, pullx/pullf series,
                          plus ENVIRONMENT (GPU/GROMACS/conda) and the MANIFEST of the authors' full archive (687 files;
                          the published zips omit mdrun logs and .tpr files, see README_curation.md inside)
data/
  raw_results.csv         80 replicates: polymer, compound, seed, Coulomb, LJ, total short-range IE (kJ/mol)
  results_table.csv, experiment_summary.json, tier.json
  pmf/                    WHAM profiles, bootstrap profiles, window histograms, reports, slab checks (see data/pmf/README.md)
  pmf_result.json         four-system PMF summary (supersedes the 2026-08 single-system entry)
  rq1_summary.json, rq3_docking.json, tfeof_meta.json, surface_summary.json
  realization2/           second slab realisation: per-replicate outputs, IE time series, ECTFE/PE slabs, re-derived topologies
  charge_model/           charges_alt.json, Coulomb/LJ under three charge models, rerun time series
```

## Reproducing

```bash
pip install -r requirements.txt
python code/rq1_massbalance.py                                   # mass-balance bound
python code/analyze_screen.py data/raw_results.csv out           # screen statistics
python code/make_fig3.py out out/fig                             # evidence tier (tier.json) + SI Fig S2
python code/make_fig4.py                                         # Fig 3 (mass balance)
python code/make_fig2.py out data/raw_results.csv out/fig data/pmf/pmf_summary.json   # Fig 2
python code/make_figS1.py out/fig PTFE=data/pmf/PTFE_prof.xvg:data/pmf/PTFE_bs.xvg:1.147:3.674 \
    PVDF=data/pmf/PVDF_prof.xvg:data/pmf/PVDF_bs.xvg:1.512:3.642 ECTFE=data/pmf/ECTFE_prof.xvg:data/pmf/ECTFE_bs.xvg:1.653:3.966 \
    PE=data/pmf/PE_prof.xvg:data/pmf/PE_bs.xvg:1.408:3.677                                   # Fig S1
python code/make_numbers_tex.py out data/pmf/pmf_summary.json    # numbers.tex macros
python code/make_si_tables.py out data/raw_results.csv data/pmf/pmf_summary.json
python code/emit_md_inputs.py                                     # .mdp / packmol inputs of the production runs
```

Fixed seeds (`20260710` Monte Carlo and conformer search, `7` bootstrap) make the numeric outputs reproducible.

## Provenance and integrity notes

- All 80 screen trajectories and the ECTFE and PE PMFs were computed in September 2026 on the authors' GPU
  workstation; the PTFE and PVDF PMF systems were rebuilt and recomputed in October 2026 on a second machine with
  the identical two-stage slab protocol after a hardware failure destroyed the first workstation. Rebuilt dry-slab
  interior densities agree with SI Table S1 to −1.0 % (PTFE) and +1.3 % (PVDF) (`data/pmf/*_dryslab.txt`).
- Interaction energies are ordinal proxies, not free energies; the absolute scale comes from the PMFs. The PMF
  protocol, valid ranges and block analysis are in `data/pmf/report_*.txt` and SI section S4.
- PMF provenance (all runs ever made, see SI): 2026-08 PVDF −13.8 kJ/mol (k = 2000, 19 windows, 5 ns, different
  protocol, superseded); 2026-09 first-workstation PTFE/PVDF runs (lost before completion; overlap criterion not met,
  not reported); 2026-10 rebuilt PTFE/PVDF runs = the reported values. `batch_run.py` also carries a PP arm that
  produced no completed trajectories and is not part of the study.

### Not archived

- **Raw trajectories** (all systems) — size; those of the first workstation were lost with it.
- **Force-field topologies and slab/cell files of the September 2026 screen** (ECTFE, PE, PFBS, HFPO-DA, 6:2 FTOH systems
  and the September PTFE/PVDF cells) were lost with the first workstation. `code/build/` archives the deterministic
  structure generation and the acpype protocol that produced them; because AM1-BCC charges depend on the input
  geometry, a re-derivation is equivalent in protocol, not byte-identical. The ECTFE and PE PMF outputs
  (`data/pmf/`) are complete; their window inputs are not archived.
- The October 2026 PTFE, PVDF and PFOA topologies, slabs, cells and all PMF window inputs **are** archived
  (`topologies/`, `slabs/`, `md_inputs/`, `archive/`).

## License

Code: MIT (`LICENSE`). Data and figures: CC BY 4.0.

## Citation

Cite the manuscript above and the archived release (`CITATION.cff`; Zenodo concept DOI 10.5281/zenodo.23838772).

## Contact

Corresponding author: Ku Kang — bisu9082@gmail.com · Repository: HanYeol Nam (@METX-915)
