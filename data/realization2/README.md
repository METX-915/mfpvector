# Second slab realisation of the screen (SI Section S11, Table S17)

All sixteen cells of the interaction-energy screen repeated on a second, independently built slab per polymer, 5 replicates x 6 ns each,
with the protocol of code/batch_run.py (identical .mdp settings, gen-seed = 1000*seed + 7, CPU -rerun with energygrps POL PFA,
energies averaged over the final 4 ns). GROMACS 2026.3 (conda-forge CUDA build) on an RTX 2070 laptop GPU, 2026-11-04 to 11-12.

- results_realization2.csv : PFOA on the rebuilt PTFE and PVDF slabs (driver code/run_realization2.py, from md_inputs/<SYS>/npt.gro)
- setAB/raw_results_realization2.csv : PFBS and 6:2 FTOH on the rebuilt PTFE and PVDF slabs (set A) and PFOA on new ECTFE and PE slabs (set B),
  driver code/realization2_setAB/batch_run_env.py (= code/batch_run.py with polymer/compound/CSV selected by environment variables and
  packmol's non-zero exit tolerated when its output exists, as in the October rebuild); run_setA.sh / run_setB.sh
- setAB/<P>_<X>_s<k>_ie.xvg : 10 ps time series of Coul-SR and LJ-SR over the analysed 4 ns (401 frames)
- slabs/ECTFE_slab.gro, slabs/PE_slab.gro, *_density.txt : second-realisation dry slabs (two-stage protocol; rho 1.441 and 0.790 g/cm3)
- topologies/{PFBS,FTOH62,ECTFE,PE} : GAFF2/AM1-BCC topologies re-derived with acpype 2023.10 from the archived structures (code/build/)
- raw_results_realization2.csv : 70 replicates of sets A-D (supersedes setAB/raw_results_realization2.csv; set C = PFBS and 6:2 FTOH on ECTFE and PE, set D = GenX on all four polymers; drivers code/realization2_setCD/run_setC.sh, run_setD.sh, run_all.sh)
- setCD_xvg/<P>_<X>_s<k>_ie.xvg, logs/ : time series and run logs of sets C and D (PE-GenX was restarted from its build after a power loss; all five replicates are from the restarted build)
- screen_both_realisations_160.csv : original screen (realisation 1, 80 rows) and second realisation (80 rows) in one table
- topologies/GenX (repository topologies/GenX_GMX.*) : re-derived with acpype 2023.10.27
- Charge-model test of these trajectories: data/charge_model/ (SI Section S12). The PTFE-PFOA and PVDF-PFOA trajectories were not retained, so that test covers 14 cells.
