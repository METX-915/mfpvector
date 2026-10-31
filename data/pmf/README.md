# Umbrella-sampling PMF data (PFOA on PTFE, PVDF, ECTFE, PE)

One protocol for all four systems (SI section S4): 8 ns lateral-fixed cell equilibration, inward SMD 0.002 nm/ps,
0.13 nm windows, k = 1000 kJ mol-1 nm-2, 8 ns per window (ECTFE 13 ns), first 4 ns discarded, gmx wham with
200 b-hist bootstrap resamples. dG_ads = minimum relative to the mean of the outer 25 % of the valid profile.

| file | content |
|---|---|
| `<POLY>_prof.xvg` | `gmx wham -o` profile (z in nm, PMF in kJ/mol) over the analysed segment |
| `<POLY>_bs.xvg` | `gmx wham -bsres` bootstrap mean and standard deviation |
| `<POLY>_hist.xvg` | `gmx wham -hist` window histograms (overlap check) |
| `report_<POLY>_*blocks.txt` | report_pmf.py output: dG_ads, plateau, valid range, min. overlap, block (first/second half) analysis |
| `PTFE/PVDF_cell.txt`, `_dryslab.txt` | rebuilt-slab checks: interior density, z_surf, interface width |
| `pmf_summary.json` | final values used in Fig. 2d and Table S7 |

Final values: PTFE -30.72 +/- 3.52, PVDF -26.62 +/- 3.18, ECTFE -32.17 +/- 2.32, PE -31.35 +/- 3.92 kJ/mol.
PTFE and PVDF systems were rebuilt (identical two-stage slab protocol) after a hardware failure destroyed the original
trajectories; ECTFE and PE are the original runs. Scripts: `code/pmf/`.
