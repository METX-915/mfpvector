#!/usr/bin/env bash
# Set B: build second-realisation ECTFE and PE slabs (two-stage packmol -> compress/anneal, ~7 ns each) and run PFOA, 5 seeds x 6 ns.
# Requires ECTFE.acpype/ and PE.acpype/ (provided) and packmol in the conda env. Run AFTER run_setA.sh finishes.
cd "$(dirname "$0")"
export MFP_ROOT="$PWD" MFP_POLYS="ECTFE,PE" MFP_PFAS="PFOA" MFP_CSV="raw_results_realization2.csv"
exec python code/batch_run.py
