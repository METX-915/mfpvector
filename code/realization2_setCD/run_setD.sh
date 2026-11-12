#!/usr/bin/env bash
# Set D: GenX on all four second-realisation slabs, 4 cells x 5 seeds. Completes the 16-cell second realisation.
cd "$(dirname "$0")"
export MFP_ROOT="$PWD" MFP_POLYS="PTFE,PVDF,ECTFE,PE" MFP_PFAS="GenX" MFP_CSV="raw_results_realization2.csv"
exec python code/batch_run.py
