#!/usr/bin/env bash
# Set A: PFBS and 6:2 FTOH on the rebuilt PTFE and PVDF slabs (second realisation), 5 seeds x 6 ns each = 20 runs.
# Resumable. Output: results/raw_results_realization2.csv, results/progress.log, systems/<P>_<X>/s<k>_ie.xvg
cd "$(dirname "$0")"
export MFP_ROOT="$PWD" MFP_POLYS="PTFE,PVDF" MFP_PFAS="PFBS,FTOH62" MFP_CSV="raw_results_realization2.csv"
exec python code/batch_run.py
