#!/usr/bin/env bash
# Set C: PFBS and 6:2 FTOH on the second-realisation ECTFE and PE slabs (slabs reused from set B), 4 cells x 5 seeds.
cd "$(dirname "$0")"
export MFP_ROOT="$PWD" MFP_POLYS="ECTFE,PE" MFP_PFAS="PFBS,FTOH62" MFP_CSV="raw_results_realization2.csv"
exec python code/batch_run.py
