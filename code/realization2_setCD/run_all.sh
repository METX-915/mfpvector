#!/usr/bin/env bash
# 1) charge-model sensitivity on the existing trajectories (~2 h, CPU)
# 2) set C (~14-27 h, GPU)  3) set D (~14-27 h, GPU)  4) charge sensitivity on the new trajectories (~3 h, CPU)
cd "$(dirname "$0")"
echo "[$(date '+%m-%d %H:%M')] STEP 1 charge sensitivity" ; python charge_sensitivity.py > cs1.out 2>&1
echo "[$(date '+%m-%d %H:%M')] STEP 2 set C" ; bash run_setC.sh > setC.out 2>&1
echo "[$(date '+%m-%d %H:%M')] STEP 3 set D" ; bash run_setD.sh > setD.out 2>&1
echo "[$(date '+%m-%d %H:%M')] STEP 4 charge sensitivity (new cells)" ; python charge_sensitivity.py > cs2.out 2>&1
echo "[$(date '+%m-%d %H:%M')] ALL STEPS DONE"
