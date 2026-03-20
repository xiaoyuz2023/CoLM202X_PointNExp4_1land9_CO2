#!/bin/bash
set -e
ulimit -s unlimited

# ════════════════════════════════════════
#  只需改这里
# ════════════════════════════════════════
export CASE_NAME="US-Ne1_4_1land9_spin200_OFF_N0"

export FERT="0"
export USE_SOILACID=".true."
export SOILACID_ER=".false."
export SOILACID_NITR=".false."
export SOILACID_VC=".false."

export SOILACID_EXCHANGE_CALIB="0.0"


echo "========================================"
echo "  开始: $(date '+%Y-%m-%d %H:%M:%S')"
echo "  CASE_NAME     : ${CASE_NAME}"
echo "  FERT          : ${FERT}"
echo "  ER / NITR / VC: ${SOILACID_ER} / ${SOILACID_NITR} / ${SOILACID_VC}"
echo "  pH_opt / sens : ${SOILACID_PH_OPT} / ${SOILACID_PH_SENS}"
echo "========================================"

echo "---------- [1/3] update_namelist ----------"
source ./update_namelist.sh

echo "---------- [2/3] model_SA_0106.py ----------"
python3 model_SA_0106.py

echo "---------- [3/3] mksrfdata ----------"
./run/mksrfdata.x "${NAMELIST_FILE}"

echo "---------- [4/4] colm ----------"
./run/colm.x "${NAMELIST_FILE}" > out_OFF_N0.log

echo "========================================"
echo "  完成: $(date '+%Y-%m-%d %H:%M:%S')"

echo "========================================"

