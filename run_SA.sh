#!/bin/bash

echo "========== 开始实验1 =========="
export CASE_NAME="US-Ne1_4_1land9_spin200_OFF_N6"
export FERT="1.5"
export USE_SOILACID=".true."
export SOILACID_VC=".false."
export SOILACID_ER=".false."
export SOILACID_NITR=".false."
source ./update_namelist.sh
python3 model_SA_0106.py
./run/mksrfdata.x "${NAMELIST_FILE}"
./run/colm.x "${NAMELIST_FILE}" > out_OFF_N6.log
wait  # 确保进程完全结束
echo "实验1完成: $(date '+%Y-%m-%d %H:%M:%S')"

# echo "========== 开始实验2 =========="
# export CASE_NAME="US-Ne1_4_1land9_nospin_ER_C3"
# export FERT="1.5"
# export USE_SOILACID=".true."
# export SOILACID_VC=".false."
# export SOILACID_ER=".true."
# export SOILACID_NITR=".false."
# source ./update_namelist.sh
# python3 model_SA_0106.py
# ./run/mksrfdata.x "${NAMELIST_FILE}"
# ./run/colm.x "${NAMELIST_FILE}" > /dev/null
# wait  # 确保进程完全结束
# echo "实验2完成: $(date '+%Y-%m-%d %H:%M:%S')"

# echo "========== 开始实验3 =========="
# export CASE_NAME="US-Ne1_4_1land9_nospin_NITR_C3"
# export FERT="1.5"
# export USE_SOILACID=".true."
# export SOILACID_VC=".false."
# export SOILACID_ER=".false."
# export SOILACID_NITR=".true."
# source ./update_namelist.sh
# python3 model_SA_0106.py
# ./run/mksrfdata.x "${NAMELIST_FILE}"
# ./run/colm.x "${NAMELIST_FILE}" > /dev/null
# wait  # 确保进程完全结束
# echo "实验3完成: $(date '+%Y-%m-%d %H:%M:%S')"