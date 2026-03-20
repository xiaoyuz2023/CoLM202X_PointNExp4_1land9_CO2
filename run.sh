#!/bin/bash
echo "========== 开始实验0: PH打开 =========="
export CASE_NAME="US-Ne1_4_1land9_nospin_Ab12_test"
export FERT="1.1"
export USE_SOILACID=".true."
export SOILACID_VC=".false."
export SOILACID_ER=".false."
export SOILACID_NITR=".true."
source ./update_namelist.sh
python3 model_SA_0106.py
./run/mksrfdata.x "${NAMELIST_FILE}"
./run/colm.x "${NAMELIST_FILE}" > Ab12_test_out.log
wait  # 确保进程完全结束
echo "实验0完成: $(date '+%Y-%m-%d %H:%M:%S')"

