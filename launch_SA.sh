#!/bin/bash
#SBATCH -J launcher
#SBATCH -p second
#SBATCH -N 1 -n 1
#SBATCH --mem=1G
#SBATCH -t 0-00:10
#SBATCH -o launcher_%j.log

# ════════════════════════════════════════
#  本次调参设置（只需改这里）
# ════════════════════════════════════════
K_REF_CO2="1.0e-3"
mkdir -p ./namelists/260316CO2/K_REF_CO2/logs

# ════════════════════════════════════════
#  情景定义（OFF → ER → ALL）
# ════════════════════════════════════════
SCENARIOS=("OFF" "ER" "ALL")

declare -A SC_ER
declare -A SC_NITR
declare -A SC_VC

SC_ER["OFF"]=".false.";  SC_NITR["OFF"]=".false.";  SC_VC["OFF"]=".false."
SC_ER["ER"]=".true.";    SC_NITR["ER"]=".false.";   SC_VC["ER"]=".false."
SC_ER["ALL"]=".true.";   SC_NITR["ALL"]=".true.";   SC_VC["ALL"]=".true."

# ════════════════════════════════════════
#  N 梯度（字面量）
# ════════════════════════════════════════
FERT_VALS=(0 1 2)
N_LABELS=("N0" "N1" "N2")

# ════════════════════════════════════════
#  提交（情景外循环，N梯度内循环）
# ════════════════════════════════════════

TOTAL=0
submitted_jobs=0  # 提交的作业计数
for sc in "${SCENARIOS[@]}"; do
    for i in "${!FERT_VALS[@]}"; do
        FERT=${FERT_VALS[$i]}
        NLABEL=${N_LABELS[$i]}
        CASE_NAME="US-Ne1_4_1land9_spin200_${sc}_${NLABEL}"

        PARAMS="ALL"
        PARAMS="${PARAMS},CASE_NAME=${CASE_NAME}"
        PARAMS="${PARAMS},FERT=${FERT}"
        PARAMS="${PARAMS},USE_SOILACID=.true."
        PARAMS="${PARAMS},SOILACID_ER=${SC_ER[$sc]}"
        PARAMS="${PARAMS},SOILACID_NITR=${SC_NITR[$sc]}"
        PARAMS="${PARAMS},SOILACID_VC=${SC_VC[$sc]}"
        PARAMS="${PARAMS},SOILACID_K_REF_CO2=${K_REF_CO2}"


        JOB_ID=$(sbatch --job-name="${CASE_NAME}" --export="${PARAMS}" job_SA_single.slurm | awk '{print $NF}')
        echo "已提交 [${JOB_ID}]  ${CASE_NAME}  FERT=${FERT}  ER=${SC_ER[$sc]}  NITR=${SC_NITR[$sc]}  VC=${SC_VC[$sc]}"
        submitted_jobs=$((submitted_jobs + 1))
        TOTAL=$((TOTAL + 1))

        # 检查是否已提交超过3个作业，如果是则等待当前作业完成
        if [ $submitted_jobs -ge 3 ]; then
            echo "已提交 3 个作业，等待它们完成..."
            while [ $(squeue -u $USER | grep -c "PENDING") -gt 0 ]; do
                sleep 5  # 等待5秒
            done
            submitted_jobs=0  # 重置计数器，准备提交下一批
        fi
    done
done

echo "========================================"
echo "共提交 ${TOTAL} 个作业  K_REF_CO2=${K_REF_CO2}"
echo "========================================"
