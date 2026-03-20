#!/bin/bash
#SBATCH -J CoLM_SA
#SBATCH -p second
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH -o logs/%x_%j.o
#SBATCH -e logs/%x_%j.e
#SBATCH -t 6-24:00

ulimit -s unlimited

echo "作业开始: $(date '+%Y-%m-%d %H:%M:%S')"
echo "作业ID: $SLURM_JOB_ID  节点: $SLURM_NODELIST"

# ── 参数（launcher 通过 --export 传入优先；未传则用此处默认值）──
: "${CASE_NAME:=Ne1_4_1land9_spin200_OFF_N0}"
: "${FERT:=0}"
: "${USE_SOILACID:=.true.}"
: "${SOILACID_ER:=.false.}"
: "${SOILACID_NITR:=.false.}"
: "${SOILACID_VC:=.false.}"
: "${SOILACID_K_REF_CO2:=1.0e-3}"

echo "========================================"
echo "  实验名: ${CASE_NAME}"
echo "  FERT  : ${FERT}"
echo "  情景  : ER=${SOILACID_ER}  NITR=${SOILACID_NITR}  VC=${SOILACID_VC}"
echo "  调参  : K_REF_CO2=${SOILACID_K_REF_CO2}"
echo "========================================"

set -e

source ./update_namelist.sh

python3 model_SA_0106.py

./run/mksrfdata.x "${NAMELIST_FILE}"

./run/colm.x "${NAMELIST_FILE}" > /dev/null

echo "作业完成: $(date '+%Y-%m-%d %H:%M:%S')"
