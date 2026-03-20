#!/bin/bash
#SBATCH -J Point
#SBATCH -p first
#SBATCH -N 1 
#SBATCH -n 1                   # 单点 → 只开 1 个任务
#SBATCH --cpus-per-task=1      # 每个任务用 1 核
#SBATCH --mem=64G              # 内存缩小，单点不需要 256G
#SBATCH -o test_o.%j
#SBATCH -e test_e.%j
##SBATCH --exclusive
#SBATCH -t 6-24:00
##SBATCH --nodelist=tms08

ulimit -s unlimited

echo "作业开始: $(date '+%Y-%m-%d %H:%M:%S')"
echo "作业ID: $SLURM_JOB_ID, 节点: $SLURM_NODELIST, 核数: $SLURM_NTASKS"

set -e  
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
./run/colm.x "${NAMELIST_FILE}" > /dev/null
wait  # 确保进程完全结束
echo "实验1完成: $(date '+%Y-%m-%d %H:%M:%S')"

echo "所有实验完成: $(date '+%Y-%m-%d %H:%M:%S')"