#!/bin/bash
#SBATCH -J PointSA
#SBATCH -p second
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

echo "作业ID: $SLURM_JOB_ID"
echo "作业名称: $SLURM_JOB_NAME"
echo "分区: $SLURM_JOB_PARTITION"
echo "节点数: $SLURM_JOB_NUM_NODES"
echo "每节点任务数: $SLURM_NTASKS_PER_NODE"
echo "每任务CPU数: $SLURM_CPUS_PER_TASK"
echo "内存: $SLURM_MEM_PER_NODE"
echo "时间限制: 6-24:00"
echo "节点列表: $SLURM_JOB_NODELIST"
echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"

set -e   # 遇到错误立即退出
# mpirun -n 48 /stu02/zhangsy24/CoLM202X_RegionNExp4_1landspin/run/mksrfdata.x Northeast_soil_acid_P_spinup.nml
# mpirun -n 48 /stu02/zhangsy24/CoLM202X_RegionNExp4_1landspin/run/mkinidata.x Northeast_soil_acid_P_spinup.nml
# mpirun -n 48 /stu02/zhangsy24/CoLM202X_RegionNExp4_1landspin/run/colm.x Northeast_soil_acid_P_spinup.nml

time python3 model_1006.py
/stu02/zhangsy24/CoLM202X_PointNExp4_1land8/run/colm.x US_Ton_soil_acid_P.nml
