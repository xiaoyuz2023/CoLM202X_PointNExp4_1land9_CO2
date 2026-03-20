#!/bin/bash
#SBATCH -J region
#SBATCH -p first
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=256G
#SBATCH -o test_o.%j
#SBATCH -e test_e.%j
#SBATCH --exclusive
#SBATCH -t 6-24:00
##SBATCH --nodelist=tms08

ulimit -s unlimited



srun /stu02/zhangsy24/CoLM202X_PointNExp4_1land8/run/mksrfdata.x US_Ne1_soil_acid_P.nml
srun /stu02/zhangsy24/CoLM202X_PointNExp4_1land8/run/mkinidata.x US_Ne1_soil_acid_P.nml
srun /stu02/zhangsy24/CoLM202X_PointNExp4_1land8/run/colm.x US_Ne1_soil_acid_P.nml
