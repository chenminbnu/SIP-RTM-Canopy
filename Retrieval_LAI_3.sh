#!/bin/zsh
#SBATCH -N 1
#SBATCH -t 10000
#SBATCH -p slurm
#SBATCH -A sif
#SBATCH -J LAI_3_2
#SBATCH -o LAI_3_2.out

source  /etc/profile.d/modules.sh
module purge
module load matlab

cd /pic/projects/sif/dalei/code/SIP
SLURM_ARRAY_TASK_ID=$((${SLURM_ARRAY_TASK_ID} + 0))
echo $SLURM_ARRAY_TASK_ID
matlab -r -nodisplay "Retrieval_LAI_3"
