#!/bin/zsh
#SBATCH -N 1
#SBATCH -t 10000
#SBATCH -p slurm
#SBATCH -A sif
#SBATCH -J SIP_update_%a
#SBATCH -o SIP_update_%a.out

source  /etc/profile.d/modules.sh
module purge
module load matlab

cd /pic/projects/sif/dalei/code/SIP
SLURM_ARRAY_TASK_ID=$((${SLURM_ARRAY_TASK_ID} + 0))
echo $SLURM_ARRAY_TASK_ID
matlab -r -nodisplay "Retrieval_LAI\($SLURM_ARRAY_TASK_ID\)"

