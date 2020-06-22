#!/bin/zsh
#SBATCH -N 1
#SBATCH -t 1000
#SBATCH -p slurm
#SBATCH -A sif
#SBATCH -J MAIAC_data_extract
#SBATCH -o  MAIAC_data_extract.out

source  /etc/profile.d/modules.sh
module purge
module load matlab

cd /pic/projects/sif/dalei/code/SIP/

matlab -r -nodesktop -nosplash  "Extract_EPIC_MAIAC_data"