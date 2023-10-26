#!/bin/bash
#SBATCH --account=hd_mh280
#SBATCH --job-name PANACEA_networkcalc
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --mem 40GB
#SBATCH --time 7-00:00:00
#SBATCH --mail-type=FAIL,END

echo "I, ${USER}, am running on host ${HOSTNAME}"
 
echo "STARTING JOB"

source ${HOME}/.bashrc
conda activate mthesis

python3 PANACEA_analysis_cluster.py

echo "FINISHED JOB"