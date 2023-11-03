#!/bin/bash
#SBATCH --job-name PANACEA_networkcalc
#SBATCH --time 7-00:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --queue=long

echo "I, ${USER}, am running on host ${HOSTNAME}"
 
echo "STARTING JOB"

source ${HOME}/.bashrc
conda activate mthesis
cd /net/data.isilon/ag-saez/bq_vpaton/mthesis

# Executing an R script
# nextflow -C PANACEA_network.config run PANACEA_network_calc.nf -resume -profile cluster --network network_collectri.sif --source_file panacea_sources.tsv --tf_file tf_activity_results.tsv --iterations 1000
nextflow -C PANACEA_network.config run PANACEA_network_eval.nf -resume -profile cluster --network network_collectri.sif --dirpath ./results/ --tf_file tf_activity_results.tsv --offtargets panacea_offtargets.tsv --iterations 1000

echo "FINISHED JOB"