#!/bin/bash
#SBATCH --job-name networkcalc
#SBATCH --time 7-00:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --queue=long

echo "I, ${USER}, am running on host ${HOSTNAME}"
 
echo "STARTING JOB"

source ${HOME}/.bashrc
conda activate mthesis
cd /net/data.isilon/ag-saez/bq_vpaton/mthesis

# PANACEA
# nextflow -C uptown_network.config run uptown_network_calc.nf -resume -profile cluster --network network_collectri.sif --source_file panacea_sources.tsv --tf_file tf_activity_results.tsv --iterations 1000
nextflow -C uptown_network.config run uptown_network_eval.nf -resume -profile cluster --network network_collectri.sif --dirpath ./results/ --tf_file tf_activity_results.tsv --offtargets panacea_offtargets.tsv --iterations 1000

# decryptm
# nextflow -C uptown_network.config run uptown_network_calc.nf -resume -profile cluster --network network_collectri.sif --source_file decryptm_sources.tsv --tf_file proteomics_targets.tsv --iterations 1000 --dataset decryptm
nextflow -C uptown_network.config run uptown_network_eval.nf -resume -profile cluster --network network_collectri.sif --dirpath results_decryptm/ --phospho phospho_prots.tsv --offtargets decryptm_offtargets.tsv --dataset decryptm --tf_file proteomics_targets.tsv --iterations 1000 --ec50_file ec50_info.tsv --sources_file panacea_sources.tsv


echo "FINISHED JOB"