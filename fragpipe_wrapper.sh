#!/bin/bash
#SBATCH --job-name fragpipe_fullprot
#SBATCH --time 2-00:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=20G

source ${HOME}/.bashrc
conda activate mthesis

/net/data.isilon/ag-saez/bq_vpaton/mthesis/fragpipe/bin/fragpipe --headless --workflow /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe_fullprot.workflow --manifest /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe-files.fp-manifest --workdir /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/all_phospho_processed_full_results
/net/data.isilon/ag-saez/bq_vpaton/mthesis/fragpipe/bin/fragpipe --headless --workflow /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe_tfs.workflow --manifest /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe-files.fp-manifest --workdir /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/all_phospho_processed_tfs_results