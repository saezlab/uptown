#!/bin/bash
#SBATCH --job-name=fragpipe
#SBATCH --time=3-24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mem=150GB

echo "I, ${USER}, am running on host ${HOSTNAME}"
 
echo "STARTING JOB"
 
# Load conda module
source ${HOME}/.bashrc
conda activate mthesis

# /net/data.isilon/ag-saez/bq_vpaton/mthesis/fragpipe/bin/fragpipe --headless --workflow /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe.workflow --manifest /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Afatinib_30min_R1.fp-manifest --workdir /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Afatinib_30min_R1 --threads 2
# /net/data.isilon/ag-saez/bq_vpaton/mthesis/fragpipe/bin/fragpipe --headless --workflow /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe.workflow --manifest /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Afatinib_30min_R2.fp-manifest --workdir /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Afatinib_30min_R2 --threads 2
# /net/data.isilon/ag-saez/bq_vpaton/mthesis/fragpipe/bin/fragpipe --headless --workflow /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe.workflow --manifest /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Afatinib_30min_R3.fp-manifest --workdir /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Afatinib_30min_R3 --threads 2
# /net/data.isilon/ag-saez/bq_vpaton/mthesis/fragpipe/bin/fragpipe --headless --workflow /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe.workflow --manifest /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Afatinib_30min_R4.fp-manifest --workdir /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Afatinib_30min_R4 --threads 2
# /net/data.isilon/ag-saez/bq_vpaton/mthesis/fragpipe/bin/fragpipe --headless --workflow /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe.workflow --manifest /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Dasatinib_30min_R1.fp-manifest --workdir /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Dasatinib_30min_R1 --threads 2
# /net/data.isilon/ag-saez/bq_vpaton/mthesis/fragpipe/bin/fragpipe --headless --workflow /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe.workflow --manifest /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Dasatinib_30min_R2.fp-manifest --workdir /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Dasatinib_30min_R2 --threads 2
# /net/data.isilon/ag-saez/bq_vpaton/mthesis/fragpipe/bin/fragpipe --headless --workflow /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe.workflow --manifest /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Dasatinib_30min_R3.fp-manifest --workdir /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Dasatinib_30min_R3 --threads 2
# /net/data.isilon/ag-saez/bq_vpaton/mthesis/fragpipe/bin/fragpipe --headless --workflow /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe.workflow --manifest /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Dasatinib_30min_R4.fp-manifest --workdir /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Dasatinib_30min_R4 --threads 2
# /net/data.isilon/ag-saez/bq_vpaton/mthesis/fragpipe/bin/fragpipe --headless --workflow /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe.workflow --manifest /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Gefitinib_30min_R1.fp-manifest --workdir /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Gefitinib_30min_R1 --threads 2
# /net/data.isilon/ag-saez/bq_vpaton/mthesis/fragpipe/bin/fragpipe --headless --workflow /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe.workflow --manifest /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Gefitinib_30min_R2.fp-manifest --workdir /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Gefitinib_30min_R2 --threads 2
# /net/data.isilon/ag-saez/bq_vpaton/mthesis/fragpipe/bin/fragpipe --headless --workflow /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe.workflow --manifest /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Gefitinib_30min_R3.fp-manifest --workdir /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Gefitinib_30min_R3 --threads 2
/net/data.isilon/ag-saez/bq_vpaton/mthesis/fragpipe/bin/fragpipe --headless --workflow /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/fragpipe.workflow --manifest /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Gefitinib_30min_R4.fp-manifest --workdir /net/data.isilon/ag-saez/bq_vpaton/mthesis/decryptm/ddPTM_A431_Gefitinib_30min_R4 --threads 2
