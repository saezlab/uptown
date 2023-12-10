#!/bin/bash
#SBATCH --job-name DECRYPTM_files
#SBATCH --time 2-00:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=16G

source ${HOME}/.bashrc
conda activate mthesis

python decryptm_downloader.py