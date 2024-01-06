#!/bin/bash
#SBATCH --job-name=GMMC_batch
#SBATCH --partition=standard
#SBATCH --account=bakulski0
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0-5:0:0
#SBATCH --output=/home/lutiffan/GMMCResults/gmmc.log

## to run, use sbatch --export=start=a,end=b,ps=0.01 /home/lutiffan/runGMMCbatch.sh

module load Rtidyverse
Rscript --no-save --no-restore "/home/lutiffan/GMMCScripts/runGMMCbatch.R"
