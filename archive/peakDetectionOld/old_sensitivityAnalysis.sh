#!/bin/bash
#SBATCH --partition=standard
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --job-name=sensitivity_analysis
#SBATCH --mem=8G
#SBATCH --time=0-1:0:0
#SBATCH --output=/home/lutiffan/sensitivityAnalysis/sa.log

## to run, use sbatch --export=prop=y,space=z,start=a,end=b /home/lutiffan/sensitivityAnalysis.sh
# R CMD BATCH --no-save --no-restore "/home/lutiffan/sensitivityAnalysis.R"
module load Rtidyverse
Rscript --no-save --no-restore "/home/lutiffan/sensitivityAnalysis.R"