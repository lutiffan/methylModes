#!/bin/bash
#SBATCH --job-name=sensitivity_analysis
#SBATCH --partition=standard
#SBATCH --account=bakulski0
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0-2:0:0
#SBATCH --output=/home/lutiffan/peakDetectionResults/pd.log

## to run, use sbatch --export=data=(child/teen),rangeStart=a,rangeEnd=b,propStart=c,propEnd=d,spaceStart=e,spaceEnd=f,spacing=z /home/lutiffan/sensitivityTests.sh

module load Rtidyverse
Rscript --no-save --no-restore "/home/lutiffan/peakDetectionScripts/sensitivityTests.R"
