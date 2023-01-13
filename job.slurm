#!/bin/bash
#####  Constructed by HPC everywhere #####
#SBATCH --mail-user=brykim@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=0-4:0:00
##SBATCH --partition=debug
#SBATCH --job-name=run
#SBATCH --output=slurm1.out

######  Module commands #####

######  Job commands go below this line #####
date
python3 run.py
date

