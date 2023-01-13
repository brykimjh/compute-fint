#!/bin/bash
#####  Constructed by HPC everywhere #####
#SBATCH --mail-user=brykim@iu.edu
#SBATCH --nodes=!NNODE
#SBATCH --ntasks-per-node=!NTASK
#SBATCH --time=!WALLTIME
##SBATCH --partition=general
#SBATCH --job-name=g16
#SBATCH --output=slurm2.out

######  Module commands #####
module load gaussian/16

######  Job commands go below this line #####
date

WRKDIR=!WRKDIR
DEST=!DEST
SCRATCH=!SCRATCH

rm -rf $SCRATCH/*

cd $SCRATCH

for (( INODE=1; INODE<=!NNODE; INODE++ )); do
    srun -N 1 -n 1 -c !NTASK python3 $WRKDIR/g16/g16.py $INODE !NCALC !NCONFIG $WRKDIR/$DEST &
done
wait

rm -rf $SCRATCH/*

date

