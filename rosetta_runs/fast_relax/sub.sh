#!/bin/bash
#
#SBATCH --job-name=relax
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err

MUT=$( sed -n "$SLURM_ARRAY_TASK_ID p" list ) 
module load rosetta 
rosetta_scripts.linuxgccrelease @flags $MUT
