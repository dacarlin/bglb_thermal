#!/bin/bash
#
#SBATCH --job-name=ddg
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err

MUT=$( sed -n "$SLURM_ARRAY_TASK_ID p" list ) 
/share/work/rosetta/source/bin/ddg_monomer.linuxgccrelease @flags -resfile resfiles/$MUT
