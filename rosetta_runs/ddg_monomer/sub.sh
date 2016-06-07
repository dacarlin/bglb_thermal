#!/bin/bash
#
#SBATCH --job-name=ddg
#SBATCH --output=logs/slurm-%A_%a.out
#SBATCH --error=logs/slurm-%A_%a.err
#SBATCH -p gc128 

cd out # sigh 
MUT=$( sed -n "$SLURM_ARRAY_TASK_ID p" ../list ) 
/share/work/rosetta/source/bin/ddg_monomer.linuxgccrelease @flags -resfile ../res/${MUT}.res
