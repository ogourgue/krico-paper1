#!/bin/bash
#SBATCH --job-name=krico_fs1
#SBATCH --qos=nf
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --output=fs1_%j.out

module load python3
python3 aggregate.py

