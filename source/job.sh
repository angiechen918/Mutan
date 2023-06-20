#!/usr/bin/env bash
# File       : thermo_stabiliy.sh
# Description: Submit a test job to check your implementation
# Copyright 2023 Harvard University. All Rights Reserved.
#SBATCH --job-name=anqi
#SBATCH --output=anqi_check%j.out
#SBATCH --mem-per-cpu=2048
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=10:00:00
#SBATCH --mail-user=anqichen@g.harvard.edu
#SBATCH --mail-type=ALL



module load python/3.9.12-fasrc01 gcc/12.1.0-fasrc01 openmpi/4.1.3-fasrc01

root=$(pwd -P)
cwd=library_check_${SLURM_JOBID}
mkdir -p ${cwd}
cd ${cwd}
cp -t . ${root}/*.csv ${root}/*.py ${root}/*.fastq.gz ${root}/*.fasta

python -m pip install biopython

python 230523_thermo_library.py

