#!/bin/bash

#SBATCH --job-name=rm
#SBATCH --nodes=4
#SBATCH --ntasks=112
#SBATCH -o /home/innocenti/MITgcm_BFM/NWMed/log/rm_dir/log.%j
#SBATCH -e /home/innocenti/MITgcm_BFM/NWMed/log/rm_dir/err.%j
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=a.innocenti@lamma.toscana.it
#SBATCH --time=144:00:00

source /home/innocenti/.bash_profile
source /home/innocenti/MITgcm_BFM/NWMed/HOST/hpc-fe1/venvs/venv1/bin/activate

mpirun -np 112 python /home/innocenti/MITgcm_BFM/NWMed/rm_run_directory.py
