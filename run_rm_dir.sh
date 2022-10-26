#!/bin/bash

#SBATCH --job-name=rm
#SBATCH --nodes=4
#SBATCH --ntasks=112
#SBATCH -o /home/innocenti/MITgcm_BFM/NWMed/log/rm_dir/log.%j
#SBATCH -e /home/innocenti/MITgcm_BFM/NWMed/log/rm_dir/err.%j
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a.innocenti@lamma.toscana.it
#SBATCH --time=144:00:00

source /home/itai/.bash_profile
source /home/itai/MITgcm_BFM/NWMed/HOST/atos/venvs/venv1/bin/activate

srun -n 112 python /home/itai/MITgcm_BFM/NWMed/rm_run_directory.py
