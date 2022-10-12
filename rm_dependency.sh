#!/bin/bash

#today=$(date +'%Y%m%d')
#jobid=$(awk -F ":" '/post/{print $NF}' /home/innocenti/MITgcm_BFM/NWMed/log/{$today}/mit.hpc-fe1.001.slurm.jobids)

if [ -f "/home/innocenti/MITgcm_BFM/NWMed/wrkdir/POSTPROC/postproc_done.txt" ]; then
    sbatch /home/innocenti/MITgcm_BFM/NWMed/run_rm_dir.sh
fi
