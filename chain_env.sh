#!/bin/bash

# MIT_HOME is your MITgcm_BFM_chain directory, the root of the git repository
export MIT_HOME=/home/itfc/MODEL/NWMed
# your machine name
export MIT_HOSTNAME=atos
# your work-launch (e.g., "scratch" area) directory in a parallel filesystem
export MIT_WORKDIR=/ec/res4/scratch/itfc/MITgcm_BFM/WORK_benchmark
export MIT_VERSION_NUMBER=1
export MIT_STAGE=devel
alias mitcd="cd $MIT_HOME ; pwd"
