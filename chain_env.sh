#!/bin/bash

# MIT_HOME is your MITgcm_BFM_chain directory, the root of the git repository
export MIT_HOME=/home/itai/MITgcm_BFM/NWMed
# your machine name
export MIT_HOSTNAME=atos
# your work-launch (e.g., "scratch" area) directory in a parallel filesystem
export MIT_WORKDIR=/ec/res4/scratch/itai/MITgcm_BFM/WORK_NWMed
export MIT_VERSION_NUMBER=1
export MIT_STAGE=devel
alias mitcd="cd $MIT_HOME ; pwd"

export PATH=$HOME/jq-1.6:$PATH
