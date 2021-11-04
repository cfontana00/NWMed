#!/bin/bash

# MIT_HOME is your MITgcm_BFM_chain directory, the root of the git repository
export MIT_HOME=/home/innocenti/MITgcm_BFM/NWMed
# your machine name
export MIT_HOSTNAME=hpc-fe1
# your work-launch (e.g., "scratch" area) directory in a parallel filesystem
export MIT_WORKDIR=/home/innocenti/MITgcm_BFM/WORK_NWMed
export MIT_VERSION_NUMBER=1
export MIT_STAGE=devel
alias mitcd="cd $MIT_HOME ; pwd"
