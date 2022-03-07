#!/bin/bash

source /home/innocenti/.bash_profile
source /home/innocenti/MITgcm_BFM/NWMed/chain_env.sh

today=$(date +'%Y%m%d')

/home/innocenti/MITgcm_BFM/NWMed/bin/mit_start.ksh --pass --job-multiple --force --try-resume --without-kill --rundate $today
