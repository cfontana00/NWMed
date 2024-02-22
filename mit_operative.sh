#!/bin/bash

source /home/itfc/.bash_profile
source /home/itfc/MITgcm_BFM/NWMed/chain_env.sh

today=$(date -d "+ 1 days" +'%Y%m%d')

/home/itfc/MITgcm_BFM/NWMed/bin/mit_start.ksh --pass --job-multiple --force --try-resume --without-kill --rundate $today
