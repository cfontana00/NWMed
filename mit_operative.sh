#!/bin/bash

source /home/itai/.bash_profile
source /home/itai/MITgcm_BFM/NWMed/chain_env.sh

today=$(date -d "+ 1 days" +'%Y%m%d')

/home/itai/MITgcm_BFM/NWMed/bin/mit_start.ksh --pass --job-multiple --force --try-resume --without-kill --rundate $today
