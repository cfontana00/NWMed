#!/bin/bash

source /home/itai/MITgcm_BFM/NWMed/chain_env.sh

DATE_FILE="/home/itai/MITgcm_BFM/NWMed/date.txt"

STOP=12

if [ -f $DATE_FILE ]; then
    while IFS=: read line;
    do
	DAY=$( date --date="$line" +%e )
	if [ "$DAY" -lt "$STOP" ]; then
	    /home/itai/MITgcm_BFM/NWMed/bin/mit_start.ksh --pass --job-multiple --force --try-resume --without-kill --rundate $line
	fi

    done < $DATE_FILE
else
    echo "file date.txt not present" > /home/itai/MITgcm_BFM/NWMed/out.txt
fi
