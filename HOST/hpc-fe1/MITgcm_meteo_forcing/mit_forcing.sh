#!/bin/bash

#SBATCH --job-name=mit_prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -o /home/innocenti/MITgcm_BFM/NWMed/HOST/hpc-fe1/MITgcm_meteo_forcing/log.%j
#SBATCH -e /home/innocenti/MITgcm_BFM/NWMed/HOST/hpc-fe1/MITgcm_meteo_forcing/err.%j
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=a.innocenti@lamma.toscana.it
#SBATCH --time=144:00:00
#SBATCH -Q

# Load common profile
#. @@(I:MIT_HOME)/bin/mit_profile.inc

# Rundate definition
#mit_set_run

# start
#mit_start
#mit_prex "source $MIT_VENV_1/bin/activate"
#export PYTHONPATH=$PYTHONPATH:$MIT_BITSEA

#DATESTART=$( date -d "${MIT_RUNDATE}  -  7 days  " +%Y%m%d-%H:%M:%S )
#DATE__END=$( date -d "${MIT_RUNDATE}  + 72 hours " +%Y%m%d-%H:%M:%S )

MIT_RUNDATE=20220110
MIT_HOSTDIR=/home/innocenti/MITgcm_BFM/NWMed/HOST/hpc-fe1

FRC_DIR=$MIT_HOSTDIR/MITgcm_meteo_forcing
FRC_DIR_GRIB=$FRC_DIR/meteo_grib
FRC_DIR_NC=$FRC_DIR/meteo_nc/$MIT_RUNDATE
FRC_DIR_TXT=$FRC_DIR/meteo_txt/$MIT_RUNDATE
#[ -d $FRC_DIR_GRIB ] || mkdir -p $FRC_DIR_GRIB;
#[ -d $FRC_DIR_NC ] || mkdir -p $FRC_DIR_NC;
#[ -d $FRC_DIR_TXT ] || mkdir -p $FRC_DIR_TXT;
wrf_grb=00_arw_ecm_3km.grb
wrf_nc=00_arw_ecm_3km.nc
WRF_DIR=/OCEANASTORE/progetti/meteop/data/wrf03ecm
#for i in {0..7}
#do
#    DATE_I=$( date -d "${MIT_RUNDATE}  -  ${i} days  " +%Y%m%d )
#    ifile=$DATE_I$wrf_grb
#    ofile=$DATE_I$wrf_nc
#    cp $WRF_DIR/$ifile $FRC_DIR_GRIB 
#    grib_filter -o $FRC_DIR/temp.grb $FRC_DIR/winds_scripts/filters/grib.filter.wrf03ecm.test $FRC_DIR_GRIB/$ifile
#    grib_to_netcdf -D NC_FLOAT $FRC_DIR/temp.grb -o $FRC_DIR_NC/$ofile

#done

#rm $FRC_DIR_GRIB/*

# 
# Sets up the MCR environment
#         
MCRROOT=/cluster/universal/matlab/opt/matlab_R2012b/MATLAB_Compiler_Runtime/v80/
MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/runtime/glnxa64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}
XAPPLRESDIR=${MCRROOT}/X11/app-defaults
export LD_LIBRARY_PATH
export XAPPLRESDIR

$FRC_DIR/make_meteo_mitgcm_fc $MIT_RUNDATE

#mit_exit "$errors"


