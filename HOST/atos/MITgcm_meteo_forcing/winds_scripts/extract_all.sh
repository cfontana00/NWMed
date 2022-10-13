#!/bin/bash

# Usage:
# ./extract_wind10m.sh gfilter ifile ofile
# ./winds_scripts/extract_all.sh 2021070900_arw_ecm_3km.grb temp.grb

ifile=$1
ofile=$2

grib_filter -o temp.grb winds_scripts/filters/grib.filter.wrf03ecm.test $ifile
grib_to_netcdf -D NC_FLOAT temp.grb -o output/20220113/$ofile
