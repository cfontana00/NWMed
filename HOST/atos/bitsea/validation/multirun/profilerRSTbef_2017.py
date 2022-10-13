#!/usr/bin/env python
# Author: Giorgio Bolzon <gbolzon@ogs.trieste.it>
# Script to generate profiles of model files in
# the same time and locations where instruments
# such as bioFloats, mooring or vessels have been found.

# When imported, this scripts only defines settings for matchup generation.
from instruments.superfloat import FloatSelector
from instruments.var_conversions import LOVFLOATVARS
from instruments.matchup_manager import Matchup_Manager
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from basins.region import Rectangle
# location of input big ave files, usually the TMP directory.
# ave files are supposed to have N3n, O2o and chl

RUN='MULTIVARIATE_24/TEST_04'

INPUTDIR='/gpfs/scratch/userexternal/ateruzzi/' + RUN + \
    '/wrkdir/POSTPROC/output/DA__FREQ_1/TMP/'

# output directory, where aveScan.py will be run.


BASEDIR='/gpfs/scratch/userexternal/ateruzzi/ELAB_DAFloat_24/VALID_float/' \
    + RUN + '/PROFILATORE_RSTbef/'


DATESTART = '20170101'
DATE__END = '20190101'

varmodel = 'P_l'

T_INT = TimeInterval(DATESTART,DATE__END, '%Y%m%d')
TL = TimeList.fromfilenames(T_INT, INPUTDIR,"RSTbefore*13:00*.nc", \
    prefix='RSTbefore.',filtervar=varmodel)

# ALL_PROFILES = FloatSelector(LOVFLOATVARS[varmodel],T_INT, Rectangle(-6,36,30,46))
ALL_PROFILES = FloatSelector(None,T_INT, Rectangle(-6,36,30,46))


vardescriptorfile="/gpfs/scratch/userexternal/ateruzzi/" + \
    "ELAB_DAFloat_24/VALID_float/bit.sea/validation/multirun/" + \
    "VarDescriptorRSTaft_2017.xml"

#This previous part will be imported in matchups setup.

# The following part, the profiler, is executed once and for all.
# It might take some time, depending on length of simulation or size of files.
if __name__ == '__main__':
    # Here instruments time and positions are read as well as model times
    M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)


    profilerscript = BASEDIR + 'jobProfiler.sh'
    aggregatedir="/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v19_3/wrkdir/POSTPROC/output/AVE_FREQ_1/TMP/"
    M.writefiles_for_profiling(vardescriptorfile, profilerscript, aggregatedir=aggregatedir) # preparation of data for aveScan

    M.dumpModelProfiles(profilerscript) # sequential launch of aveScan
