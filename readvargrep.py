import os
import pandas as pd
import re

def readvargrep(var,filename):

    varname_dynstat=['eta','uvel','vvel','wvel','theta','salt','sst','sss']
    varname_forcing=['qnet','qsw','empmr','fu','fv']
    varname_ke=['ke']
    varname_exf=['ustress','vstress','hflux','sflux','uwind','vwind',
                 'wspeed','atemp','aqh','lwflux','precip','swflux','evap',
                 'swdown','lwdown','apressure','runoff']
    varname_trc=['ptracer14','ptracer19','ptracer23','ptracer27',
                 'ptracer01','ptracer02','ptracer03','ptracer04',
                 'ptracer52','ptracer53','ptracer54','ptracer55',
                 'ptracer56','ptracer57','ptracer58','ptracer59',
                 'ptracer60','ptracer61','ptracer62','ptracer63']

    if any(var in s for s in varname_dynstat): # checkvar physics
        vargrep='dynstat_'+var+'_mean' # var physics
        # vargrep=['dynstat_' var '_min']; # var physics
        # vargrep=['dynstat_' var '_max']; # var physics
    elif any(var in s for s in varname_forcing): # checkvar forcing
        vargrep='forcing_'+var+'_mean' # var forcing
    elif any(var in s for s in varname_ke): # checkvar ke
        vargrep=var+'_mean' # var ke
    elif any(var in s for s in varname_exf): # checkvar exf
        vargrep='exf_'+var+'_mean' # var exf
    # vargrep=['exf_' var '_max']; # var exf
    elif any(var in s for s in varname_trc):  # checkvar tracer
        vargrep='trcstat_'+var+'_mean' # var tracer
    # vargrep=['trcstat_' var '_max']; # var tracer

    #command=''grep ' vargrep ' ' filename ' | awk ''{print $6}'' > tmp.txt''
    f = open('tmp.txt','w')
    data=[]
    with open(filename, 'r') as f2:
        for line in f2.readlines():
            if re.search(vargrep,line):
                line_i = line.split()
                #print(line)
                f.write(line_i[5]+'\n')
                data.append(float(line_i[5]))
    f.close()
    #os.system('grep vargrep STDOUT.0000_g100 > temp.txt')
    #os.system('awk ''{print $6}'' temp.txt > tmp.txt')
    #data = pd.read_csv('tmp.txt')
    return data

