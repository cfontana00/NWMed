import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates Netcdf from MIT .data files
    ''')


    parser.add_argument(   '--outputdir',"-o",
                                type = str,
                                required = True,
                                help = '/some/path/')

    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                default = None,
                                required = True,
                                help ='''Directory containg outputs of MITgcm
                                '''
                                )

    parser.add_argument(   '--rundate',"-d",
                                type = str,
                                default = None,
                                required = True
                                )

    parser.add_argument(   '--maskfile',"-m",
                                type = str,
                                default = None,
                                required = True
                                )
    parser.add_argument(   '--varlist',"-v",
                                type = str,
                                default = None,
                                required = True
                                )
    parser.add_argument(   '--length',"-l",
                                type = str,
                                default = '71',
                                required = False
                                )
    return parser.parse_args()

args = argument()

from commons import genUserDateList as DL 
from commons import netcdf4
from commons.mask import Mask
import numpy as np
from datetime import datetime
from commons.utils import addsep, file2stringlist



try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
except:
    rank   = 0
    nranks = 1
    isParallel = False

def readFrame_from_file(filename,Frame,shape):
    jpk,jpj,jpi=shape
    domain_size=jpi*jpj*jpk
    fid=open(filename,'rb')
    fid.seek(domain_size*Frame*4)
    A=np.fromfile(fid,dtype=np.float32,count=domain_size)
    fid.close()
    return A.reshape(jpk,jpj,jpi)


INPUTDIR=addsep(args.inputdir)
OUTDIR=addsep(args.outputdir)
RUNDATE=args.rundate
dateformat="%Y%m%d-%H:%M:%S"

fc_length=int(args.length)

rundate_dt = datetime.strptime(RUNDATE,"%Y%m%d")
#datestart = (rundate_dt - DL.relativedelta(  days=7)).strftime(dateformat)
datestart = (rundate_dt).strftime(dateformat)
#dateend   = (rundate_dt + DL.relativedelta(hours=71)).strftime(dateformat)
dateend   = (rundate_dt + DL.relativedelta(hours=fc_length)).strftime(dateformat)
#dateend   = (rundate_dt - DL.relativedelta(days=7) + DL.relativedelta(hours=24)).strftime(dateformat)
#dateend = (rundate_dt - DL.relativedelta(  days=1)).strftime(dateformat)

TheMask = Mask(args.maskfile)
VARLIST= file2stringlist(args.varlist)

timelist=DL.getTimeList(datestart, dateend, hours=1)
timestep = 120 #s, hardcoded
offset = 24*7 # hardcoded, number of outputs to skip (it depends on datestart; if datestart=rundate - 7 => offset = 0)

TimeSteps_in_h = 3600/timestep
#TimeSteps_in_h = 1
ALL_INDEXES = np.arange((len(timelist)))
local_INDEXES = ALL_INDEXES[rank::nranks]


for var in VARLIST:    
    for it in local_INDEXES:
        t = timelist[it]
        inputfile = "%s%s.%010d.data" %(INPUTDIR,var, (it+1+offset)*TimeSteps_in_h)
        outfile   = "%save.%s.%s.nc"  %(OUTDIR,t.strftime(dateformat),var)
        print(outfile)
        M3d = readFrame_from_file(inputfile, 0, TheMask.shape)
        M3d[~TheMask.mask] = 1.e+20
        netcdf4.write_3d_file(M3d, var, outfile, TheMask, compression=True)

