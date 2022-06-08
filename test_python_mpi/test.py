import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates netCDF4 compressed files

    Parallel executable, can be called by mpirun.    
   ''',formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(   '--filelist',"-l",
                                type = str,
                                default = "*.nc",
                                help = 'ave*.N1p.nc')
    
    return parser.parse_args()


args = argument()

import os
import glob

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nranks = comm.size
except:
    rank = 0
    nranks = 1

PATH_NAME = args.filelist
FILELIST=glob.glob(PATH_NAME)
#FILELIST.sort
print(FILELIST)

for filename in FILELIST[rank::nranks]:
    basename=os.path.basename(filename)
    outfile= basename
    print(outfile,rank,nranks)
