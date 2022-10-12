#!/bin/python

import os
import shutil
import numpy as np

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

INPUTDIR="/home/innocenti/MITgcm_BFM/WORK_NWMed/V1/devel/wrkdir/MODEL/run/output0"

ALL_INDEXES = np.arange(112)
local_INDEXES = ALL_INDEXES[rank::nranks]

for it in local_INDEXES:
    print(it,rank)
    rm_dir = INPUTDIR + str(it).zfill(3)
    if os.path.isdir(rm_dir):
        shutil.rmtree(rm_dir)
#        for filename in os.listdir(rm_dir):
#            file_path = os.path.join(rm_dir , filename)
#            os.remove(file_path)
#        os.rmdir(rm_dir)
    print("rm ", rm_dir)

print("done")
