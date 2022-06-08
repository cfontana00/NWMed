#!/opt/rh/rh-python36/root/usr/bin/python
# read STDOUT and plot infomation

import matplotlib.pyplot as plt
import numpy as np

from mds import rdmds
from utils import addsep

T,its,meta=rdmds('output0000/T',18,machineformat='l',returnmeta=True)
XC=rdmds('output*/XC',machineformat='l')

print(XC.shape)
fig,ax=plt.subplots(1,1)
cp = ax.contourf(T[0][:][:])
cbar=fig.colorbar(cp) 
#plt.show()
plt.savefig('fig1.png')
OUTDIR='output'
OUTPUTDIR=addsep(OUTDIR)
print(OUTPUTDIR)
