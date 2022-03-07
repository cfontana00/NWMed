#!/opt/rh/rh-python36/root/usr/bin/python
# read STDOUT and plot infomation

from readvargrep import readvargrep
import matplotlib.pyplot as plt
import numpy as np

nd=10
tlim=nd*24
tim = np.linspace(0,240,num=240)

# read STDOUT
filename='wrkdir/MODEL/run/STDOUT.0000'
#filename2='/home/innocenti/MITgcm_BFM/products/6_nodi_nodebug/STDOUT.0000'
filename2='/home/innocenti/MITgcm_BFM/WORK/V1/devel/wrkdir/MODEL/run/STDOUT.0000'

list1=['eta','uvel','vvel','wvel','theta','salt','sst','sss']
plt.figure(1)
for i in range(0,8):
    var=list1[i]
    tmp=readvargrep(var,filename)
    tmp2=readvargrep(var,filename2)
    plt.subplot(4,2,i+1)
    plt.plot(tim,tmp[0:tlim], label=var+'LON=5')
    plt.plot(tim,tmp2[0:tlim], label=var+'LON=3')
    plt.legend()
    plt.savefig('fig1.png')

list2=['qnet','qsw','empmr']
plt.figure(2)
for i in range(0,3):
    var=list2[i]
    tmp=readvargrep(var,filename)
    tmp2=readvargrep(var,filename2)
    plt.subplot(3,1,i+1)
    plt.plot(tim,tmp[0:tlim], label=var+'LON=5')
    plt.plot(tim,tmp2[0:tlim], label=var+'LON=3')
    plt.legend()
    plt.savefig('fig2.png')

list3=['ke']
plt.figure(3)
for i in range(0,1):
    var=list3[i]
    tmp=readvargrep(var,filename)
    tmp2=readvargrep(var,filename2)
    plt.subplot(1,1,i+1)
    plt.plot(tim,tmp[0:tlim], label=var+'LON=5')
    plt.plot(tim,tmp2[0:tlim], label=var+'LON=3')
    plt.legend()
    plt.savefig('fig3')

#% list4={'hflux','sflux','lwflux','precip','swflux','evap','ustress','vstress'};
#% list4={'hflux','sflux','lwdown','precip','swdown','wspeed','atemp','aqh'};
list4=['wspeed','hflux','sflux','lwflux','swflux','atemp','evap','precip']
plt.figure(4)
for i in range(0,8):
    var=list4[i]
    tmp=readvargrep(var,filename)
    tmp2=readvargrep(var,filename2)
    plt.subplot(4,2,i+1)
    plt.plot(tim,tmp[0:tlim], label=var+'LON=5')
    plt.plot(tim,tmp2[0:tlim], label=var+'LON=3')
    plt.legend()
    plt.savefig('fig4')

#plt.show(block=False)
#plt.pause(0.001)         
#input("hit[enter] to end.")                                                               
#plt.close('all') 

list5=['ptracer01','ptracer02','ptracer03','ptracer04'];
plt.figure(5)
for i in range(0,4):
    var=list5[i]
    tmp=readvargrep(var,filename)
    tmp2=readvargrep(var,filename2)
    plt.subplot(3,2,i+1)
    plt.plot(tim,tmp[0:tlim], label=var+'LON=5')
#    plt.plot(tim,tmp2[0:tlim], label=var+'LON=3')
    plt.legend()
    plt.savefig('fig6')

list6=['trAdv_CFL_u_max','trAdv_CFL_v_max','trAdv_CFL_w_max'];
plt.figure(6)
for i in range(0,3):
    var=list6[i]
    tmp=readvargrep(var,filename)
    tmp2=readvargrep(var,filename2)
    plt.subplot(3,1,i+1)
    plt.plot(tim,tmp[0:tlim], label=var+'LON=5')
    plt.plot(tim,tmp2[0:tlim], label=var+'LON=3')
    plt.legend()
    plt.savefig('fig7')


#% list5={'ptracer14','ptracer19','ptracer23','ptracer27'};
#% list5={'ptracer52','ptracer53','ptracer54','ptracer55',...
#%     'ptracer56','ptracer57','ptracer58','ptracer59',...
#%     'ptracer60','ptracer61','ptracer62','ptracer63'};
#figure
#for i=1:1:length(list5)
#    var=list5{i};     tmp=readvargrep(var,filename);
#    subplot(2,2,i); plot(tmp(1:tlim),'-k'); hold on; legend(var,0);
#    % subplot(4,3,i); plot(tmp,'-k'); hold on; legend(var,0);
#end
#%}
