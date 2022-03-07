import argparse
def arguments():
    parser = argparse.ArgumentParser(description = '''
   Generates maskfile by reading bathymetry
    ''', formatter_class=argparse.RawTextHelpFormatter)
 
 
    parser.add_argument(   '--bathymetry','-b',
                                type = str,
                                required = True,
                                help = 'Path of the bathymetry file')
    parser.add_argument(   '--outputfile','-o',
                                type = str,
                                required = True,
                                help = 'Path of the maskfile')                                

    return parser.parse_args()

args = arguments()
import numpy as np
import scipy.io.netcdf as NC

bathyfile = args.bathymetry
maskfile  = args.outputfile



#delZ= np.array([1.500,  1.501,  3.234,  3.483,  3.750,  4.035,  4.339,  4.665,  5.012,  5.384, 
#       5.781,  6.206,  6.659,  7.144,  7.661,  8.215,  8.806,  9.437, 10.112, 10.833, 
#      11.603, 12.426, 13.305, 14.244, 15.247, 16.319, 17.463])

delZ= np.array([1.500,  1.501,  3.234,  3.477,  3.737,  4.018,  4.319,  4.643,  4.991,  5.365,
       5.768,  6.200,  6.665,  7.165,  7.703,  8.280,  8.901,  9.569, 10.287, 11.058,
      11.888, 12.779, 13.738, 14.768, 15.875, 17.066, 18.346, 19.722, 21.201, 22.791,
      24.501, 26.338, 28.314, 30.437, 32.720, 35.174, 37.812, 40.648, 43.696, 46.974,
      50.497, 54.284, 58.355, 62.732, 67.437, 72.494, 77.931, 83.776, 90.059, 96.814,
     104.075,111.881,120.272,129.292,138.989,149.413,160.619,172.665,185.615,199.537])

CellBottoms=np.cumsum(delZ)
Depth = CellBottoms - delZ/2

jpi = 784
jpj = 336
jpk = Depth.size

Lat = 41.88515625 + np.arange(jpj)*1./128
Lon = 6.06640625 + np.arange(jpi)*1./128

# tmask construction

fid=open(bathyfile,'rb')
domain_size=jpi*jpj
A=np.fromfile(fid,dtype=np.float32,count=domain_size).astype(np.float64)
fid.close()
Bathy = -A.reshape(jpj,jpi)

#ii=Bathy>0
#Bathy[ii] = Bathy[ii]-1.4e-07

LEVELS=np.zeros((jpj,jpi),np.int32)

for ji in range(jpi):
    for jj in range(jpj):
        if Bathy[jj,ji] == 0:
            LEVELS[jj,ji] = 0;
        else:
            for jk in range(jpk):
                if CellBottoms[jk] >= Bathy[jj,ji]:
                    break
            LEVELS[jj,ji]=jk+1

tmask = np.zeros((jpk,jpj,jpi), np.bool);

for ji in range(jpi):
    for jj in range(jpj):
        for jk in range(LEVELS[jj,ji]):
            tmask[jk, jj,ji] = True



ncOUT = NC.netcdf_file(maskfile,"w")
ncOUT.createDimension("lon",jpi)
ncOUT.createDimension("lat",jpj)
ncOUT.createDimension("depth",jpk)

ncvar    = ncOUT.createVariable("lon", 'f', ("lon",))
ncvar[:] = Lon
ncvar    = ncOUT.createVariable("lat", 'f', ("lat",))
ncvar[:] = Lat

ncvar    = ncOUT.createVariable("depth", 'f', ("depth",))
ncvar[:] = Depth

ncvar    = ncOUT.createVariable("CellBottoms", 'f', ("depth",))
ncvar[:] = CellBottoms

ncvar    = ncOUT.createVariable("tmask", 'b', ("depth","lat","lon"))
ncvar[:] = tmask


ncOUT.close()

