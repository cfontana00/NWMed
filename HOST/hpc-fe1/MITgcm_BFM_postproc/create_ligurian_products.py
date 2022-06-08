#!/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 13:58:04 2022

@author: Vincent FAURE (vincent.faure@ensta.org)
"""
import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Create products for geoserver NWMED : concatenate all the days of the simulation
     
    ''')
    
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                default = None,
                                required = True,
                                help ='''/some/path/  Directory containg files to interpolate. 
                                '''
                                )
    return parser.parse_args()

args = argument()


import xarray as xr
import dask
import glob
import os

# PARAMETERS
#gulf of lions
name_project='sharemod'
name_place='sharels'
longname_place='Ligurian Sea'
liste_depth=[0,1.5,10,20,50,100,1000]
#daily_hourly='sm'
daily_hourly='fc'
dir_input=args.inputdir
date_run=dir_input[-9:-1]

format_name_file="{date_of_file}_h-OGS--{name_var_file}-MITgcmBFM-pilot8-b{date_run}_{daily_hourly}-v01.nc"

# VARIABLES
list_store=[]
list_store.append(dict(name_store='pftc',name_var=['chl','zooc'],longname_var=['Modelled Chl.a ('+longname_place+' - MITgcmBFM)', \
                 'Modelled Zooplancton ('+longname_place+' - MITgcmBFM)'], \
                palette=['algae','turbid'],round_val=3,digit=3,units=['mg/m3','mmol/m3']))
list_store.append(dict(name_store='rfvl',name_var=['Speed','Direction'],longname_var=['Modelled currents speed ('+longname_place+' - MITgcmBFM)', \
                'Modelled currrents direction ('+longname_place+' - MITgcmBFM)'], \
                palette=['speed','amp'],round_val=2,digit=2,units=['m/s','degree']))
list_store.append(dict(name_store='nutr',name_var=['no3','po4','nh4','si'],longname_var=['Modelled nitrate ('+longname_place+' - MITgcmBFM)', \
                'Modelled phosphate ('+longname_place+' - MITgcmBFM)','Modelled ammonium ('+longname_place+' - MITgcmBFM)', \
                       'Modelled silica ('+longname_place+' - MITgcmBFM)'],\
                palette=['matter','matter','matter','matter'],round_val=2,digit=2,units=['mmol/m3','mmol/m3','mmol/m3','mmol/m3']))
list_store.append(dict(name_store='temp',name_var=['thetao'],longname_var=['Modelled temperature ('+longname_place+' - MITgcmBFM)'], \
                  palette=['thermal'],round_val=2,digit=3, units=['&#176;C']))
list_store.append(dict(name_store='psal',name_var=['so'],longname_var=['Modelled salinity ('+longname_place+' - MITgcmBFM)'], \
                 palette=['haline'],round_val=2,digit=2, units=['']))
    
    
for istore in list_store:
    if os.path.exists(dir_input):
        file_format=format_name_file.format(date_of_file="*",name_var_file=istore['name_store'].upper(),date_run=date_run,daily_hourly=daily_hourly)
        print(dir_input+"/"+file_format)
        liste_file=glob.glob(dir_input+"/"+file_format)
        print('liste_file',liste_file)
        if (len(liste_file))>0:
            #liste_file.sort(key=os.path.getctime)            
            liste_file.sort()
            # Merge files on time
            print('Start opening file, file 0',liste_file[0])
            df_tot=xr.open_dataset(liste_file[0],chunks={})
            df_tot=df_tot.sel(depth=liste_depth,method='nearest')
            for file in liste_file[1:]:
                print('open ',file)
                df_here=xr.open_dataset(file,chunks={})
                df_here=df_here.sel(depth=liste_depth,method='nearest')
                print('concat')
                df_tot=xr.concat([df_tot,df_here],dim='time')
            #write nc file
            outputfile=dir_input+'/'+name_project+'_'+name_place+'_'+istore['name_store']+'.nc'
            print('Start print file')
            df_tot.to_netcdf(dir_input+'/'+name_project+'_'+name_place+'_'+istore['name_store']+'.nc',encoding={'time':{'units':'days since 1900-01-01'}})        
            df_tot.close()
            print('--File created: ',dir_input+'/'+name_project+'_'+name_place+'_'+istore['name_store']+'.nc')
        
