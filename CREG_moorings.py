#!/usr/bin/env python
#### Enthought Python Distribution (free for Academic) recommended

import matplotlib
matplotlib.use('Agg')
import sys
import matplotlib.pylab as plt
import numpy as npy
from checkfile import *
from CREG_moorings_cont import *
from CREG_moorings_func import *
import xarray as xr
from fsspec.implementations.local import LocalFileSystem
fs = LocalFileSystem()

s_year=XXSYEAXX
e_year=XXEYEAXX
lgTS_ys=XXLGTSSXX
lgTS_ye=XXLGTSEXX
xiosfreq=XXXIOSFREQXX
NCDF_OUT=XXNCDFOUTXX
m_alpha=1.
main_dir='./'
CONFIG='XXCONFXX'   ; CASE='XXCASEXX'     ;  CASE2='XXCASE2XX'
CONFCASE=CONFIG+'-'+CASE
data_dir=main_dir+CONFIG+'/'+CONFCASE+'-MEAN/'
grid_dir=main_dir+CONFIG+'/GRID/'

DIR_FIG_OUT='./'

monthly_data = 1  # Read monthly mean field for temporal evolution 
depth_time_plot=True   # Do the plot or not 
do_Kprofile=True     # Read a previous profile to compare the current one

########################################
# Read GRID 
########################################
#------------------------------------------------------------------------------------------------------------------------
# Read only once for the same section type
print()
print()
print('      >>>>>>>>>>>>>>>>>>>>>>>>>>>>Perform the treatment of the '+CONFCASE+' experiment ')
print()
print()
print('         Start reading geographical coordinates, scale factors & masks ')
print()

locpath=grid_dir
locfile=CONFCASE+'_mask.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
    ds_grd=xr.open_dataset(locpath+locfile)
    tmask = ds_grd['tmask'].isel(time_counter=0)
    tmask = tmask.rename({'nav_lev':'z'}) 

locpath=grid_dir
locfile=CONFCASE+'_mesh_hgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
    ds_mes=xr.open_dataset(locpath+locfile)
    lon = ds_mes['glamt'].squeeze()
    lat = ds_mes['gphit'].squeeze()
    e1t = ds_mes['e1t'].squeeze()
    e2t = ds_mes['e2t'].squeeze()

locpath=grid_dir
locfile=CONFCASE+'_mesh_zgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
    ds_zgr=xr.open_dataset(locpath+locfile)
    e3t = ds_zgr['e3t_0'].squeeze()
    e3t = e3t.rename({'nav_lev':'z'})

name_z='deptht'
locpath=data_dir+'/'+str(s_year)+'/'
locfile=CONFCASE+'_y'+str(s_year)+'.'+xiosfreq+'_gridT.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
    ds_fld=xr.open_dataset(locpath+locfile)
    z1D = ds_fld[name_z]
z2dt=npy.reshape(z1D.values,(z1D.size,1))
print()
print('         End reading geographical coordinates, scale factors & masks ')

# Define time axis
if monthly_data == 1:
        t_months=(npy.arange(12)*30.+15.)/365.
        c_year=s_year   ; start = 1
        while c_year <= e_year:
                y_years=npy.tile(c_year,12)+t_months
                if start == 1:
                        time_axis=y_years
                        start=0
                else:
                        time_axis=npy.append(time_axis,y_years)
                c_year+=1
else:
        time_axis=npy.arange(e_year-s_year+1)+s_year+0.5

if monthly_data == 1 and e_year-s_year == 0:
    time_axis=npy.arange(12)+1.5

# 2D fields to use for the plotting step 
xplt= npy.tile(time_axis,(z1D.size,1))
zplt = npy.repeat(z2dt,time_axis.shape[0],axis=1)

#------------------------------------------------------------------------------------------------------------------------
########################################
# Read input  data 
########################################
c_year=s_year
while c_year <= e_year:
      # List files to be read
      locpath=data_dir+'/'+str(s_year)+'/'
      locfile=CONFIG+'-'+CASE+'_y'+str(s_year)+'m??.'+xiosfreq+'_gridT.nc'
      TS_files = [f for f in fs.glob(locpath+locfile)]
      
      drp_var=["time_centered", "deptht_bounds","time_centered_bounds","time_counter_bounds"]
      if len(TS_files) == 12 :
         print()
         print('  Start to read T & S data ')
         print() 
         ds_TSdata = xr.open_mfdataset(locpath+locfile, engine="netcdf4", concat_dim=["time_counter"], combine='nested', parallel=True, drop_variables=drp_var)
         ds_TSdata = ds_TSdata.rename({'deptht':'z'})
         ## Fill NaN with zero
         #Sdata_read = Sdata_read.fillna(0)
         ## Mask land grid points
         #Sdata_read = Sdata_read * tmask

      c_year = c_year + 1 

# Read the initial state 
locpath=data_dir+'/'
locfile=CONFCASE+'_init_gridT.nc'
if chkfile(locpath+locfile) :
    print()
    print('  Start to read initial T & S data ')
    print() 
    drp_var=["time_centered", "deptht_bounds","time_centered_bounds","time_counter_bounds"]
    ds_TSinit = xr.open_dataset(locpath+locfile, engine="netcdf4", drop_variables=drp_var)

#------------------------------------------------------------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
###########################################   LONG TIME SERIES PROFILES   ##########################################
####################################################################################################################
####################################################################################################################
#if False: 
if lgTS_ye-lgTS_ys+1 > 1 :
    plt.clf()
    plt.subplots_adjust(bottom=0.1, left=0.1, right=0.95, top=0.9)
    num_fram=410
    All_box=XXALLBOXXX
    All_var=DEF_INFO_VAR('gridT','votemper')
    cumul_box=''
    moor_n=1
    for zmybox in All_box :
        zbox=DEF_MOOR_BOX(CONFIG,zmybox)
        hsct_lev,Red_My_varinit=DEF_REDVAR(CONFIG,CASE,ds_TSdata,ds_TSinit,zbox,All_var,s_year,e_year,tmask,e1t,e2t,e3t)
        fram=num_fram+moor_n
        DEF_ZTIME(CONFIG,CASE,lgTS_ys,lgTS_ye,zbox,z1D,z2dt,hsct_lev,zfram=fram,zoutNC=NCDF_OUT) 
        cumul_box=cumul_box+'-'+zbox['box']
        moor_n+=2

    # Save year time series to use for a longer time series
    timeadd='_month' if monthly_data == 1 else ''

    file_ext1='_TemSalEvol'+timeadd
    file_ext=file_ext1+cumul_box+'_LGTS_y'+str(lgTS_ys)+'LASTy'
    file_name=CONFCASE+file_ext

    plt.tight_layout()
    plt.savefig(DIR_FIG_OUT+file_name+'.png',dpi=300)


#------------------------------------------------------------------------------------------------------------------------
####################################################################################################################
####################################################################################################################
###########################################   VERTICAL PROFILES   ##################################################
####################################################################################################################
####################################################################################################################
if do_Kprofile :
    print()
    print('             ##################################################################  ')
    print('             ##################################################################  ')
    print('             ######### PLOT VERTICAL PROFILES AT MOORINGS LOCATION ############  ')
    print('             ##################################################################  ')
    print('             ##################################################################  ')
    print()
    plt.clf()
    num_fram=220
    All_box=XXALLBOXXX
    All_var=DEF_INFO_VAR('gridT','votemper')
    cumul_box=''
    moor_n=1
    for zmybox in All_box :
        zbox=DEF_MOOR_BOX(CONFIG,zmybox)
        fram=num_fram+moor_n
        hsct_lev,Red_My_varinit=DEF_REDVAR(CONFIG,CASE,ds_TSdata,ds_TSinit,zbox,All_var,s_year,e_year,tmask,e1t,e2t,e3t)
        DEF_ZPROFILE(CONFIG,CASE,str(s_year),hsct_lev,Red_My_varinit,zbox,All_var,zplt,zfram=fram,zs_year=s_year) 
        cumul_box=cumul_box+'-'+zbox['box']
        moor_n+=2

    pfile_ext='_TemSaliProfile'+cumul_box+'_y'+str(s_year)+'_z0-500'
    file_name=CONFIG+'-'+CASE+pfile_ext
    plt.tight_layout()
    plt.savefig(DIR_FIG_OUT+file_name+'.png',dpi=300)
