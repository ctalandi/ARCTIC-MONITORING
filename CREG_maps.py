#!/usr/bin/env python

import sys 
import matplotlib
matplotlib.use('Agg')
import numpy as npy
from CREG_maps_func import *
from checkfile import *
import matplotlib.pylab as plt
import matplotlib as mpl
import subprocess
from datetime import datetime
import xarray as xr 
from fsspec.implementations.local import LocalFileSystem
fs = LocalFileSystem()

s_year=XXSYEAXX
e_year=XXEYEAXX
xiosfreq=XXXIOSFREQXX
main_dir='./'
CONFIG='XXCONFXX'   ; CASE='XXCASEXX'     ;  CASE2='XXCASE2XX'
CONFCASE=CONFIG+'-'+CASE   ; CONFCASE2=CONFIG+'-'+CASE2
data_dir=main_dir+CONFIG+'/'+CONFCASE+'-MEAN/'    ; data_dir2=main_dir+CONFIG+'/'+CONFCASE2+'-MEAN/'
grid_dir=main_dir+CONFIG+'/GRID/'
obs_dir='XXOBS_DIRXX'

DIR_FIG_OUT='./'

# Set the considered period
if e_year-s_year == 0 :
        climyear=str(s_year)
else:
        climyear=str(s_year)+str(e_year)

AW_Tmax_maps=XXAW_TMAXXX
FWC_maps=XXFWC_MAPSXX
ICE_maps=XXICE_MAPSXX
MLD_maps=XXMLD_MAPSXX
DYN_maps=XXDYN_MAPSXX
TSD_maps=XXTSD_MAPSXX
ATL_maps=XXATL_MAPSXX
MOC_maps=XXMOC_MAPSXX
MTS_maps=XXMTS_MAPSXX

NCDF_OUT=XXNCDFOUTXX

print()
print('                                Configuration :' + CONFCASE)
print('                                Period        :' + str(s_year)+' - '+str(e_year))
print()


########################################
# Read GRID 
########################################
#------------------------------------------------------------------------------------------------------------------------
# Read only once for the same section type
locpath=grid_dir
locfile=CONFCASE+'_mask.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
    ds_grd = xr.open_dataset(locpath+locfile)[['glamt','gphit','tmask','fmask']]
    lon = ds_grd['glamt'].squeeze()
    lat = ds_grd['gphit'].squeeze()
    tmask = ds_grd['tmask'].isel(time_counter=0)
    tmask = tmask.rename({'nav_lev':'z'}) 
    fmask = ds_grd['fmask'].isel(time_counter=0)
    fmask = fmask.rename({'nav_lev':'z'}) 

locpath=grid_dir
locfile=CONFCASE+'_mesh_hgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
    ds_mes = xr.open_dataset(locpath+locfile)[['e1t','e2t']]
    e1t = ds_mes['e1t'].squeeze()
    e2t = ds_mes['e2t'].squeeze()

locpath=grid_dir
locfile=CONFCASE+'_mesh_zgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
    ds_zgr = xr.open_dataset(locpath+locfile)[['e3t_0','gdept_1d','gdepw_0']]
    ze3 = ds_zgr['e3t_0'].squeeze()
    ze3 = ze3.rename({'nav_lev':'z'})
    gdept1d = ds_zgr['gdept_1d'].squeeze()
    gdepw_0 = ds_zgr['gdepw_0'].squeeze()
    gdepw_0 = gdepw_0.rename({'nav_lev':'z'})

e1te2t = e1t * e2t * tmask[0,:,:]

########################################
# Read DATA 
########################################
#------------------------------------------------------------------------------------------------------------------------
# Read the whole time series
c_year=s_year
while c_year <= e_year:
        print()
        print('         The concerned year :' + str(c_year))
        print()
        locpath=data_dir+'/'+str(c_year)+'/'+xiosfreq+'/'

        #########################################################################################################################################
        if FWC_maps or ATL_maps :
             print('                    Read SSH variable')

             # List files to be read
             locfile=CONFCASE+'_y'+str(c_year)+'.'+xiosfreq+'_gridT.nc'
             SSH_files = [f for f in fs.glob(locpath+locfile)]
             
             if chkfile(locpath+locfile) :
                ds_ssh = xr.open_dataset(locpath+locfile, engine="netcdf4")[['ssh']]
                My_var1ssh = ds_ssh['ssh']
                My_var1ssh = xr.where( My_var1ssh >= 1e15, 0. , My_var1ssh )
                My_sshmean1 = (My_var1ssh * e1te2t).sum(dim=['y','x'])/e1te2t.sum(dim=['y','x'])
                My_var1ssh = My_var1ssh - My_sshmean1
                if c_year == e_year : 
                        My_var1ssh = xr.where( tmask[0,:,:] < 1, npy.nan, My_var1ssh ).squeeze()

        #########################################################################################################################################
        if ICE_maps :

             print('                    Read Ice variables ')

             zMyvar='sivolu'
             # List files to be read
             locfile=CONFCASE+'_y'+str(c_year)+'.'+xiosfreq+'_icemod.nc'
             ICE_files = [f for f in fs.glob(locpath+locfile)]
             
             if chkfile(locpath+locfile)   :      
                ds_Idata = xr.open_dataset(locpath+locfile, engine="netcdf4")[[zMyvar]]
                My_var1 = ds_Idata[zMyvar]
                if c_year == e_year : 
                        My_var1 = xr.where( tmask[0,:,:] < 1, npy.nan, My_var1 ).squeeze()
             
             zMyvar='siconc'
             # List files to be read
             locfile=CONFCASE+'_y'+str(c_year)+'m*.'+xiosfreq+'_icemod.nc'
             ICE_files = [f for f in fs.glob(locpath+locfile)]
             
             if len(ICE_files) == 12 :
                ds_Idata = xr.open_mfdataset(locpath+locfile, engine="netcdf4", concat_dim=["time_counter"], combine='nested', parallel=True)[[zMyvar]]
                My_var1frld_SeasM = ds_Idata[zMyvar].isel(time_counter=2).squeeze()
                My_var1frld_SeasS = ds_Idata[zMyvar].isel(time_counter=8).squeeze()
                if c_year == e_year : 
                        My_var1frld_SeasM = xr.where( tmask[0,:,:] < 1, npy.nan, My_var1frld_SeasM ).squeeze()
                        My_var1frld_SeasS = xr.where( tmask[0,:,:] < 1, npy.nan, My_var1frld_SeasS ).squeeze()

        #########################################################################################################################################
        if MLD_maps or MTS_maps or ATL_maps :

             print('                         Read Mixed layer depth variable  ')

             zMyvar='mldr10_1'
             # List files to be read
             locfile=CONFCASE+'_y'+str(c_year)+'m*.'+xiosfreq+'_gridT.nc'
             TS_files = [f for f in fs.glob(locpath+locfile)]
             
             if len(TS_files) == 12 :
                ds_TSdata = xr.open_mfdataset(locpath+locfile, engine="netcdf4", concat_dim=["time_counter"], combine='nested', parallel=True)[[zMyvar,'votemper','vosaline','deptht']]
                ds_TSdata = ds_TSdata.rename({'deptht':'z'})
                # Keep the MLD only in March and September
                My_var1SeasM = ds_TSdata[zMyvar].isel(time_counter=2).squeeze()
                My_var1SeasS = ds_TSdata[zMyvar].isel(time_counter=8).squeeze()
                # Keep the MLD over the 12 months 
                Mdata_read = ds_TSdata[zMyvar]
                # Keep both temperature/salinity in March and September
                My_varTSM = ds_TSdata.isel(time_counter=2).squeeze() 
                My_varTSS = ds_TSdata.isel(time_counter=8).squeeze() 

                if c_year == e_year : 
                        My_var1SeasM = xr.where( tmask[0,:,:] < 1, npy.nan, My_var1SeasM ).squeeze()
                        My_var1SeasS = xr.where( tmask[0,:,:] < 1, npy.nan, My_var1SeasS ).squeeze()
                        My_varTSM = xr.where( tmask < 1, npy.nan, My_varTSM ).squeeze()
                        My_varTSS = xr.where( tmask < 1, npy.nan, My_varTSS ).squeeze()

        #########################################################################################################################################
        if DYN_maps :

             zMyvar = 'sobarstf' 
             print('                    Read PSI variable ')

             # List files to be read
             # Read annual mean
             locfile=CONFCASE+'_y'+str(c_year)+'.'+xiosfreq+'_PSI.nc'
             if chkfile(locpath+locfile) :
                ds_psi = xr.open_dataset(locpath+locfile, engine="netcdf4")[[zMyvar]]
                My_var1 = ds_psi[zMyvar]
                if c_year == e_year : 
                        My_var1 = xr.where( fmask[0,:,:] < 1, npy.nan, My_var1 ).squeeze()

             zMyvar = 'voeke' 
             zd1 = npy.round(gdept1d.isel(nav_lev=0).values,2).item()  ; zd2 = npy.round(gdept1d.isel(nav_lev=23).values,2).item()
             print('                    Read EKE variable @ depth ' + str(zd1) + ' m & ' + str(zd2) + ' m' )
             
             # List files to be read
             # Read annual mean
             locfile=CONFCASE+'_y'+str(c_year)+'.'+xiosfreq+'_EKE.nc'
             if chkfile(locpath+locfile) : 
                ds_eke = xr.open_dataset(locpath+locfile, engine="netcdf4")[[zMyvar]]
                My_var1SeasM = ds_eke[zMyvar].isel(z=0)
                My_var1SeasS = ds_eke[zMyvar].isel(z=23)
                if c_year == e_year : 
                        My_var1SeasM = xr.where( tmask[0,:,:]  < 1, npy.nan, My_var1SeasM ).squeeze()
                        My_var1SeasS = xr.where( tmask[23,:,:] < 1, npy.nan, My_var1SeasS ).squeeze()
             else : 
                My_var1SeasM = xr.full_like(My_var1, fill_value=npy.nan)
                My_var1SeasS = xr.full_like(My_var1, fill_value=npy.nan)

        #########################################################################################################################################
        if ( AW_Tmax_maps or FWC_maps or TSD_maps or ATL_maps ) :

                   print('                      Read full T/S 3D spatial variables to compute the AWTmax or FWC ')

                   # List files to be read
                   # Read annual mean 
                   locfile=CONFCASE+'_y'+str(c_year)+'.'+xiosfreq+'_gridT.nc'
                   if chkfile(locpath+locfile) :
                          ds_TSdata = xr.open_dataset(locpath+locfile, engine="netcdf4")[['votemper','vosaline','deptht']]
                          ds_TSdata = ds_TSdata.rename({'deptht':'z'})
                          My_var1T = ds_TSdata['votemper'].squeeze()
                          My_var1S = ds_TSdata['vosaline'].squeeze()
                          if c_year == e_year :
                                  My_var1T = xr.where( tmask  < 1, npy.nan, My_var1T ).squeeze()
                                  My_var1S = xr.where( tmask  < 1, npy.nan, My_var1S ).squeeze()

        #########################################################################################################################################
        if MOC_maps :

                   print('                      Read full MOC 2D spatial variables ')

                   # Read annual mean 
                   locfile=CONFCASE+'_y'+str(c_year)+'.'+xiosfreq+'_MOC.nc'
                   if chkfile(locpath+locfile) : 
                        ds_moc = xr.open_dataset(locpath+locfile, engine="netcdf4")[['zomsfglo']]
                        My_MOC = ds_moc['zomsfglo'].squeeze()
        
        c_year = c_year + 1 

# Read the initial state to compare with
if AW_Tmax_maps or FWC_maps or TSD_maps or ATL_maps : My_varTinit, My_varSinit = CREG_INIT( CONFIG, CASE )

###########################################
# CALL FUNCTIONS TO PERFORM DIAGS AND PLOTS 
###########################################
#------------------------------------------------------------------------------------------------------------------------
# In the following call the appropriate function depending on the diagnostic

# To plot the mean T/S in the ML
if MTS_maps : 
	MTS_maps( lon, lat, CONFIG, CASE, My_var1SeasM, My_var1SeasS, My_varTSM, My_varTSS, gdepw_0, ze3, climyear, data_dir, grid_dir, NCDF_OUT )

# To plot the Atlantic Water maximum temperature as the associated depth
# Use the salinity criteria S < 33.5
if AW_Tmax_maps : 
        zMyvar='votemper'
        AWT_maps( lon, lat, My_var1T, My_var1S, gdept1d, zMyvar, CONFIG, CASE, climyear, NCDF_OUT )

# To plot SSH and FWC (based on a salinity ref of 34.8 PSU)
if FWC_maps : 
        FWC_maps( lon, lat, My_var1S, My_varSinit, My_var1ssh, CONFIG, CASE, climyear, ze3, tmask, NCDF_OUT )

# To plot ICE variables
if ICE_maps : 
        ICE_maps( lon, lat, My_var1, My_var1frld_SeasM, My_var1frld_SeasS, CONFIG, CASE, climyear, s_year, c_year, NCDF_OUT )

# To plot MLD variable
if MLD_maps : 
        MLD_maps( lon, lat, My_var1SeasM, My_var1SeasS, CONFIG, CASE, climyear, NCDF_OUT )

# To plot DYN variables PSI and EKE 
if DYN_maps : 
        DYN_maps( lon, lat, My_var1, My_var1SeasM, My_var1SeasS, gdept1d, CONFIG, CASE, climyear, s_year, NCDF_OUT )

# To plot T/S drift at the surface, ~100m, ~200m & ~300m
if TSD_maps : 
        TSD_maps( lon, lat, My_var1T, My_var1S, My_varTinit, My_varSinit, gdept1d, CONFIG, CASE, climyear, NCDF_OUT )

# To plot ATL variables such as MLD, SSH in differents areas GIN, LAB, IRM seas 
if ATL_maps : 
        ATL_maps( lon, lat, My_var1SeasM, Mdata_read, My_var1T, My_varTinit, My_var1ssh, gdept1d, CONFIG, CASE, climyear, s_year, e_year, NCDF_OUT )

# To plot AMOC and Time series 
if MOC_maps : 
        MOC_maps( lon, lat, My_MOC, gdept1d, CONFIG, CASE, climyear, s_year, e_year, NCDF_OUT )
