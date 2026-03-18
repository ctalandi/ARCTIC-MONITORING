#!/usr/bin/env python

import sys 
import matplotlib
matplotlib.use('Agg')
import numpy as npy
import CREG_maps_func
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

# Infos concernant les climatologies sur la periode de la simulation
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

time_dim=e_year-s_year+1

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
    ds_grd=xr.open_dataset(locpath+locfile)[['glamt','gphit','tmask','fmask']]
    lon = ds_grd['glamt'].squeeze()
    lat = ds_grd['gphit'].squeeze()
    tmask = ds_grd['tmask'].isel(time_counter=0)
    tmask = tmask.rename({'nav_lev':'z'}) 
    fmask = ds_grd['fmask'].isel(time_counter=0)
    fmask = fmask.rename({'nav_lev':'z'}) 

locpath=grid_dir
locfile=CONFCASE+'_mesh_hgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
    ds_mes=xr.open_dataset(locpath+locfile)
    e1t = ds_mes['e1t'].squeeze()
    e2t = ds_mes['e2t'].squeeze()

locpath=grid_dir
locfile=CONFCASE+'_mesh_zgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
    ds_zgr=xr.open_dataset(locpath+locfile)[['e3t_0','gdept_1d','gdepw_0']]
    ze3 = ds_zgr['e3t_0'].squeeze()
    ze3 = ze3.rename({'nav_lev':'z'})
    z = ds_zgr['gdept_1d'].squeeze()
    gdepw_0 = ds_zgr['gdepw_0'].squeeze()
    gdepw_0 = gdepw_0.rename({'nav_lev':'z'})

e1te2t = e1t * e2t * tmask[0,:,:]

########################################
# Read DATA 
########################################
#------------------------------------------------------------------------------------------------------------------------
time_dim=(e_year-s_year+1)*12
t_months=(npy.arange(12)*30.+15.)/365.   ;   start = 1
time_grid=npy.arange(s_year,e_year+2,1.,dtype=int)
newlocsx  = npy.array(time_grid,'f')
newlabelsx = npy.array(time_grid,'i')

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

             # Set the time axis
             y_years=npy.tile(c_year,12)+t_months
             if start == 1:
                     time_axis=y_years
                     start=0
             else:
                     time_axis=npy.append(time_axis,y_years)

        #########################################################################################################################################
        if DYN_maps :

             zMyvar = 'sobarstf' 
             print('                    Read PSI variable ')

             # List files to be read
             # Read annual mean
             locfile=CONFCASE+'_y'+str(c_year)+'.'+xiosfreq+'_PSI.nc'
             drp_var=["time_centered", "deptht_bounds","time_centered_bounds","time_counter_bounds"]
             if chkfile(locpath+locfile) :
                ds_psi = xr.open_dataset(locpath+locfile, engine="netcdf4")[[zMyvar]]
                My_var1 = ds_psi[zMyvar]
                if c_year == e_year : 
                        My_var1 = xr.where( fmask[0,:,:] < 1, npy.nan, My_var1 ).squeeze()

             zMyvar = 'voeke' 
             zd1 = npy.round(z.isel(nav_lev=0).values,2).item()  ; zd2 = npy.round(z.isel(nav_lev=23).values,2).item()
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
if AW_Tmax_maps or FWC_maps or TSD_maps or ATL_maps :
   locpath=data_dir+'/'
   locfile=CONFCASE+'_init_gridT.nc'
   if chkfile(locpath+locfile) :
       print('                      Read initial state  ')
       ds_TSinit = xr.open_dataset(locpath+locfile, engine="netcdf4")[['votemper','vosaline']]
       ds_TSinit = ds_TSinit.rename({'nav_lev':'z'})
       My_varTinit = ds_TSinit['votemper'].squeeze()
       My_varSinit = ds_TSinit['vosaline'].squeeze()

###########################################
# Call FUNCTIONS TO PERFORM DIAGS AND PLOTS 
###########################################
#------------------------------------------------------------------------------------------------------------------------
# In the following call the appropriate function depending on the diagnostic

# To plot the mean T/S in the ML
if MTS_maps : CREG_maps_func.MTS_maps(lon, lat, CONFIG, CASE, My_var1SeasM, My_var1SeasS, My_varTSM, My_varTSS, gdepw_0, ze3, climyear, data_dir, grid_dir, NCDF_OUT )

# To plot the Atlantic Water maximum temperature as the associated depth
# Use the salinity criteria S < 33.5
if AW_Tmax_maps : 
        zMyvar='votemper'
        CREG_maps_func.AWTmax_maps( lon, lat, My_var1T, My_var1S, z, zMyvar, CONFIG, CASE, climyear, NCDF_OUT)

# To plot SSH and FWC (based on a salinity ref of 34.8 PSU)
if FWC_maps : 
        CREG_maps_func.FWC_maps( lon, lat, My_var1S, My_varSinit, My_var1ssh, CONFIG, CASE, climyear, ze3, tmask, time_dim, NCDF_OUT )

# To plot ICE variables
if ICE_maps : 
        num_fram=320
        # Annual mean Ice thickness
        zMyvar='sivolu'   ; fram=num_fram+1
        My_var1 = xr.where( My_var1 == 0., npy.nan, My_var1 )
        #My_var1=npy.ma.masked_where(My_var1 == 0., My_var1).squeeze()
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, My_var1, zMyvar, climyear, zfram=fram)
        # March mean Ice fraction 
        zMyvar='siconc'   ; fram=num_fram+3
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, My_var1frld_SeasM, zMyvar, climyear, seas='m03', zfram=fram )
        # September mean Ice fraction 
        zMyvar='siconc'   ; fram=num_fram+5
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, My_var1frld_SeasS, zMyvar, climyear, seas='m09', zfram=fram )

        # Annual mean Ice thickness PIOMASS observations
        # WARNING this part has been interpolated on CREG025 grid directly 
        # Read lat, lon from CREG025 grid since the PIOMASS data is not available on CREG12.L75 grid
        if CONFIG == 'CREG12.L75' : 
                locpath=grid_dir
                locfile='CREG025.L75_mask.nc'
                if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
                        ds_msk=xr.open_dataset(locpath+locfile)[['glamt','gphit','tmask']]
                        lon = ds_msk['glamt'].squeeze()
                        lat = ds_msk['gphit'].squeeze()
                        zmask= ds_msk['tmask'].squeeze()
        else: 
                zmask = tmask.copy()

        zMyvar='sivolu'   ; fram=num_fram+2
        obs_thick = CREG_maps_func.ICE_THICK_OBS(zconfig=CONFIG,t_year=s_year)
        obs_thick = xr.where( npy.squeeze(zmask[0,:,:]) < 1., npy.nan, obs_thick)
        obs_thick = xr.where( obs_thick == 0., npy.nan, obs_thick)
        #obs_thick = npy.ma.masked_where(npy.squeeze(zmask[0,:,:]) == 0., obs_thick)
        #obs_thick = npy.ma.masked_where(obs_thick == 0., obs_thick)
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, obs_thick, zMyvar, s_year, zfram=fram, plot_obs=1)

        obs_conc_m03,obs_conc_m09,obs_lon,obs_lat=CREG_maps_func.ICE_CONCE_OBS(t_year=c_year)   
        # March mean Ice fraction 
        zMyvar='siconc'   ; fram=num_fram+4
        CREG_maps_func.simple_maps( obs_lon, obs_lat, CONFIG, CASE, obs_conc_m03, zMyvar, s_year, seas='m03', zfram=fram, plot_obs=1 )
        # September mean Ice fraction 
        zMyvar='siconc'   ; fram=num_fram+6
        CREG_maps_func.simple_maps( obs_lon, obs_lat, CONFIG, CASE, obs_conc_m09, zMyvar, s_year, seas='m09', zfram=fram, plot_obs=1 )

        zfile_ext='_ICEClim_'
        plt.tight_layout()
        plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.png',dpi=300)

        if NCDF_OUT:
                ds_out = xr.Dataset()
                
                # ICE fields
                #######################
                ds_out['IceThick_mod']= (('y','x'), My_var1.values.astype('float32')) 
                ds_out['IceThick_mod'].attrs['long_name']='Model annual mean Ice thickness'
                ds_out['IceThick_mod'].attrs['units']='m'
                
                ds_out['IceConceM03_mod']= (('y','x'), My_var1frld_SeasM.values.astype('float32')) 
                ds_out['IceConceM03_mod'].attrs['long_name']='Model monthly mean Ice concentration in march'
                ds_out['IceConceM03_mod'].attrs['units']='-'
                
                ds_out['IceConceM09_mod']= (('y','x'), My_var1frld_SeasS.values.astype('float32')) 
                ds_out['IceConceM09_mod'].attrs['long_name']='Model monthly mean Ice concentration in september'
                ds_out['IceConceM09_mod'].attrs['units']='-'
                
                ds_out['IceThick_obs']= (('y','x'), obs_thick.values.astype('float32')) 
                if s_year >= 1979 and s_year <= 2013 :
                        ds_out['IceThick_obs'].attrs['long_name']='PIOMAS annual mean Ice thickness over '+str(s_year)
                else :
                        ds_out['IceThick_obs'].attrs['long_name']='PIOMAS climatological mean Ice thickness over 1979-2013'
                ds_out['IceThick_obs'].attrs['units']='m'
                
                ds_out['IceConceM03_obs']= (('yobs','xobs'), obs_conc_m03.values.astype('float32')) 
                if s_year >= 1979 and s_year <= 2015 :
                        ds_out['IceConceM03_obs'].attrs['long_name']='NSDIC monthly mean Ice concentration in march '+str(s_year)
                else :
                        ds_out['IceConceM03_obs'].attrs['long_name']='NSDIC climatological mean Ice concentration in march over 1979-2015'
                ds_out['IceConceM03_obs'].attrs['units']='-'
                
                ds_out['IceConceM09_obs']= (('yobs','xobs'), obs_conc_m09.values.astype('float32')) 
                if s_year >= 1979 and s_year <= 2015 :
                        ds_out['IceConceM09_obs'].attrs['long_name']='NSDIC monthly mean Ice concentration in september '+str(s_year)
                else :
                        ds_out['IceConceM09_obs'].attrs['long_name']='NSDIC climatological mean Ice concentration in september over 1979-2015'
                ds_out['IceConceM09_obs'].attrs['units']='-'
                
                ds_out['lat_obs']= (('yobs','xobs'), obs_lat.values.astype('float32')) 
                ds_out['lat_obs'].attrs['long_name']='Degrees north'
                ds_out['lat_obs'].attrs['units']='Deg'
                
                ds_out['lon_obs']= (('yobs','xobs'), obs_lon.values.astype('float32')) 
                ds_out['lon_obs'].attrs['long_name']='Degrees east'
                ds_out['lon_obs'].attrs['units']='Deg'
                
                ds_out['lat_mod']= (('y','x'), lat.values.astype('float32')) 
                ds_out['lat_mod'].attrs['long_name']='Degrees north'
                ds_out['lat_mod'].attrs['units']='Deg'
                
                ds_out['lon_mod']= (('y','x'), lon.values.astype('float32')) 
                ds_out['lon_mod'].attrs['long_name']='Degrees east'
                ds_out['lon_mod'].attrs['units']='Deg'

                ds_out = ds_out.set_coords(['lat_obs','lon_obs','lat_mod','lon_mod'])

                # Write the NetCDF file 
                ds_out.attrs['History'] = 'Diagnostics have been calculated using the Arctic monitoring tool '
                ds_out.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
                nc_f = './NETCDF/'+CONFIG+'-'+CASE+'_ICEClim_'+'y'+climyear+'.nc'
                ds_out.to_netcdf(nc_f,engine='netcdf4')

# To plot MLD variable
if MLD_maps : 
        num_fram=220
        # March mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+1
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, My_var1SeasM, zMyvar, climyear, seas='m03', zfram=fram )
        # September mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+3
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, My_var1SeasS, zMyvar, climyear, seas='m09', zfram=fram )

        # MLD from observation
        mld_obs_m03, mld_obs_m09, lon_obs, lat_obs = CREG_maps_func.MLD_OBS()
        # March mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+2
        CREG_maps_func.simple_maps( lon_obs, lat_obs, CONFIG, CASE, mld_obs_m03, zMyvar, climyear, seas='m03', zfram=fram, plot_obs=1 )
        # September mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+4
        CREG_maps_func.simple_maps( lon_obs, lat_obs, CONFIG, CASE, mld_obs_m09, zMyvar, climyear, seas='m09', zfram=fram, plot_obs=1 )

        zfile_ext='_MLDClim_'
        plt.tight_layout()
        plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.png',dpi=300)

        if NCDF_OUT:
                ds_out = xr.Dataset()
                
                # MLD fields
                #######################
                ds_out['MLDd001M03_mod']= (('y','x'), My_var1SeasM.values.astype('float32')) 
                ds_out['MLDd001M03_mod'].attrs['long_name']='Model monthly mean MLD in march based on a density criteria of 0.01 kg/m^3'
                ds_out['MLDd001M03_mod'].attrs['units']='m'

                ds_out['MLDd001M09_mod']= (('y','x'), My_var1SeasS.values.astype('float32')) 
                ds_out['MLDd001M09_mod'].attrs['long_name']='Model monthly mean MLD in september based on a density criteria of 0.01 kg/m^3'
                ds_out['MLDd001M09_mod'].attrs['units']='m'

                ds_out['MLDd001M03_obs']= (('yobs','xobs'), mld_obs_m03.values.astype('float32')) 
                ds_out['MLDd001M03_obs'].attrs['long_name']='MIMOC climatological mean in march over 2003-2014 based on a density criteria of 0.01 kg/m^3'
                ds_out['MLDd001M03_obs'].attrs['units']='m'

                ds_out['MLDd001M09_obs']= (('yobs','xobs'), mld_obs_m09.values.astype('float32')) 
                ds_out['MLDd001M09_obs'].attrs['long_name']='MIMOC climatological mean in september over 2003-2014 based on a density criteria of 0.01 kg/m^3'
                ds_out['MLDd001M09_obs'].attrs['units']='m'

                ds_out['lat_obs']= (('yobs','xobs'), lat_obs.astype('float32')) 
                ds_out['lat_obs'].attrs['long_name']='Degrees north'
                ds_out['lat_obs'].attrs['units']='Deg'
                
                ds_out['lon_obs']= (('yobs','xobs'), lon_obs.astype('float32')) 
                ds_out['lon_obs'].attrs['long_name']='Degrees east'
                ds_out['lon_obs'].attrs['units']='Deg'
                
                ds_out['lat_mod']= (('y','x'), lat.values.astype('float32')) 
                ds_out['lat_mod'].attrs['long_name']='Degrees north'
                ds_out['lat_mod'].attrs['units']='Deg'
                
                ds_out['lon_mod']= (('y','x'), lon.values.astype('float32')) 
                ds_out['lon_mod'].attrs['long_name']='Degrees east'
                ds_out['lon_mod'].attrs['units']='Deg'
                
                ds_out = ds_out.set_coords(['lat_obs','lon_obs','lat_mod','lon_mod'])

                # Write the NetCDF file 
                ds_out.attrs['History'] = 'Diagnostics have been calculated using the Arctic monitoring tool '
                ds_out.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
                nc_f = './NETCDF/'+CONFIG+'-'+CASE+'_MLDClim_'+'y'+climyear+'.nc'
                ds_out.to_netcdf(nc_f,engine='netcdf4')

# To plot DYN variable
if DYN_maps : 
        num_fram=220
        # mean PSI
        zMyvar='sobarstf'   ; fram=num_fram+4
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, My_var1, zMyvar, climyear, zfram=fram )
        # Surface mean EKE
        zMyvar='voeke'   ; fram=num_fram+1
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, My_var1SeasM, zMyvar, climyear, slev=str(int(z[0])) , zfram=fram )
        # ~100m depth mean EKE
        zMyvar='voeke'   ; fram=num_fram+3
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, My_var1SeasS, zMyvar, climyear, slev=str(int(z[23])), zfram=fram )
        # EKE from observation
        zMyvar='voeke'   ; fram=num_fram+2
        obs_eke,lon_obs,lat_obs=CREG_maps_func.EKE_OBS(t_year=s_year)
        obs_eke=npy.ma.masked_where(obs_eke >= 9e20, obs_eke)
        CREG_maps_func.simple_maps( lon_obs, lat_obs, CONFIG, CASE, obs_eke, zMyvar, climyear, slev=str(int(z[0])), zfram=fram, plot_obs=1 )

        zfile_ext='_DYNClim_'
        plt.tight_layout()
        plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.png',dpi=300)

        if NCDF_OUT:
                ds_out = xr.Dataset()
                
                # DYN fields
                #######################
                ds_out['PSI_mod']= (('y','x'), My_var1.values.astype('float32')) 
                ds_out['PSI_mod'].attrs['long_name']='Model annual mean barotropic streamfunction '
                ds_out['PSI_mod'].attrs['units']='Sv'

                ds_out['EKESurf_mod']= (('y','x'), My_var1SeasM.values.astype('float32')) 
                ds_out['EKESurf_mod'].attrs['long_name']='Model annual mean EKE at the surface'
                ds_out['EKESurf_mod'].attrs['units']='cm2/s2'

                ds_out['EKEz100_mod']= (('y','x'), My_var1SeasS.values.astype('float32')) 
                ds_out['EKEz100_mod'].attrs['long_name']='Model annual mean EKE @ ~100m depth'
                ds_out['EKEz100_mod'].attrs['units']='cm2/s2'

                ds_out['EKESurf_obs']= (('yobs','xobs'), obs_eke.astype('float32')) 
                if s_year >= 2003 and s_year <= 2014 :
                        ds_out['EKESurf_obs'].attrs['long_name']='EKE annual mean derived from DOT field (Armitage et al. 2017) in '+str(s_year)
                else :
                        ds_out['EKESurf_obs'].attrs['long_name']='EKE climatological mean derived from DOT field (Armitage et al. 2017) over 2003-2014'
                ds_out['EKESurf_obs'].attrs['units']='cm2/s2'

                ds_out['lat_obs']= (('yobs','xobs'), lat_obs.astype('float32')) 
                ds_out['lat_obs'].attrs['long_name']='Degrees north'
                ds_out['lat_obs'].attrs['units']='Deg'
                
                ds_out['lon_obs']= (('yobs','xobs'), lon_obs.astype('float32')) 
                ds_out['lon_obs'].attrs['long_name']='Degrees east'
                ds_out['lon_obs'].attrs['units']='Deg'
                
                ds_out['lat_mod']= (('y','x'), lat.values.astype('float32')) 
                ds_out['lat_mod'].attrs['long_name']='Degrees north'
                ds_out['lat_mod'].attrs['units']='Deg'
                
                ds_out['lon_mod']= (('y','x'), lon.values.astype('float32')) 
                ds_out['lon_mod'].attrs['long_name']='Degrees east'
                ds_out['lon_mod'].attrs['units']='Deg'

                ds_out = ds_out.set_coords(['lat_obs','lon_obs','lat_mod','lon_mod'])

                # Write the NetCDF file 
                ds_out.attrs['History'] = 'Diagnostics have been calculated using the Arctic monitoring tool '
                ds_out.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
                nc_f = './NETCDF/'+CONFIG+'-'+CASE+'_DYNClim_'+'y'+climyear+'.nc'
                ds_out.to_netcdf(nc_f,engine='netcdf4')


# To plot T/S variables
if TSD_maps : 
        num_fram=220
        # Surface temperature
        zMyvar='votemper'   ; fram=num_fram+1
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, npy.squeeze(My_var1T[0,:,:]-My_varTinit[0,:,:]), zMyvar, climyear, slev=str(int(z[0])) , zfram=fram, ano=1 )
        # ~100m temperature
        zMyvar='votemper'   ; fram=num_fram+2
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, npy.squeeze(My_var1T[23,:,:]-My_varTinit[23,:,:]), zMyvar, climyear, slev=str(int(z[23])), zfram=fram, ano=1 )
        # Surface salinity
        zMyvar='vosaline'   ; fram=num_fram+3
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, npy.squeeze(My_var1S[0,:,:]-My_varSinit[0,:,:]), zMyvar, climyear, slev=str(int(z[0])), zfram=fram, ano=1 )
        # ~100m  salinity
        zMyvar='vosaline'   ; fram=num_fram+4
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, npy.squeeze(My_var1S[23,:,:]-My_varSinit[23,:,:]), zMyvar, climyear, slev=str(int(z[23])), zfram=fram, ano=1 )

        zfile_ext='_TSDIffClim_@'+str(int(z[0]))+'m@'+str(int(z[23]))+'m_'
        plt.tight_layout()
        plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.png',dpi=300)

        plt.clf()
        num_fram=220
        # ~200m temperature
        zMyvar='votemper'   ; fram=num_fram+1
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, npy.squeeze(My_var1T[30,:,:]-My_varTinit[30,:,:]), zMyvar, climyear, slev=str(int(z[30])) , zfram=fram, ano=1 )
        # ~300m temperature
        zMyvar='votemper'   ; fram=num_fram+2
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, npy.squeeze(My_var1T[34,:,:]-My_varTinit[34,:,:]), zMyvar, climyear, slev=str(int(z[34])), zfram=fram, ano=1 )
        # ~200m salinity
        zMyvar='vosaline'   ; fram=num_fram+3
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, npy.squeeze(My_var1S[30,:,:]-My_varSinit[30,:,:]), zMyvar, climyear, slev=str(int(z[30])), zfram=fram, ano=1 )
        # ~300m  salinity
        zMyvar='vosaline'   ; fram=num_fram+4
        CREG_maps_func.simple_maps( lon, lat, CONFIG, CASE, npy.squeeze(My_var1S[34,:,:]-My_varSinit[34,:,:]), zMyvar, climyear, slev=str(int(z[34])), zfram=fram, ano=1 )

        zfile_ext='_TSDIffClim_@'+str(int(z[30]))+'m@'+str(int(z[34]))+'m_'
        plt.tight_layout()
        plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.png',dpi=300)



# To plot ATL variables
if ATL_maps : 
        # MLD IN THE LABRADOR SEA IN MARCH
        ###################################
        plt.figure()
        num_fram=110
        # March mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+1
        my_cblab=r'MLD (m)'   ;   my_cmap=plt.get_cmap('Blues')
        ztitle=CASE +' mean MLD01 over \n'+climyear+'  m03'
        vmin=0. ; vmax=2400. ; vint=100.   ;   contours=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]
        limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]

        plt.subplot(fram)
        zoutmap=CREG_maps_func.Atl_Bat(ztype='isol1000',zarea='labsea')
        My_var1SeasM.shape
        CREG_maps_func.Atl_plot(lon, lat, My_var1SeasM, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='labsea')

        zfile_ext='_LAB_MLDClimm03_'
        plt.tight_layout()
        plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.png',dpi=300)

        # Plot the Time-serie for MLD at a specific location K1 mooring in the Labrador Sea and in Irminger Sea
        # After Schott et al. DSRI2009 56.33N, -52.40W
        plt.clf()
        plt.figure()
        i_K1=173   ;   j_K1=168   # CREG025.L75 C-type indices
        #i_K1=518   ;   j_K1=502   # CREG12.L75 C-type indices
        ax=plt.subplot(211)
        # In Lab. Sea
        plt.plot(time_axis,-1.*npy.squeeze(Mdata_read[:,j_K1:j_K1+1,i_K1:i_K1+1]),linewidth=0.8, color='k', label='Lab Sea K1')
        # Plot obs. MLD in March
        year_obs=npy.arange(1995,2006,1)+0.20547945   ; mld_obs=[-2300.,-1300.,-1400.,-1000.,-1000.,-1100.,-1100.,-1200.,-1400.,-700.,-1300.]
        plt.scatter(year_obs,mld_obs)

        # In Irm. Sea
        i_K1=232   ;   j_K1=192   # CREG025.L75 C-type indices geo loc   60.88N  -36.99W
        #i_K1=697   ;   j_K1=577   # CREG12.L75 C-type indices geo loc   60.88N  -36.99W
        plt.plot(time_axis,-1.*npy.squeeze(Mdata_read[:,j_K1:j_K1+1,i_K1:i_K1+1]),linewidth=0.8, color='g', label='Irm Sea ')
        plt.title(CASE+' MLD 0.01 in Lab. & Irm. Seas \n '+str(s_year)+str(e_year),size=9)
        plt.ylabel('Mean depth \n'+r'(m)', size=7)
        plt.ylim([-2500.,0.])
        plt.xticks(newlocsx,newlabelsx,size=5)
        plt.setp(ax.get_xticklabels(),rotation=90)
        plt.yticks(size=6)
        plt.grid(True, linewidth=0.7,linestyle='--',alpha=0.7,color='grey')
        plt.legend(loc='lower center',ncol=2)
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=5.)
        
        zfile_ext='_LAB-IRM_MLDClim_LGTS_'
        plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.png',dpi=300)

        # Add an artificial mooring within the deepest convection area 
        # -54W 58N
        #  dl_dis=    1.634 km
        #      507       507       541       541
        # -54.0272  -54.0272   57.9970   57.9970
        plt.clf()
        plt.figure()
        i_K1=169   ;   j_K1=181   # CREG025.L75 C-type indices
        #i_K1=506   ;   j_K1=540   # CREG12.L75 C-type indices
        ax=plt.subplot(211)
        # In Lab. Sea
        plt.plot(time_axis,-1.*npy.squeeze(Mdata_read[:,j_K1:j_K1+1,i_K1:i_K1+1]),linewidth=0.8, color='k', label='Lab Sea DeepConv')
        # Plot obs. MLD in March
        year_obs=npy.arange(1995,2006,1)+0.20547945   ; mld_obs=[-2300.,-1300.,-1400.,-1000.,-1000.,-1100.,-1100.,-1200.,-1400.,-700.,-1300.]
        plt.scatter(year_obs,mld_obs)
        plt.title(CASE+' MLD 0.01 in Lab. @ -54W,58N \n '+str(s_year)+str(e_year),size=9)
        plt.ylabel('Mean depth \n'+r'(m)', size=7)
        plt.ylim([-3500.,0.])
        plt.xticks(newlocsx,newlabelsx,size=5)
        plt.setp(ax.get_xticklabels(),rotation=90)
        plt.yticks(size=6)
        plt.grid(True, linewidth=0.7,linestyle='--',alpha=0.7,color='grey')
        plt.legend(loc='lower center',ncol=2)
        leg = plt.gca().get_legend()
        ltext = leg.get_texts()
        plt.setp(ltext, fontsize=5.)

        zfile_ext='_LABM52W58N-MLDClim_LGTS_'
        plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.png',dpi=300)


        plt.clf()
        plt.figure()
        # MLD IN THE IRMINGER SEA IN MARCH
        ###################################
        num_fram=110
        # March mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+1
        my_cblab=r'MLD (m)'   ;   my_cmap=plt.get_cmap('Blues')
        ztitle=CASE +' mean MLD01 over \n'+climyear+'  m03'
        vmin=0. ; vmax=1600. ; vint=100.   ;   contours=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
        limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
        plt.subplot(fram)
        zoutmap=CREG_maps_func.Atl_Bat(ztype='isol1000',zarea='irmsea')
        CREG_maps_func.Atl_plot(lon, lat, My_var1SeasM, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='irmsea')

        zfile_ext='_IRM_MLDClimm03_'
        plt.tight_layout()
        plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.png',dpi=300)

        plt.clf()
        plt.figure()
        # MLD IN THE GIN SEAS IN MARCH
        ###################################
        num_fram=110
        # March mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+1
        my_cblab=r'MLD (m)'   ;   my_cmap=plt.get_cmap('Blues')
        ztitle=CASE +' mean MLD01 over \n'+climyear+'  m03'
        vmin=0. ; vmax=1600. ; vint=100.   ;   contours=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
        limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
        plt.subplot(fram)
        zoutmap=CREG_maps_func.Atl_Bat(ztype='isol1000',zarea='ginsea')
        CREG_maps_func.Atl_plot(lon, lat, My_var1SeasM, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='ginsea')

        zfile_ext='_GIN_MLDClimm03_'
        plt.tight_layout()
        plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.png',dpi=300)

        # PLOT ISOTHERM 17 Deg OFF CAPE HATTERAS
        #################################
        num_fram=110
        zMyvar='votemper'   ; fram=num_fram+1
        my_cblab=r'ISO 17 (DegC)'   ;   my_cmap=plt.get_cmap('jet')
        ztitle=CASE +' mean Iso 17 DegC over \n'+climyear
        vmin=0. ; vmax=2400. ; vint=100.   ;   contours=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]
        limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]

        plt.subplot(fram)
        klev=29
        zoutmap=CREG_maps_func.Atl_Bat(ztype='isol1000', zarea='GulfS')
        lon = lon.values   ; lat = lat.values 
        CREG_maps_func.Atl_plot(lon, lat, npy.squeeze(My_var1T[klev,:,:])   , contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='GulfS')
        CREG_maps_func.Atl_plot(lon, lat, npy.squeeze(My_varTinit[klev,:,:]), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='GulfS', data_ref=True)

        zfile_ext='_ISO17Clim_'
        plt.tight_layout()
        plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.png',dpi=300)

        # PLOT SSH OVER THE ATLANTIC AREA
        #################################
        num_fram=110
        zMyvar='sossheig'   ; fram=num_fram+1
        my_cblab=r'SSH (cm)'   ;   my_cmap=plt.get_cmap('coolwarm')
        ztitle=CASE +' mean SSH over \n'+climyear
        vmin=-100. ; vmax=100. ; vint=5.  ;   contours=npy.arange(vmin,vmax+vint,vint)
        limits=[vmin,vmax,vint]           ;   myticks=npy.arange(vmin,vmax+vint,vint)

        plt.figure()
        #plt.subplot(fram)
        zoutmap=CREG_maps_func.Atl_Bat(ztype='isol1000')
        CREG_maps_func.Atl_plot(lon, lat, My_var1ssh*100.  , contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar)
        #CREG_maps_func.Atl_plot(lon, lat, npy.squeeze(My_varTinit[klev,:,:]), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='GulfS', data_ref=True)

        zfile_ext='_SSHClim_'
        #plt.tight_layout()
        plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.png',dpi=300)


# To plot AMOC 
if MOC_maps : 
        plt_AMOCTS=True

        if  plt_AMOCTS: 
        #if not plt_AMOCTS: 
                # AAMOC 
                #######
                plt.figure()
                num_fram=210
                fram=num_fram+1
                my_cblab=r'AMOC (Sv)'   ;   my_cmap=plt.get_cmap('jet')
                ztitle=CASE +' mean AAMOC over \n'+climyear
                vmin=-15. ; vmax=15. ; vint=1.    ;   contours=npy.arange(vmin,vmax+vint,vint)
                limits=[vmin,vmax,vint]  ;             myticks=npy.arange(vmin,vmax+vint,vint)
                norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
                plt.subplot(fram)

                lat2Dz=npy.reshape(lat,(lat.size,1)).T
                ypltz = npy.repeat(lat2Dz,z.shape[0],axis=0)
                locpath='./'
                locfile='Bathymetry.nc'
                if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
                        ds_msk = xr.open_dataset(locpath+locfile)[['nav_lat']]
                        full_lat = ds_msk['nav_lat']
                if CONFIG == 'CREG12.L75' : 
                	select_ylat=full_lat[:,862]
                else: 
                	select_ylat=full_lat[:,288]
                select_ylat_reshape=npy.reshape(select_ylat.values,(select_ylat.size,1))
                ypltz=npy.repeat(select_ylat_reshape,z.shape[0],axis=1).T

                z2dt=npy.reshape(z.values,(z.size,1))
                zplt = npy.repeat(z2dt,lat.shape[0],axis=1)

                plt.title(ztitle,fontsize=6)
                plt.contourf(ypltz,zplt*(-1.e-3),npy.squeeze(My_MOC),contours,cmap=my_cmap,norm=norm,extend='both')
                contours=npy.arange(vmin,vmax+vint,2.*vint)
                C=plt.contour(ypltz,zplt*(-1.e-3),npy.squeeze(My_MOC),linewidths=0.5,levels=contours, colors='k')
                plt.clabel(C, C.levels, inline=True, fmt='%3.0f', fontsize=6)

                zfile_ext='_AAMOCClim_'
                #plt.tight_layout()
                plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.png',dpi=300)


        if plt_AMOCTS :
                # Max. AAMOC Time-series @ 40N or 43N
                #######################################

                time_grid=npy.arange(s_year,e_year+1,1.,dtype=int)
                newlocsx  = npy.array(time_grid,'f')
                newlabelsx = npy.array(time_grid,'i')

                jloc= 205 # Lat 40N
                #jloc= 252 # Lat 43N
                if jloc == 205   : zlat='L40N'
                elif jloc == 252 : zlat='L43N'

                plt.clf()
                ax=plt.subplot(211)
                plt.plot( time_grid, My_MOC.isel(y=jloc).max(dim='depthw'), 'g',linewidth=1., label=CONFIG)
                #plt.plot( time_grid, npy.max(My_MOC[:,:,jloc],axis=1)   , 'g',linewidth=1., label=CONFIG)
                plt.title(CASE+' Max AAMOCz @ '+zlat+'\n '+str(s_year)+str(e_year),size=9)
                plt.ylabel('Max AAMOCz  \n'+r'(Sv)', size=7)
                plt.ylim([10.,15.])
                plt.xticks(newlocsx,newlabelsx,size=5)
                plt.setp(ax.get_xticklabels(),rotation=90)
                plt.yticks(size=6)
                plt.grid(True, linewidth=0.7,linestyle='--',alpha=0.7,color='grey')

                ax=plt.subplot(212)
                plt.plot( time_grid, My_MOC.max().item(), 'k-*',linewidth=1., label=CONFIG)
                #plt.plot( time_grid, npy.max(My_MOC[:,:,:],axis=(1,2))   , 'k',linewidth=1., label=CONFIG)
                #plt.title(CASE+' Max AAMOCz '+'\n '+str(s_year)+str(e_year),size=9)
                plt.ylabel('Global Max AAMOCz  \n'+r'(Sv)', size=7)
                plt.ylim([10.,20.])
                plt.xticks(newlocsx,newlabelsx,size=5)
                plt.setp(ax.get_xticklabels(),rotation=90)
                plt.yticks(size=6)
                plt.grid(True, linewidth=0.7,linestyle='--',alpha=0.7,color='grey')


                zfile_ext='_MaxAAMOCz_'+zlat+'_TiSe_'
                #plt.tight_layout()
                plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.png',dpi=300)
