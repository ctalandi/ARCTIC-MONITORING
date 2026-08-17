"""
CREG_maps_func.py

Description:
This module defines a set of functions dedicated to :
- ATL_MAPSF   : plot Atlantic variables such as MLD, SSH in differents areas such as GIN, LAB, IRM seas
- AWT_MAPSF   : plot the Atlantic Water maximum temperature as the associated depth
- BFG_MAPSF   : plot the Beaufort Gyre center based on SSH
- BFG_COMPUTE : compute metrics such as SSH closed contours and mass center location
- CONV_SA2PS  : convert Absolute Salinity to Practical Salinity units using the GSW package 
- CREG_MSK    : read the CREG configuration mask 
- CREG_INIT   : read CREG initial state temperature/salinity 
- DO_MAPS     : intermediate function before calling the projection itself 
- DYN_MAPSF   : plot the barotropic stream function and the EKE at surface ~69m and 508m model depths
- EKE_CALC    : compute the EKE using monthly mean and annual mean velocities 
- FND_CEDGE   : algorithm to compute the closed contours using SSH 
- FWC_MAPSF   : plot the SSH and FWC (based on a salinity ref of 34.8 PSU)
- ICE_MAPSF   : plot ICE concentration & thickness
- BATHY_MAP   : plot few iso-bathymetric contours 
- MLD_CALC    : compture de the MLD using a 0.1 kg m-3 density criteria
- MLD_MAPSF   : plot the Mixed Layer Depth 
- MOC_MAPSF   : plot the AMOC and its maximum time series
- MTS_MAPSF   : plot the Mixed Layer mean T/S 
- PROJ_PLOT   : apply the right projection over the Arctic for the plot 
- P3Dzyx      : extend the pressure vector to spatial 3D 
- TOPOS_CALC  : compute the annual mean topostrophy 
- TSD_MAPSF   : plot the T/S drift at the surface, ~100m, ~200m & ~300m

Author:
Claude Talandier (claude.talandier@cnrs.fr)
"""

import matplotlib
matplotlib.use('Agg')
import sys
import numpy as npy
from CREG_maps_cont import *
from CREG_maps_obse import *
from checkfile import *
import subprocess
import xarray as xr 
from datetime import datetime
from cartopy import crs as ccrs
import cartopy
import matplotlib.path as mpath
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import pandas as pd
import calendar
import time
import gsw as gsw 
from pathlib import Path 
from xnemogcm import open_domain_cfg, open_nemo, process_nemo, open_namelist, open_nemo_and_domain_cfg
from xnemogcm import __version__ as xnemogcm_version
import xgcm
from fsspec.implementations.local import LocalFileSystem
fs = LocalFileSystem()


# Matplotlib
try:
	import matplotlib.pylab as plt
	import matplotlib as mpl
	from matplotlib import rcParams
except:
	print('  matplotlib is not available on your machine')
	print('  check python path or install this package') ; exit()

# Basemap
try:
	from mpl_toolkits.basemap import Basemap
except:
	print('  Basemap is not available on your machine')
	print('  check python path or install this package') ; exit()

################################################################################################################################
def BFG_MAPSF( zlon, zlat, zvar_ssh, zbathy, zarea, zCONF, zCASE, zs_year, zlgTS_ys, zlgTS_ye, zncout ) :
################################################################################################################################
	"""
	Function dedicated to compute the Beaufort Gyre mass center as closed contours using the SSH 
	
	History: This code has been developed by Heather Regan and slightly adapted to be included into the MONARC
		 See Regan et al. JPO2020 ; https://doi.org/10.1175/JPO-D-19-0234.1
	
	Input:
	    zlon     : longitude 2D 
	    zlat     : latitude 2D 
	    zvar_ssh : monbthly SSH over one year 
	    zbathy   : ocean floor depth in meters 
	    zarea    : surface of each grid cell 
	    zCONF    : configuration name 
	    zCASE    : experiment name associated to the configuration
	    zs_year  : current year 
	    zlgTS_ys : first year for long time-series 
	    zlgTS_ye : last year for long time-series 
	    zncout   : logical to perform (or not) outputs into a netcdf file 
	
	Output:
	    LongTS_OBS_icevol          : PIOMAS ice thickness time-series 
	    LongTS_OBS_iceext          : NSIDC v6 ice extent time-series 
	    LongTS_OBS_Septiceext      : NSIDC v6 september ice extent time-series 
	    LongTS_OBS_iceare          : NSIDC v6 ice area time-series
	    IABPObservations           : IABP ice drift time-series 
	    time_axis_IABP             : IABP data time axis 
	    time_axis_PIO              : PIOMAS data time axis 
	    time_axis_NSIDC            : NSIDC data time axis 
	    LongTS_OBS_Septiceext_time : NSIDC data time axis for September minimum ice extent 
	"""
	# ----------------------------------------------------------------------

	## This is the increment that will be iterated on for identifying closed contours
	## smaller value = smaller increments, so takes longer, but less likely to miss small variations in contour edges so better for higher resolutions
	increment_in = 0.05; incstr = '5cm' #TOEDIT

	######################### IMPORTANT NOTES ######################
	# This script assumes depth values are greater than zero. Be sure to multiply the depth array be -1 if this is the opposite!
	# This script assumes -180 < longitude < 180. If it goes from 0-360, be sure to adjust for this!
	# This script needs four input variables: SSH, latitude, longitude, and depth. All are required in netcdf format
	# Parts that need to be edited or verified by the user can be found by #TOEDIT
	######################### IMPORTANT NOTES ######################

	#------------------------------------------------------------------------------------------------------------------------
	########################################
	# Compute the Beaufort Gyre center location and closed contour
	########################################
	#------------------------------------------------------------------------------------------------------------------------

	timetype = 'yearly'

	# Compute the annual mean SSH
	zssh_ym = zvar_ssh.mean(dim='time_counter').values.squeeze()
	zincr = increment_in*1

	start_time = time.time()
	msk_ym, BG_maxval_ym, BG_minval_ym, BG_maxlat_ym, BG_maxlon_ym, BG_area_ym = BFG_COMPUTE( zlon, zlat, zssh_ym, zbathy, 'SSH', zincr, zarea )
	print('			Computing wall time (s) : ',time.time() - start_time)

	#------------------------------------------------------------------------------------------------------------------------

	msk = xr.full_like(zvar_ssh, fill_value=0.)
	BG_maxval =  xr.DataArray(npy.zeros(12), dims=['time'])
	BG_minval = xr.DataArray(npy.zeros(12), dims=['time'])
	BG_maxlat = xr.DataArray(npy.zeros(12), dims=['time'])
	BG_maxlon = xr.DataArray(npy.zeros(12), dims=['time'])
	BG_area = xr.DataArray(npy.zeros(12), dims=['time'])

	timetype = 'monthly'

	# For monthly mean 
	for zmm in range(0,12):
	
		print()
		print(' Considered month : ',calendar.month_name[zmm+1])

		zincr = increment_in*1

		# Select one month 
		zssh = zvar_ssh.isel(time_counter=zmm).values.squeeze()

		start_time = time.time()
		msk[zmm,:,:], BG_maxval[zmm], BG_minval[zmm], BG_maxlat[zmm], BG_maxlon[zmm], BG_area[zmm] = BFG_COMPUTE( zlon, zlat, zssh, zbathy, 'SSH', zincr, zarea )
		print('			Computing wall time (s) : ',time.time() - start_time)

	#------------------------------------------------------------------------------------------------------------------------

	# Plot the yearly closed contours as the BFG center as well 
	plt.clf()
	fig=plt.figure()
	projection = ccrs.Orthographic(central_longitude=-160, central_latitude=60)
	
	fram=211
	ax = fig.add_subplot(fram, projection=projection)
	zoutmap = BATHY_MAP( ztype='isol1000', zarea='cassis_BGZoom', ax=ax )	

	if npy.nansum(msk_ym) > 0:
		msk_plot = xr.where( npy.isnan(msk_ym), 0., msk_ym*1 )
		CS2 = ax.contour( zlon, zlat, msk_plot, linewidths=0.5, colors='k', transform=ccrs.PlateCarree() )
	# Get indices of the BFG center 
	[r,c] = npy.nonzero( msk_ym*zssh_ym == npy.nanmax(msk_ym*zssh_ym) )
	
	clat = [zlat[r.item(),c.item()].values,]
	clon = [zlon[r.item(),c.item()].values,]
	ax.scatter(clon[0],clat[0], s=10, marker='o', color='k', transform=ccrs.PlateCarree())
	plt.title( zCASE+' BFG SSH contours \n yearly mean SSH '+str(zs_year), fontsize=6 )


        # Plot the monthly closed contours as the BFG center as well 
	cmap = plt.get_cmap('Spectral_r')
	colors = [cmap(i) for i in npy.linspace(0, 1, 12)]

	fram=212
	ax = fig.add_subplot(fram, projection=projection)
	zoutmap = BATHY_MAP( ztype='isol1000', zarea='cassis_BGZoom', ax=ax )	

	for zmm in range(0,12):
		if npy.nansum(msk[zmm,:,:]) > 0:
			msk_plot = xr.where( npy.isnan(msk[zmm,:,:]), 0., msk[zmm,:,:]*1 )
			CS2 = ax.contour( zlon, zlat, msk_plot, linewidths=0.5, colors=colors[zmm], transform=ccrs.PlateCarree() )
		# Get indices of the BFG center 
		[r,c] = npy.nonzero( msk[zmm,:,:]*zvar_ssh.isel(time_counter=zmm) == npy.nanmax(msk[zmm,:,:]*zvar_ssh.isel(time_counter=zmm)) )
		
		clat = [zlat[r.values.item(),c.values.item()],]
		clon = [zlon[r.values.item(),c.values.item()],]
		ax.scatter(clon[0],clat[0], s=10, marker='o', color=colors[zmm], transform=ccrs.PlateCarree())

	legend_elements = [ Line2D([0], [0], color=colors[i], lw=1, label=calendar.month_name[i+1]) for i in range(len(colors)) ]
	plt.legend(handles=legend_elements, bbox_to_anchor=(0.95, 0.85), ncol=1, fontsize=6, handlelength=1.5, handletextpad=0.5, borderpad=0.2, labelspacing=0.2)
	plt.title( ' monthly mean SSH '+str(zs_year), fontsize=6 )
	plt.tight_layout()

	zfile_ext='_BFGCenter_'+'y'+str(zs_year)
	plt.savefig(zCONF+'-'+zCASE+zfile_ext+'.png',dpi=300)

	plt.close()

	#------------------------------------------------------------------------------------------------------------------------

	if zncout:
		ds_out = xr.Dataset()
		
		ds_out.coords['BGlat'] = (('y','x'), zlat.values.astype('float32')) 
		ds_out.coords['BGlat'].attrs['long_name'] = 'latitude'
		ds_out.coords['BGlat'].attrs['units'] = 'degrees_north'
		
		ds_out.coords['BGlon'] = (('y','x'), zlon.values.astype('float32')) 
		ds_out.coords['BGlon'].attrs['long_name'] = 'longitude'
		ds_out.coords['BGlon'].attrs['units'] = 'degrees_east'
 
		# Save diags. based on yearly mean SSH
		timevalue_ym = pd.to_datetime(str(zs_year)+'-06-30')
		ds_out.coords['time_ym'] = (('time_ym'), [timevalue_ym] )
		
		ds_out['BGmask_ym'] = (('time_ym','y','x'), [msk_ym.astype('float32')]) 
		ds_out['BGmask_ym'].attrs['long_name'] = 'Beaufort_Gyre_mask yearly mean'
		ds_out['BGmask_ym'].attrs['units'] = '-'
		
		ds_out['BGmax_ym'] = (('time_ym'), [BG_maxval_ym.astype('float32')]) 
		ds_out['BGmax_ym'].attrs['long_name'] = 'Maximum_ssh_in_gyre yearly mean'
		ds_out['BGmax_ym'].attrs['units'] = 'metres'
		
		ds_out['BGmin_ym'] = (('time_ym'), [BG_minval_ym.astype('float32')]) 
		ds_out['BGmin_ym'].attrs['long_name'] = 'Minimum_ssh_in_gyre yearly mean'
		ds_out['BGmin_ym'].attrs['units'] = 'metres'
		
		ds_out['BGmaxhlat_ym'] = (('time_ym'), [BG_maxlat_ym.astype('float32')]) 
		ds_out['BGmaxhlat_ym'].attrs['long_name'] = 'latitude_of_max_ssh_in_gyre yearly mean'
		ds_out['BGmaxhlat_ym'].attrs['units'] = 'degrees'
		
		ds_out['BGmaxhlon_ym'] = (('time_ym'), [BG_maxlon_ym.astype('float32')]) 
		ds_out['BGmaxhlon_ym'].attrs['long_name'] = 'longitude_of_max_ssh_in_gyre yearly mean'
		ds_out['BGmaxhlon_ym'].attrs['units'] = 'degrees'
		
		ds_out['BGarea_ym'] = (('time_ym'), [BG_area_ym.astype('float32')]) 
		ds_out['BGarea_ym'].attrs['long_name'] = 'area_of_gyre yearly mean'
		ds_out['BGarea_ym'].attrs['units'] = 'metres_squared'
		
		# Save diags. based on monthly mean SSH
		timevalue = pd.date_range(start=str(zs_year)+'-01',end=str(zs_year)+'-12',freq='MS')+ pd.DateOffset(days=14)
		ds_out.coords['time'] = (('time'), timevalue )
		
		ds_out['BGmask'] = (('time','y','x'), msk.values.astype('float32')) 
		ds_out['BGmask'].attrs['long_name'] = 'Beaufort_Gyre_mask'
		ds_out['BGmask'].attrs['units'] = '-'
		
		ds_out['BGmax'] = (('time'), BG_maxval.values.astype('float32')) 
		ds_out['BGmax'].attrs['long_name'] = 'Maximum_ssh_in_gyre'
		ds_out['BGmax'].attrs['units'] = 'metres'
		
		ds_out['BGmin'] = (('time'), BG_minval.values.astype('float32')) 
		ds_out['BGmin'].attrs['long_name'] = 'Minimum_ssh_in_gyre'
		ds_out['BGmin'].attrs['units'] = 'metres'
		
		ds_out['BGmaxhlat'] = (('time'), BG_maxlat.values.astype('float32')) 
		ds_out['BGmaxhlat'].attrs['long_name'] = 'latitude_of_max_ssh_in_gyre'
		ds_out['BGmaxhlat'].attrs['units'] = 'degrees'
		
		ds_out['BGmaxhlon'] = (('time'), BG_maxlon.values.astype('float32')) 
		ds_out['BGmaxhlon'].attrs['long_name'] = 'longitude_of_max_ssh_in_gyre'
		ds_out['BGmaxhlon'].attrs['units'] = 'degrees'
		
		ds_out['BGarea'] = (('time'), BG_area.values.astype('float32')) 
		ds_out['BGarea'].attrs['long_name'] = 'area_of_gyre'
		ds_out['BGarea'].attrs['units'] = 'metres_squared'

		# Write the NetCDF file 
		ds_out.attrs['History'] = 'Diagnostics have been calculated using the Arctic monitoring tool '
		ds_out.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
		nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_BFGfromSSH_inc'+incstr+'_y'+str(zs_year)+'.nc'
		ds_out.to_netcdf(nc_f,engine='netcdf4')
	
	#------------------------------------------------------------------------------------------------------------------------
	########################################
	# Plot LONG TIME-SERIES
	########################################
	#------------------------------------------------------------------------------------------------------------------------
	
	if zlgTS_ye-zlgTS_ys+1 > 1 :

		print()
		print('				##################################################################  ')
		print('				##################################################################  ')
		print('				######## PLOT BFG CENTER CARATERISTICS LONG TIME-SERIES ##########  ')
		print('				##################################################################  ')
		print('				##################################################################  ')
		print()
		
		# Read Obs. data set 
		#######################################################################################################
		obs_datafile = './DATA/BGmask_2003to2014.nc'
		ds_obsrssh = xr.open_dataset( obs_datafile, engine="netcdf4" ) 
		# Rebuild a proper new time axis 
		timevalue = pd.date_range(start='2003-01',end='2014-12',freq='MS')+ pd.DateOffset(days=14)
		obs_time_axis = timevalue.year + (timevalue.month - 1 + 0.5) / 12
		ds_obs = xr.Dataset()
		ds_obs.coords['time'] = (('time'), timevalue)
		ds_obs['maxheight'] = (('time'), ds_obsrssh["maxheight"].values.squeeze())
		ds_obs['area_m2'] = (('time'), ds_obsrssh["area_m2"].values.squeeze())
		ds_obs['maxh_lat'] = (('time'), ds_obsrssh["maxh_lat"].values.squeeze())
		ds_obs['maxh_lon'] = (('time'), ds_obsrssh["maxh_lon"].values.squeeze())
		
		# Read model data set 
		#######################################################################################################
		locpath = './NETCDF/'
		locfile = zCONF+'-'+zCASE+'_BFGfromSSH_inc'+incstr+'_y????.nc'
		ds_rssh = xr.open_mfdataset(locpath+locfile, engine='netcdf4', concat_dim=['time'], combine='nested', parallel=True)
		mod_time_axis = ds_rssh.time.dt.year.values + (ds_rssh.time.dt.month.values - 1 + 0.5) / 12
	
		# Make plots 
		#######################################################################################################
		plt.clf()
		xwind=410

		time_grid = npy.arange(1979,2025.,1.,dtype=int)
		newlocsx = npy.array(time_grid,'f')
		newlabelsx = npy.array(time_grid,'i')
		
		# Max SSH 
		################
		ax=plt.subplot(xwind+1)
		ax.set_title(zCASE,size=7)
		plt.plot(mod_time_axis,ds_rssh['BGmax'].values,color='k', linewidth=0.6, label='model')
		plt.plot(obs_time_axis,ds_obs['maxheight'],color='g', linewidth=0.6, label='obs')
		plt.xlim([1978,2025])
		plt.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
		plt.ylabel(' max ssh at gyre centre \n (m)',size=6)
		plt.xticks(newlocsx,newlabelsx,size=5)
		plt.setp(ax.get_xticklabels(),visible=False)
		plt.yticks(size=6)
		plt.legend(fontsize='small', ncol=2)
		leg = plt.gca().get_legend()
		ltext = leg.get_texts()
		plt.setp(ltext, fontsize=4)

		# Gyre area 
		################
		ax=plt.subplot(xwind+2)
		plt.plot(mod_time_axis,(ds_rssh['BGarea'].values)*1e-12,color='k', linewidth=0.6, label='model')
		plt.plot(obs_time_axis,(ds_obs['area_m2'].values)*1e-12,color='g', linewidth=0.6, label='obs')
		plt.xlim([1978,2025])
		plt.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
		ax.yaxis.set_major_locator(MultipleLocator(0.25))
		plt.ylabel('gyre area \n'+r'(x$10^6$ $km^{2}$)',size=6)
		plt.xticks(newlocsx,newlabelsx,size=5)
		plt.setp(ax.get_xticklabels(),visible=False)
		plt.yticks(size=6)
	
		# Latitude of Max SSH 
		######################
		ax=plt.subplot(xwind+3)
		plt.plot(mod_time_axis,ds_rssh["BGmaxhlat"].values,color='k', linewidth=0.6, label='model')
		plt.plot(obs_time_axis,ds_obs['maxh_lat'].values,color='g', linewidth=0.6, label='obs')
		plt.xlim([1978,2025])
		plt.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
		ax.yaxis.set_major_locator(MultipleLocator(2))
		plt.ylabel('latitude of max ssh',size=6)
		plt.xticks(newlocsx,newlabelsx,size=5)
		plt.setp(ax.get_xticklabels(),visible=False)
		plt.yticks(size=6)
	
		# Longitude of Max SSH 
		######################
		ax=plt.subplot(xwind+4)
		plt.plot(mod_time_axis,ds_rssh["BGmaxhlon"].values,color='k', linewidth=0.6, label='model')
		plt.plot(obs_time_axis,ds_obs['maxh_lon'].values,color='g', linewidth=0.6, label='obs')
		plt.xlim([1978,2025])
		plt.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
		ax.yaxis.set_major_locator(MultipleLocator(5))
		plt.xlabel('years' ,size=5)
		plt.ylabel('longitude of max ssh',size=6)
		plt.xticks(newlocsx,newlabelsx,size=5)
		plt.setp(ax.get_xticklabels(),rotation=90, fontsize=5)
		plt.yticks(size=6)
	
		plt.tight_layout()

		zfile_ext='_BFG_metrics_LGTS_y'+str(zlgTS_ys)+'LASTy'
		plt.savefig(zCONF+'-'+zCASE+zfile_ext+'.png',dpi=300)

	return

################################################################################################################################
def ICE_MAPSF( zlon, zlat, zIceThic, zIceConc_M, zIceConc_S, zCONF, zCASE, zc_year, zncout ) :
################################################################################################################################
        """
        Function dedicated to plot ice thickness and March/September concentration over the Arctic 

        Input:
            zlon       : longitude 2D 
            zlat       : latitude 2D 
            zIceThic   : ice thickness 
            zIceConc_M : March ice concentration 
            zIceConc_S : September ice concentration
            zCONF      : the configuration name 
            zCASE      : the experiment name associated to the configuration
            zc_year    : the current year 
            zncout     : logical to outputs (or not) results into a netcdf file 
        
        Output:
            None 
        """
        # ----------------------------------------------------------------------

        plt.clf()
        fig, axes = plt.subplots( 3, 2, figsize=(8.3,11.7),  subplot_kw={'projection': ccrs.NorthPolarStereo(central_longitude=-60, true_scale_latitude=65)} )
        axes = axes.flatten()

        # Annual mean Ice thickness
        zMyvar='sivolu'
        zIceThic = xr.where( zIceThic == 0., npy.nan, zIceThic )
        DO_MAPS( zlon, zlat, zCONF, zCASE, zIceThic, zMyvar, zc_year, ax=axes[0] )
        # March mean Ice fraction 
        zMyvar='siconc'
        DO_MAPS( zlon, zlat, zCONF, zCASE, zIceConc_M, zMyvar, zc_year, seas='m03', ax=axes[2] )
        # September mean Ice fraction 
        zMyvar='siconc'
        DO_MAPS( zlon, zlat, zCONF, zCASE, zIceConc_S, zMyvar, zc_year, seas='m09', ax=axes[4] )

        # Annual mean Ice thickness PIOMASS observations
        # WARNING this part has been interpolated on CREG025 grid directly 
        # Read lat, lon from CREG025 grid since the PIOMASS data is not available on CREG12.L75 grid
        zmask, plon, plat = CREG_MSK( zCONF, zCASE )

        obs_thick = ICE_THICK_OBS( zconfig=zCONF, t_year=zc_year )
        obs_thick = xr.where( npy.squeeze(zmask[0,:,:]) < 1., npy.nan, obs_thick )
        obs_thick = xr.where( obs_thick == 0., npy.nan, obs_thick )

        zMyvar='sivolu'
        DO_MAPS( plon, plat, zCONF, zCASE, obs_thick, zMyvar, zc_year, plot_obs=0, ax=axes[1] )

        # Read NSIDC obs. data 
        obs_conc_m03, obs_conc_m09, obs_lon, obs_lat = ICE_CONCE_OBS( t_year=zc_year )

        # March mean Ice fraction 
        zMyvar='siconc'
        DO_MAPS( obs_lon, obs_lat, zCONF, zCASE, obs_conc_m03, zMyvar, zc_year, seas='m03', plot_obs=1, ax=axes[3] )
        # September mean Ice fraction 
        zMyvar='siconc'
        DO_MAPS( obs_lon, obs_lat, zCONF, zCASE, obs_conc_m09, zMyvar, zc_year, seas='m09', plot_obs=1, ax=axes[5] )
        plt.tight_layout()

        zfile_ext='_ICEClim_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zc_year)+'.png',dpi=300)

        if zncout:
                ds_out = xr.Dataset()
                
                # ICE fields
                #######################
                ds_out['IceThick_mod']= (('y','x'), zIceThic.values.astype('float32')) 
                ds_out['IceThick_mod'].attrs['long_name']='Model annual mean Ice thickness'
                ds_out['IceThick_mod'].attrs['units']='m'
                
                ds_out['IceConceM03_mod']= (('y','x'), zIceConc_M.values.astype('float32')) 
                ds_out['IceConceM03_mod'].attrs['long_name']='Model monthly mean Ice concentration in march'
                ds_out['IceConceM03_mod'].attrs['units']='-'
                
                ds_out['IceConceM09_mod']= (('y','x'), zIceConc_S.values.astype('float32')) 
                ds_out['IceConceM09_mod'].attrs['long_name']='Model monthly mean Ice concentration in september'
                ds_out['IceConceM09_mod'].attrs['units']='-'
                
                ds_out['IceThick_obs']= (('y','x'), obs_thick.values.astype('float32')) 
                if zc_year >= 1979 and zc_year <= 2013 :
                        ds_out['IceThick_obs'].attrs['long_name']='PIOMAS annual mean Ice thickness over '+str(zc_year)
                else :
                        ds_out['IceThick_obs'].attrs['long_name']='PIOMAS climatological mean Ice thickness over 1979-2013'
                ds_out['IceThick_obs'].attrs['units']='m'
                
                ds_out['IceConceM03_obs']= (('yobs','xobs'), obs_conc_m03.values.astype('float32')) 
                if zc_year >= 1979 and zc_year <= 2015 :
                        ds_out['IceConceM03_obs'].attrs['long_name']='NSDIC monthly mean Ice concentration in march '+str(zc_year)
                else :
                        ds_out['IceConceM03_obs'].attrs['long_name']='NSDIC climatological mean Ice concentration in march over 1979-2015'
                ds_out['IceConceM03_obs'].attrs['units']='-'
                
                ds_out['IceConceM09_obs']= (('yobs','xobs'), obs_conc_m09.values.astype('float32')) 
                if zc_year >= 1979 and zc_year <= 2015 :
                        ds_out['IceConceM09_obs'].attrs['long_name']='NSDIC monthly mean Ice concentration in september '+str(zc_year)
                else :
                        ds_out['IceConceM09_obs'].attrs['long_name']='NSDIC climatological mean Ice concentration in september over 1979-2015'
                ds_out['IceConceM09_obs'].attrs['units']='-'
                
                ds_out['lat_obs']= (('yobs','xobs'), obs_lat.values.astype('float32')) 
                ds_out['lat_obs'].attrs['long_name']='Degrees north'
                ds_out['lat_obs'].attrs['units']='Deg'
                
                ds_out['lon_obs']= (('yobs','xobs'), obs_lon.values.astype('float32')) 
                ds_out['lon_obs'].attrs['long_name']='Degrees east'
                ds_out['lon_obs'].attrs['units']='Deg'
                
                ds_out['lat_mod']= (('y','x'), zlat.values.astype('float32')) 
                ds_out['lat_mod'].attrs['long_name']='Degrees north'
                ds_out['lat_mod'].attrs['units']='Deg'
                
                ds_out['lon_mod']= (('y','x'), zlon.values.astype('float32')) 
                ds_out['lon_mod'].attrs['long_name']='Degrees east'
                ds_out['lon_mod'].attrs['units']='Deg'

                ds_out = ds_out.set_coords(['lat_obs','lon_obs','lat_mod','lon_mod'])

                # Write the NetCDF file 
                ds_out.attrs['History'] = 'Diagnostics have been calculated using the Arctic monitoring tool '
                ds_out.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
                nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_ICEClim_'+'y'+str(zc_year)+'.nc'
                ds_out.to_netcdf(nc_f,engine='netcdf4')

        return

################################################################################################################################
def MLD_MAPSF( zlon, zlat, zMLD_M, zMLD_S, zMLDI_M, zMLDI_S, zCONF, zCASE, zc_year, zncout ) :
################################################################################################################################
        """
        Function dedicated to plot the Mixed Layer Depth 

        Input:
            zlon    : longitude 2D 
            zlat    : latitude 2D 
            zMLD_M  : March MLD 
            zMLD_S  : September MLD 
            zCONF   : the configuration name 
            zCASE   : the experiment name associated to the configuration
            zc_year : the current year 
            zncout  : logical to outputs (or not) results into a netcdf file 
        
        Output:
            None
        """
        # ----------------------------------------------------------------------

        plt.clf()
        fig, axes = plt.subplots( 2, 3, figsize=(11.7,8.3),  subplot_kw={'projection': ccrs.NorthPolarStereo(central_longitude=-60, true_scale_latitude=65)} )
        axes = axes.flatten()

	# Model MLD 
        # March mean MLD
        zMyvar='mldr10_1'
        DO_MAPS( zlon, zlat, zCONF, zCASE, zMLD_M, zMyvar, zc_year, seas='m03', ax=axes[0] )
        # September mean MLD
        zMyvar='mldr10_1'
        DO_MAPS( zlon, zlat, zCONF, zCASE, zMLD_S, zMyvar, zc_year, seas='m09', ax=axes[3] )

        # MLD from observation
        mld_obs_m03, mld_obs_m09, lon_obs, lat_obs = MLD_OBS()

        # March mean MLD
        zMyvar='mldr10_1'
        DO_MAPS( lon_obs, lat_obs, zCONF, zCASE, mld_obs_m03, zMyvar, zc_year, seas='m03', plot_obs=1, ax=axes[1] )
        # September mean MLD
        zMyvar='mldr10_1'
        DO_MAPS( lon_obs, lat_obs, zCONF, zCASE, mld_obs_m09, zMyvar, zc_year, seas='m09', plot_obs=1, ax=axes[4] )

	# Initial state MLD 
        # March mean MLD
        zMyvar='mldr10_1'
        DO_MAPS( zMLDI_M.lon, zMLDI_M.lat, zCONF, zCASE, zMLDI_M.values, zMyvar, zc_year, seas='m03', plot_obs=1, ax=axes[2], over_title='Init state \n m03' )
        # September mean MLD
        zMyvar='mldr10_1'
        DO_MAPS( zMLDI_S.lon, zMLDI_S.lat, zCONF, zCASE, zMLDI_S.values, zMyvar, zc_year, seas='m09', plot_obs=1, ax=axes[5], over_title='Init state \n m09' )
        plt.tight_layout()

        zfile_ext='_MLDClim_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zc_year)+'.png',dpi=300)

        if zncout:
                ds_out = xr.Dataset()
                
                # MLD fields
                #######################
                ds_out['MLDd01M03_mod']= (('y','x'), zMLD_M.values.astype('float32')) 
                ds_out['MLDd01M03_mod'].attrs['long_name']='Model monthly mean MLD in march based on a density criteria of 0.1 kg/m^3'
                ds_out['MLDd01M03_mod'].attrs['units']='m'

                ds_out['MLDd01M09_mod']= (('y','x'), zMLD_S.values.astype('float32')) 
                ds_out['MLDd01M09_mod'].attrs['long_name']='Model monthly mean MLD in september based on a density criteria of 0.1 kg/m^3'
                ds_out['MLDd01M09_mod'].attrs['units']='m'

                ds_out['MLDd01M03_woa']= (('ywoa','xwoa'), zMLDI_M.values.astype('float32')) 
                ds_out['MLDd01M03_woa'].attrs['long_name']='Climatological March MLD based on a density criteria of 0.1 kg/m^3'
                ds_out['MLDd01M03_woa'].attrs['units']='m'

                ds_out['MLDd01M09_woa']= (('ywoa','xwoa'), zMLDI_S.values.astype('float32')) 
                ds_out['MLDd01M09_woa'].attrs['long_name']='Climatological September MLD based on a density criteria of 0.1 kg/m^3'
                ds_out['MLDd01M09_woa'].attrs['units']='m'

                ds_out['MLDd01M03_obs']= (('yobs','xobs'), mld_obs_m03.values.astype('float32')) 
                ds_out['MLDd01M03_obs'].attrs['long_name']='MIMOC climatological mean in march over 2003-2014 based on a density criteria of 0.1 kg/m^3'
                ds_out['MLDd01M03_obs'].attrs['units']='m'

                ds_out['MLDd01M09_obs']= (('yobs','xobs'), mld_obs_m09.values.astype('float32')) 
                ds_out['MLDd01M09_obs'].attrs['long_name']='MIMOC climatological mean in september over 2003-2014 based on a density criteria of 0.1 kg/m^3'
                ds_out['MLDd01M09_obs'].attrs['units']='m'

                ds_out['lat_obs']= (('yobs','xobs'), lat_obs.astype('float32')) 
                ds_out['lat_obs'].attrs['long_name']='Degrees north'
                ds_out['lat_obs'].attrs['units']='Deg'
                
                ds_out['lon_obs']= (('yobs','xobs'), lon_obs.astype('float32')) 
                ds_out['lon_obs'].attrs['long_name']='Degrees east'
                ds_out['lon_obs'].attrs['units']='Deg'

                ds_out['lat_woa']= (('ywoa','xwoa'), zMLDI_S.lat.values.astype('float32')) 
                ds_out['lat_woa'].attrs['long_name']='Degrees north'
                ds_out['lat_woa'].attrs['units']='Deg'
                
                ds_out['lon_woa']= (('ywoa','xwoa'), zMLDI_S.lon.values.astype('float32')) 
                ds_out['lon_woa'].attrs['long_name']='Degrees east'
                ds_out['lon_woa'].attrs['units']='Deg'
                
                ds_out['lat_mod']= (('y','x'), zlat.values.astype('float32')) 
                ds_out['lat_mod'].attrs['long_name']='Degrees north'
                ds_out['lat_mod'].attrs['units']='Deg'
                
                ds_out['lon_mod']= (('y','x'), zlon.values.astype('float32')) 
                ds_out['lon_mod'].attrs['long_name']='Degrees east'
                ds_out['lon_mod'].attrs['units']='Deg'
                
                ds_out = ds_out.set_coords(['lat_obs','lon_obs','lat_mod','lon_mod'])

                # Write the NetCDF file 
                ds_out.attrs['History'] = 'Diagnostics have been calculated using the Arctic monitoring tool '
                ds_out.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
                nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_MLDClim_'+'y'+str(zc_year)+'.nc'
                ds_out.to_netcdf(nc_f,engine='netcdf4')

        return

################################################################################################################################
def DYN_MAPSF( zlon, zlat, zBSF, ds_eke, ds_topos, zdepth, zCONF, zCASE, zs_year, zncout ) :
################################################################################################################################
        """
        Function dedicated to plot the barotropic stream function as the EKE at the surface, 69m and 508m in the model

        Input:
            zlon       : longitude 2D 
            zlat       : latitude 2D 
            zBSF       : Barotropic function 
            ds_eke     : Dataset for the EKE 3D 
            ds_topos   : Dataset for 2D topostrophy
            zdepth     : 1D depth (used to get the depth at a given level)
            zCONF      : the configuration name 
            zCASE      : the experiment name associated to the configuration
            zs_year    : the current year 
            zncout     : logical to outputs (or not) results into a netcdf file 
        
        Output:
            None
        """
        # ----------------------------------------------------------------------

        ###########################################################################################
        # Plot barotropic stream function & topostrophy
        ###########################################################################################

        plt.clf()
        fig, axes = plt.subplots( 1, 2, figsize=(11.7,8.3),  subplot_kw={'projection': ccrs.NorthPolarStereo(central_longitude=-60, true_scale_latitude=65)} )
        axes = axes.flatten()

        # Annual mean PSI
        zMyvar='sobarstf'
        DO_MAPS( zlon, zlat, zCONF, zCASE, zBSF, zMyvar, str(zs_year), ax=axes[0] )

	# Annual mean topostrophy 
        zMyvar='topos'
        DO_MAPS( zlon, zlat, zCONF, zCASE, ds_topos['topos'], zMyvar, zs_year, plot_obs=0, ax=axes[1] )
        plt.tight_layout()

        zfile_ext='_DYNPSIClim_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zs_year)+'.png',dpi=300)

        ###########################################################################################
        # Plot EKE maps at the surface, ~70m and 500m
        ###########################################################################################

        plt.clf()
        fig, axes = plt.subplots( 3, 2, figsize=(8.3,11.7),  subplot_kw={'projection': ccrs.NorthPolarStereo(central_longitude=-60, true_scale_latitude=65)} )
        axes = axes.flatten()

        # Get the model depth at the surface , in the halocline and and the Atlantic water 
        lev0 =  0 # ~  0m in the model
        lev1 = 20 # ~ 69m in the model; corresponds to the halocline depth 
        lev2 = 39 # ~503m in the model; corresponds to the AW depth 
        zd0 = int(zdepth.isel(nav_lev=lev0).values.item()) ;  zd1 = int(zdepth.isel(nav_lev=lev1).values.item())  ; zd2 = int(zdepth.isel(nav_lev=lev2).values.item())

        # -----------  Get the EKE data from Obs.
        # -----------------------------------------------------------------------------------------
        # EKE from Von Appen et al. 
        obs_VonAppeneke = VONAPPEN_EKE_OBS( )

        # EKE from DOT observations (Armitage et al. 2017)
        obs_eke, lon_obs, lat_obs = EKE_OBS( t_year=zs_year )
        obs_eke = xr.where( obs_eke >= 9e20, npy.nan, obs_eke )

        # ----------- Mean EKE at the surface 
        # -----------------------------------------------------------------------------------------
        zMyvar='voeke'
        DO_MAPS( zlon, zlat, zCONF, zCASE, npy.log10(ds_eke[zMyvar].isel(z=lev0)), zMyvar, zs_year, slev=str(zd0) , ax=axes[0] )

        zMyvar='voeke'
        DO_MAPS( lon_obs, lat_obs, zCONF, zCASE, npy.log10(obs_eke), zMyvar, zs_year, slev=str(zd0), plot_obs=1, ax=axes[1] )

        # ----------- Mean EKE at 69m 
        # -----------------------------------------------------------------------------------------
        zMyvar='voeke'
        DO_MAPS( zlon, zlat, zCONF, zCASE, npy.log10(ds_eke[zMyvar].isel(z=lev1)), zMyvar, zs_year, slev=str(zd1) , ax=axes[2] )
	# Select data in the halocline depth
        mask_halo = npy.logical_and(obs_VonAppeneke['Mean depth'] >= 50, obs_VonAppeneke['Mean depth'] <= 100)
        lons_halo = obs_VonAppeneke.Longitude.where(mask_halo)
        lats_halo = obs_VonAppeneke.Latitude.where(mask_halo)
        eke_halo = obs_VonAppeneke.EKE.where(mask_halo)

        vmin = -6; vmax = -2
        cmap = plt.get_cmap('RdYlBu_r')
        axes[2].scatter(lons_halo, lats_halo, s = 10, c = npy.log10(eke_halo), cmap = cmap, vmin = vmin, vmax = vmax, edgecolors = 'k', linewidths = 0.5, transform=ccrs.PlateCarree() )

        # ----------- Mean EKE at 508m 
        # -----------------------------------------------------------------------------------------
        zMyvar='voeke'
        m = DO_MAPS( zlon, zlat, zCONF, zCASE, npy.log10(ds_eke[zMyvar].isel(z=lev2)), zMyvar, zs_year, slev=str(zd2), ax=axes[4] )
        # Select data in the AW layer 
        mask_aw = ~npy.isnan(obs_VonAppeneke['EKE at depth'])
        lons_aw = obs_VonAppeneke.Longitude.where(mask_aw)
        lats_aw = obs_VonAppeneke.Latitude.where(mask_aw)
        eke_aw = obs_VonAppeneke['EKE at depth'].where(mask_aw)

        vmin = -6; vmax = -2
        cmap = plt.get_cmap('RdYlBu_r')
        axes[4].scatter(lons_aw, lats_aw, s = 10, c = npy.log10(eke_aw), cmap = cmap, vmin = vmin, vmax = vmax, edgecolors = 'k', linewidths = 0.5, transform=ccrs.PlateCarree() )

        fig.delaxes(axes[3])
        fig.delaxes(axes[5])

        plt.tight_layout()

        zfile_ext='_DYNEKEClim_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zs_year)+'.png',dpi=300)

        if zncout:
                ds_out = xr.Dataset()
                
                # DYN fields
                #######################
                ds_out['PSI_mod']= (('y','x'), zBSF.values.astype('float32')) 
                ds_out['PSI_mod'].attrs['long_name']='Model annual mean barotropic streamfunction '
                ds_out['PSI_mod'].attrs['units']='Sv'

                ds_out['EKESurf_mod']= (('y','x'), ds_eke['voeke'].isel(z=lev0).values.astype('float32')) 
                ds_out['EKESurf_mod'].attrs['long_name']='Model annual mean EKE at the surface'
                ds_out['EKESurf_mod'].attrs['units']='m2/s2'

                ds_out['EKEz69_mod']= (('y','x'), ds_eke['voeke'].isel(z=lev1).values.astype('float32')) 
                ds_out['EKEz69_mod'].attrs['long_name']='Model annual mean EKE @ ~69m depth'
                ds_out['EKEz69_mod'].attrs['units']='m2/s2'

                ds_out['EKEz508_mod']= (('y','x'), ds_eke['voeke'].isel(z=lev2).values.astype('float32')) 
                ds_out['EKEz508_mod'].attrs['long_name']='Model annual mean EKE @ ~508m depth'
                ds_out['EKEz508_mod'].attrs['units']='m2/s2'

                ds_out['topos']= (('y','x'), ds_topos['topos'].values.astype('float32')) 
                ds_out['topos'].attrs['long_name']='Model annual mean topostrophy '
                ds_out['topos'].attrs['units']='m/s2'

                ds_out['EKESurf_obs']= (('yobs','xobs'), obs_eke.values.astype('float32')) 
                if zs_year >= 2003 and zs_year <= 2014 :
                        ds_out['EKESurf_obs'].attrs['long_name']='EKE annual mean derived from DOT (Armitage et al. 2017) in '+str(zs_year)
                else :
                        ds_out['EKESurf_obs'].attrs['long_name']='EKE climatological mean derived from DOT (Armitage et al. 2017) over 2003-2014'
                ds_out['EKESurf_obs'].attrs['units']='m2/s2'

                ds_out['lat_obs']= (('yobs','xobs'), lat_obs.values.astype('float32')) 
                ds_out['lat_obs'].attrs['long_name']='Degrees north'
                ds_out['lat_obs'].attrs['units']='Deg'
                
                ds_out['lon_obs']= (('yobs','xobs'), lon_obs.values.astype('float32')) 
                ds_out['lon_obs'].attrs['long_name']='Degrees east'
                ds_out['lon_obs'].attrs['units']='Deg'
                
                ds_out['lat_mod']= (('y','x'), zlat.values.astype('float32')) 
                ds_out['lat_mod'].attrs['long_name']='Degrees north'
                ds_out['lat_mod'].attrs['units']='Deg'
                
                ds_out['lon_mod']= (('y','x'), zlon.values.astype('float32')) 
                ds_out['lon_mod'].attrs['long_name']='Degrees east'
                ds_out['lon_mod'].attrs['units']='Deg'

                ds_out = ds_out.set_coords(['lat_obs','lon_obs','lat_mod','lon_mod'])

                # Write the NetCDF file 
                ds_out.attrs['History'] = 'Diagnostics have been calculated using the Arctic monitoring tool '
                ds_out.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
                nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_DYNClim_'+'y'+str(zs_year)+'.nc'
                ds_out.to_netcdf(nc_f,engine='netcdf4')
        return

################################################################################################################################
def TSD_MAPSF( zlon, zlat, zTemp, zSali, zTemp_IS, zSali_IS, zdepth, zCONF, zCASE, zc_year, zncout ) :
################################################################################################################################
        """
        Function dedicated to plot model temperature/salinity drift at the surface, ~100m, ~200m and ~300m

        Input:
            zlon     : longitude 2D  
            zlat     : latitude 2D  
            zTemp    : Temperature 3D
            zSali    : Salinity 3D 
            zTemp_IS : Initial temperature 3D 
            zSali_IS : Initial salinity 3D 
            zdepth   : 1D depth (used to get the depth at a given level)
            zCONF    : configuration name 
            zCASE    : experiment name associated to the configuration
            zc_year  : current year 
            zncout   : logical to outputs (or not) results into a netcdf file 
        
        Output:
            None
        """
        # ----------------------------------------------------------------------

        # Get the model depth at the surface and at ~100m 
        zd1 = int(zdepth.isel(nav_lev=0).values.item())  ; zd2 = int(zdepth.isel(nav_lev=23).values.item())
        #zd1 = npy.round(zdepth.isel(nav_lev=0).values).item()  ; zd2 = npy.round(zdepth.isel(nav_lev=23).values).item()

        plt.clf()
        fig = plt.figure()
        projection = ccrs.NorthPolarStereo(central_longitude=-60, true_scale_latitude=65)

        num_fram=220
        # Surface temperature
        zMyvar='votemper'   ; fram=num_fram+1
        ax = fig.add_subplot(fram, projection=projection)
        DO_MAPS( zlon, zlat, zCONF, zCASE, npy.squeeze(zTemp[0,:,:]-zTemp_IS[0,:,:]), zMyvar, zc_year, slev=str(zd1) , ano=1, ax=ax )
        # ~100m temperature
        zMyvar='votemper'   ; fram=num_fram+2
        ax = fig.add_subplot(fram, projection=projection)
        DO_MAPS( zlon, zlat, zCONF, zCASE, npy.squeeze(zTemp[23,:,:]-zTemp_IS[23,:,:]), zMyvar, zc_year, slev=str(zd2), ano=1, ax=ax )
        # Surface salinity
        zMyvar='vosaline'   ; fram=num_fram+3
        ax = fig.add_subplot(fram, projection=projection)
        DO_MAPS( zlon, zlat, zCONF, zCASE, npy.squeeze(zSali[0,:,:]-zSali_IS[0,:,:]), zMyvar, zc_year, slev=str(zd1), ano=1, ax=ax )
        # ~100m  salinity
        zMyvar='vosaline'   ; fram=num_fram+4
        ax = fig.add_subplot(fram, projection=projection)
        DO_MAPS( zlon, zlat, zCONF, zCASE, npy.squeeze(zSali[23,:,:]-zSali_IS[23,:,:]), zMyvar, zc_year, slev=str(zd2), ano=1, ax=ax )
        plt.tight_layout()

        zfile_ext='_TSDIffClim_@'+str(zd1)+'m@'+str(zd2)+'m_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zc_year)+'.png',dpi=300)

        plt.clf()

        # Get the model depth at ~200m and at ~300m 
        zd3 = int(zdepth.isel(nav_lev=30).values.item())  ; zd4 = int(zdepth.isel(nav_lev=34).values.item())

        num_fram=220
        # ~200m temperature
        zMyvar='votemper'   ; fram=num_fram+1
        ax = fig.add_subplot(fram, projection=projection)
        DO_MAPS( zlon, zlat, zCONF, zCASE, npy.squeeze(zTemp[30,:,:]-zTemp_IS[30,:,:]), zMyvar, zc_year, slev=str(zd3) , ano=1, ax=ax )
        # ~300m temperature
        zMyvar='votemper'   ; fram=num_fram+2
        ax = fig.add_subplot(fram, projection=projection)
        DO_MAPS( zlon, zlat, zCONF, zCASE, npy.squeeze(zTemp[34,:,:]-zTemp_IS[34,:,:]), zMyvar, zc_year, slev=str(zd4), ano=1, ax=ax )
        # ~200m salinity
        zMyvar='vosaline'   ; fram=num_fram+3
        ax = fig.add_subplot(fram, projection=projection)
        DO_MAPS( zlon, zlat, zCONF, zCASE, npy.squeeze(zSali[30,:,:]-zSali_IS[30,:,:]), zMyvar, zc_year, slev=str(zd3), ano=1, ax=ax )
        # ~300m  salinity
        zMyvar='vosaline'   ; fram=num_fram+4
        ax = fig.add_subplot(fram, projection=projection)
        DO_MAPS( zlon, zlat, zCONF, zCASE, npy.squeeze(zSali[34,:,:]-zSali_IS[34,:,:]), zMyvar, zc_year, slev=str(zd4), ano=1, ax=ax )
        plt.tight_layout()

        zfile_ext='_TSDIffClim_@'+str(zd3)+'m@'+str(zd4)+'m_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zc_year)+'.png',dpi=300)

        return

################################################################################################################################
def ATL_MAPSF( zlon, zlat, zMLD_M, zMLD_year, zTemp, zTemp_IS, zSSH, zdepth, zCONF, zCASE, zc_year, zs_year, ze_year, zncout ) :
################################################################################################################################
        """
        Function dedicated to plot Atlantic MLD, SSH in differents areas such as GIN, LAB, IRM seas

        Input:
            zlon      : longitude 2D  
            zlat      : latitude 2D  
            zMLD_M    : March MLD  
            zMLD_year : Full year MLD (for time-series at a mooring location)
            zTemp     : Temperature 3D  
            zTemp_IS  : Initial temperature 3D  
            zSSH      : Initial salinity 3D 
            zdepth    : 1D depth (used to get the depth at a given level)
            zCONF     : the configuration name 
            zCASE     : the experiment name associated to the configuration
            zc_year   : the current year 
            zs_year   : the first year 
            ze_year   : the last year 
            zncout    : logical to outputs (or not) results into a netcdf file 
        
        Output:
            None
        """
        # ----------------------------------------------------------------------

        # Define time axis for plot
        time_grid = npy.arange(zs_year,ze_year+2,1.,dtype=int)
        newlocsx = npy.array(time_grid,'f')
        newlabelsx = npy.array(time_grid,'i')

        # Set the time axis
        t_months = (npy.arange(12)*30.+15.)/365.
        time_axis = npy.tile(zs_year,12)+t_months

        # MLD IN THE LABRADOR SEA IN MARCH
        ###################################
        plt.figure()

        num_fram=110
        # March mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+1
        my_cblab=r'MLD (m)'   ;   my_cmap=plt.get_cmap('Blues')
        ztitle=zCASE +' mean MLD01 \n'+str(zc_year)+'  m03'
        vmin=0. ; vmax=2400. ; vint=100.   ;   contours=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]
        limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]

        plt.subplot(fram)
        zoutmap = BATHY_MAP( ztype='isol1000',zarea='labsea' )
        PROJ_PLOT( zlon, zlat, zMLD_M, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='labsea' )
        plt.tight_layout()

        zfile_ext='_LAB_MLDClimm03_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zc_year)+'.png',dpi=300)

        # Plot the Time-serie for MLD at a specific location K1 mooring in the Labrador Sea and in Irminger Sea
        # After Schott et al. DSRI2009 56.33N, -52.40W
        plt.clf()
        plt.figure()
        if zCONF == 'CREG12.L75' : 
                i_K1=518   ;   j_K1=502   # CREG12.L75 C-type indices
        else :
                i_K1=173   ;   j_K1=168   # CREG025.L75 C-type indices
        ax=plt.subplot(211)
        # In Lab. Sea
        plt.plot( time_axis, -1.*npy.squeeze(zMLD_year[:,j_K1:j_K1+1,i_K1:i_K1+1]), linewidth=0.8, color='k', label='Lab Sea K1' )
        # Plot obs. MLD in March
        year_obs = npy.arange(1995,2006,1)+0.20547945   ; mld_obs=[-2300.,-1300.,-1400.,-1000.,-1000.,-1100.,-1100.,-1200.,-1400.,-700.,-1300.]
        plt.scatter(year_obs,mld_obs)

        # In Irm. Sea
        if zCONF == 'CREG12.L75' : 
                i_K1=697   ;   j_K1=577   # CREG12.L75 C-type indices geo loc   60.88N  -36.99W
        else :
                i_K1=232   ;   j_K1=192   # CREG025.L75 C-type indices geo loc   60.88N  -36.99W
        plt.plot( time_axis, -1.*npy.squeeze(zMLD_year[:,j_K1:j_K1+1,i_K1:i_K1+1]), linewidth=0.8, color='g', label='Irm Sea ')
        plt.title( zCASE+' MLD 0.01 in Lab. & Irm. Seas' , size=9 )
        plt.ylabel( 'Mean depth \n'+r'(m)', size=7 )
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
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zs_year)+'LASTy.png',dpi=300)

        # Add an artificial mooring within the deepest convection area 
        # -54W 58N
        #  dl_dis=    1.634 km
        #      507       507       541       541
        # -54.0272  -54.0272   57.9970   57.9970
        plt.clf()
        plt.figure()
        if zCONF == 'CREG12.L75' : 
                i_K1=506   ;   j_K1=540   # CREG12.L75 C-type indices
        else : 
                i_K1=169   ;   j_K1=181   # CREG025.L75 C-type indices
        ax=plt.subplot(211)
        # In Lab. Sea
        plt.plot( time_axis,-1.*npy.squeeze(zMLD_year[:,j_K1:j_K1+1,i_K1:i_K1+1]),linewidth=0.8, color='k', label='Lab Sea DeepConv' )
        # Plot obs. MLD in March
        year_obs = npy.arange(1995,2006,1)+0.20547945   ; mld_obs=[-2300.,-1300.,-1400.,-1000.,-1000.,-1100.,-1100.,-1200.,-1400.,-700.,-1300.]
        plt.scatter(year_obs,mld_obs)
        plt.title(zCASE+' MLD 0.01 in Lab. @ -54W,58N ',size=9)
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
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zs_year)+'LASTy.png',dpi=300)

        plt.clf()
        plt.figure()
        # MLD IN THE IRMINGER SEA IN MARCH
        ###################################
        num_fram=110
        # March mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+1
        my_cblab=r'MLD (m)'   ;   my_cmap=plt.get_cmap('Blues')
        ztitle=zCASE +' mean MLD01 \n'+str(zc_year)+'  m03'
        vmin=0. ; vmax=1600. ; vint=100.   ;   contours=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
        limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
        plt.subplot(fram)
        zoutmap = BATHY_MAP( ztype='isol1000',zarea='irmsea' )
        PROJ_PLOT( zlon, zlat, zMLD_M, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='irmsea' )
        plt.tight_layout()

        zfile_ext='_IRM_MLDClimm03_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zc_year)+'.png',dpi=300)

        plt.clf()
        plt.figure()
        # MLD IN THE GIN SEAS IN MARCH
        ###################################
        num_fram=110
        # March mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+1
        my_cblab=r'MLD (m)'   ;   my_cmap=plt.get_cmap('Blues')
        ztitle=zCASE +' mean MLD01 \n'+str(zc_year)+'  m03'
        vmin=0. ; vmax=1600. ; vint=100.   ;   contours=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
        limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
        plt.subplot(fram)
        zoutmap = BATHY_MAP( ztype='isol1000',zarea='ginsea' )
        PROJ_PLOT( zlon, zlat, zMLD_M, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='ginsea' )
        plt.tight_layout()

        zfile_ext='_GIN_MLDClimm03_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zc_year)+'.png',dpi=300)

        plt.clf()
        plt.figure()
        # PLOT ISOTHERM 17 Deg OFF CAPE HATTERAS
        ########################################
        num_fram=110
        zMyvar='votemper'   ; fram=num_fram+1
        my_cblab=r'ISO 17 (DegC)'   ;   my_cmap=plt.get_cmap('jet')
        ztitle=zCASE +' mean Iso 17 DegC \n'+str(zc_year)
        vmin=0. ; vmax=2400. ; vint=100.   ;   contours=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]
        limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]

        plt.subplot(fram)
        klev=29
        zoutmap = BATHY_MAP( ztype='isol1000', zarea='GulfS' )
        zzlon = zlon.values   ; zzlat = zlat.values 
        PROJ_PLOT( zzlon, zzlat, npy.squeeze(zTemp[klev,:,:]), contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='GulfS' )
        PROJ_PLOT( zzlon, zzlat, npy.squeeze(zTemp_IS[klev,:,:]), contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='GulfS', data_ref=True )
        plt.tight_layout()

        zfile_ext='_ATL_ISO17Clim_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zc_year)+'.png',dpi=300)

        plt.clf()
        plt.figure()
        # PLOT SSH OVER THE ATLANTIC AREA
        #################################
        num_fram=110
        zMyvar='sossheig'   ; fram=num_fram+1
        my_cblab=r'SSH (cm)'   ;   my_cmap=plt.get_cmap('coolwarm')
        ztitle=zCASE +' mean SSH \n'+str(zc_year)
        vmin=-100. ; vmax=100. ; vint=5.  ;   contours=npy.arange(vmin,vmax+vint,vint)
        limits=[vmin,vmax,vint]           ;   myticks=npy.arange(vmin,vmax+vint,vint)

        zoutmap = BATHY_MAP( ztype='isol1000', zarea='natl' )
        PROJ_PLOT( zzlon,  zzlat, zSSH*100. , contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='natl' )
        #plt.tight_layout()

        zfile_ext='_ATL_SSHClim_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zc_year)+'.png',dpi=300)

        return

################################################################################################################################
def MOC_MAPSF( zlon, zlat, zAMOC, zdepth, zCONF, zCASE, zc_year, zs_year, ze_year, zncout ) :
################################################################################################################################
        """
        Function dedicated to plot AMOC and the maximum time-series 

        Input:
            zlon      : longitude 2D  
            zlat      : latitude 2D  
            zAMOC     : AMOC yearly mean 
            zdepth    : 1D depth (used to get the depth at a given level)
            zCONF     : the configuration name 
            zCASE     : the experiment name associated to the configuration
            zc_year   : the current year 
            zs_year   : the first year 
            ze_year   : the last year 
            zncout    : logical to outputs (or not) results into a netcdf file 
        
        Output:
            None
        """
        # ----------------------------------------------------------------------

        plt_AMOCTS=True

        if  plt_AMOCTS: 
                # AMOC 
                #######

                # Prepare 2 dimnsional (lat,depth) array for plotting   
                lat2Dz = npy.reshape(zlat,(zlat.size,1)).T
                ypltz = npy.repeat(lat2Dz,zdepth.shape[0],axis=0)
                locpath='./'
                locfile='Bathymetry.nc'
                if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
                        ds_msk = xr.open_dataset(locpath+locfile)[['nav_lat']]
                        full_lat = ds_msk['nav_lat']
                if zCONF == 'CREG12.L75' : 
                	select_ylat=full_lat[:,862]
                else: 
                	select_ylat=full_lat[:,288]
                select_ylat_reshape = npy.reshape(select_ylat.values,(select_ylat.size,1))
                ypltz = npy.repeat(select_ylat_reshape,zdepth.shape[0],axis=1).T

                z2dt = npy.reshape(zdepth.values,(zdepth.size,1))
                zplt = npy.repeat(z2dt,zlat.shape[0],axis=1)

                # Make the plot 
                plt.figure()

                my_cblab=r'AMOC (Sv)'   ;   my_cmap=plt.get_cmap('jet')
                ztitle=zCASE +' mean AMOC \n'+str(zc_year)
                vmin=-15. ; vmax=15. ; vint=1.    ;   contours=npy.arange(vmin,vmax+vint,vint)
                limits=[vmin,vmax,vint]  ;             myticks=npy.arange(vmin,vmax+vint,vint)
                norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])

                num_fram=210
                fram=num_fram+1
                plt.subplot(fram)
                plt.title( ztitle, fontsize=6 )
                plt.contourf( ypltz, zplt*(-1.e-3), npy.squeeze(zAMOC), contours, cmap=my_cmap, norm=norm, extend='both' )
                plt.xlabel('Latitude',size=6)
                plt.ylabel('Depth (kms)',size=6)
                plt.xticks(fontsize=6)
                plt.yticks(fontsize=6)
                contours = npy.arange(vmin,vmax+vint,2.*vint)
                C = plt.contour( ypltz, zplt*(-1.e-3), npy.squeeze(zAMOC), linewidths=0.5, levels=contours, colors='k' )
                plt.clabel( C, C.levels, inline=True, fmt='%3.0f', fontsize=6 )
                #plt.tight_layout()

                zfile_ext='_AMOCClim_'
                plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zc_year)+'.png',dpi=300)

                if zncout:
                	ds_out = xr.Dataset()
                	
                	ds_out.coords['latitude'] = (('z','y'), ypltz.astype('float32')) 
                	ds_out.coords['latitude'].attrs['long_name'] = 'latitude'
                	ds_out.coords['latitude'].attrs['units'] = 'degrees_north'
                	
                	ds_out.coords['Depth'] = (('z','y'), -1*zplt.astype('float32')) 
                	ds_out.coords['Depth'].attrs['long_name'] = 'Depth'
                	ds_out.coords['Depth'].attrs['units'] = 'm'
                
                	# Save diags. AMOC
                	timevalue = pd.to_datetime(str(zc_year)+'-06-30')
                	ds_out.coords['time'] = (('time'), [timevalue] )
                	
                	ds_out['AMOC'] = (('z','y'), zAMOC.values.astype('float32')) 
                	ds_out['AMOC'].attrs['long_name'] = 'Atlantic Meridional Overturning Circulation'
                	ds_out['AMOC'].attrs['units'] = 'Sv'
                	
                	# Write the NetCDF file 
                	ds_out.attrs['History'] = 'Diagnostics have been calculated using the Arctic monitoring tool '
                	ds_out.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
                	nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_AMOC_y'+str(zc_year)+'.nc'
                	ds_out.to_netcdf(nc_f,engine='netcdf4')
                
                #------------------------------------------------------------------------------------------------------------------------
                ########################################
                # Plot LONG TIME-SERIES
                ########################################
                #------------------------------------------------------------------------------------------------------------------------
                
                if ze_year-zs_year+1 > 1 :
                
                        print()
                        print('				##################################################################  ')
                        print('				##################################################################  ')
                        print('				##################      PLOT AMOC LONG TIME-SERIES      ##########  ')
                        print('				##################################################################  ')
                        print('				##################################################################  ')
                        print()
                        
                        # Read model data set 
                        #######################################################################################################
                        locpath = './NETCDF/'
                        locfile = zCONF+'-'+zCASE+'_AMOC_y????.nc'
                        ds_amoc = xr.open_mfdataset(locpath+locfile, engine='netcdf4', concat_dim=['time'], combine='nested', parallel=True)
                        mod_time_axis = ds_amoc.time.dt.year.values + (ds_amoc.time.dt.month.values - 1 + 0.5) / 12
                
                        # Max. AAMOC Time-series @ 40N or 43N
                        #######################################
                        time_grid = npy.arange(zs_year,ze_year+1,1.,dtype=int)
                        newlocsx  = npy.array(time_grid,'f')
                        newlabelsx = npy.array(time_grid,'i')

                        if zCONF == 'CREG025.L75' :
                              jloc= 205 # Lat 40N
                              #jloc= 252 # Lat 43N
                              if jloc == 205   : zlat='L40N'
                              elif jloc == 252 : zlat='L43N'
                        else :
                              print(' For the configuration '+zCONF+' the AMOC latitude index must be hard coded in MOC_MAPSF function')

                        plt.clf()
                        ax=plt.subplot(211)
                        plt.plot( time_grid, ds_amoc['AMOC'].isel(y=jloc).max(dim='z').values, 'k',linewidth=0.5, label=zCONF )
                        plt.title( zCASE+' Max AMOCz ', size=8 )
                        plt.ylabel( 'Max AMOCz (Sv) \n'+' @ '+zlat, size=7 )
                        plt.xlim([1978.,2025.])
                        plt.ylim([5.,15.])
                        plt.xticks(newlocsx,newlabelsx,size=5)
                        #plt.setp(ax.get_xticklabels(),rotation=90)
                        plt.setp(ax.get_xticklabels(),visible=False)
                        plt.yticks(size=6)
                        plt.grid(True, linewidth=0.7,linestyle='--',alpha=0.7,color='grey')

                        ax=plt.subplot(212)
                        plt.plot( time_grid, ds_amoc['AMOC'].max(dim=('z','y')), 'k',linewidth=0.5, label=zCONF )
                        plt.ylabel( 'Max AMOCz (Sv)', size=7 )
                        plt.xlim([1978.,2025.])
                        plt.ylim([10.,20.])
                        plt.xticks(newlocsx,newlabelsx,size=5)
                        plt.setp(ax.get_xticklabels(),rotation=90)
                        plt.yticks(size=6)
                        plt.grid(True, linewidth=0.7,linestyle='--',alpha=0.7,color='grey')
                        plt.tight_layout()

                        zfile_ext='_MaxAMOCz_'+zlat+'_LGTS_'
                        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zs_year)+'LASTy.png',dpi=300)

        return

################################################################################################################################
def MTS_MAPSF( zlon, zlat, zCONF, zCASE, zMLD_M, zMLD_S, zTS_M, zTS_S, zgdepw_0, ze3t_0, zc_year, zncout ) :
################################################################################################################################
	"""
	Function dedicated to plot mean temperature/salinity in the MLD 
	
	Input:
	    zlon     : longitude 2D  
	    zlat     : latitude 2D  
	    zCONF    : configuration name 
	    zCASE    : experiment name associated to the configuration
	    zMLD_M   : March MLD  
	    zMLD_S   : September MLD  
	    zTS_M    : March temperature/Salinity Dataset 3D 
	    zTS_S    : September temperature/Salinity Dataset 3D 
	    zgdepw_0 : 3D depth  
	    ze3t_0   : 3D vertical scale factor 
	    zc_year  : current year 
	    zncout   : logical to outputs (or not) results into a netcdf file 
	
	Output:
	    None
	"""
	# ----------------------------------------------------------------------

	# Mask all levels below the MLD 
	e3t_0msk_SeasM= xr.where( zgdepw_0 < zMLD_M, ze3t_0, 0.)
	e3t_0msk_SeasS= xr.where( zgdepw_0 < zMLD_S, ze3t_0, 0.)

	# Sum all e3t scale factors over the vertcal axis
	e3t_0sum_SeasM = e3t_0msk_SeasM.sum(dim='z').squeeze()
	e3t_0sum_SeasS = e3t_0msk_SeasS.sum(dim='z').squeeze()

	# Compute the T/S mean within the ML in March/September
	T_mldM = (e3t_0msk_SeasM * zTS_M['votemper']).sum(dim='z').squeeze()/e3t_0sum_SeasM
	S_mldM = (e3t_0msk_SeasM * zTS_M['vosaline']).sum(dim='z').squeeze()/e3t_0sum_SeasM
                                                       
	T_mldS = (e3t_0msk_SeasS * zTS_S['votemper']).sum(dim='z').squeeze()/e3t_0sum_SeasS
	S_mldS = (e3t_0msk_SeasS * zTS_S['vosaline']).sum(dim='z').squeeze()/e3t_0sum_SeasS

	# MLTS from MIMOC observations
	mlT_obs, mlS_obs, lon_obs, lat_obs = MLTS_OBS()

	# Plots Temperature maps 
	########################
	plt.clf()
	fig = plt.figure()
	projection = ccrs.NorthPolarStereo(central_longitude=-60, true_scale_latitude=65)

	vmin=-2. ; vmax=3. ; vint=0.1
	contours=npy.arange(vmin,vmax+vint,vint)  
	limits=[vmin,vmax,vint]			 
	myticks=npy.arange(vmin,vmax+vint,vint) 
	my_cblab=r'($^\circ$C)'
	my_cmap= plt.get_cmap('Spectral_r')
	
	fram=231
	ax = fig.add_subplot(fram, projection=projection)
	ztitle=zCASE +' March MLT \n'+str(zc_year)
	zoutmap = BATHY_MAP( ztype='isol1000', ax=ax )
	PROJ_PLOT( zlon, zlat, T_mldM, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, ax=ax )

	fram=234
	ax = fig.add_subplot(fram, projection=projection)
	ztitle=zCASE +' September MLT \n'+str(zc_year)
	zoutmap = BATHY_MAP( ztype='isol1000', ax=ax )
	PROJ_PLOT(zlon, zlat, T_mldS, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, ax=ax )

	fram=232
	ax = fig.add_subplot(fram, projection=projection)
	ztitle='MIMOC March MLT'
	zoutmap = BATHY_MAP( ztype='isol1000', ax=ax )
	PROJ_PLOT( lon_obs, lat_obs, mlT_obs[2,:,:].squeeze(), contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zplot_obs=1, ax=ax )

	fram=235
	ax = fig.add_subplot(fram, projection=projection)
	ztitle='MIMOC September MLT'
	zoutmap = BATHY_MAP( ztype='isol1000', ax=ax )
	PROJ_PLOT( lon_obs, lat_obs, mlT_obs[8,:,:].squeeze(), contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zplot_obs=1, ax=ax )

	plt.tight_layout()

	zfile_ext='_MTSClimT_'
	plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zc_year)+'.png',dpi=300)

	# Plots Salinity maps 
	##################### 
	plt.clf()
	vmin=26. ; vmax=36. ; vint=0.5
	contours=npy.arange(vmin,vmax+vint,vint)  
	limits=[vmin,vmax,vint]			 
	myticks=npy.arange(vmin,vmax+2.*vint,2.*vint) 
	my_cblab=r'(PSU)'
	my_cmap= plt.get_cmap('Spectral_r')

	fram=231
	ax = fig.add_subplot(fram, projection=projection)
	ztitle=zCASE +' March MLS \n'+str(zc_year)
	zoutmap = BATHY_MAP( ztype='isol1000', ax=ax )
	PROJ_PLOT( zlon, zlat, S_mldM, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar='MLTSS', ax=ax )

	fram=234
	ax = fig.add_subplot(fram, projection=projection)
	ztitle=zCASE +' September MLS \n'+str(zc_year)
	zoutmap = BATHY_MAP( ztype='isol1000', ax=ax )
	PROJ_PLOT( zlon, zlat, S_mldS, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar='MLTSS', ax=ax )

	fram=232
	ax = fig.add_subplot(fram, projection=projection)
	ztitle='MIMOC March MLS'
	zoutmap = BATHY_MAP( ztype='isol1000', ax=ax )
	PROJ_PLOT( lon_obs, lat_obs, mlS_obs[2,:,:].squeeze(), contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar='MLTSS', zplot_obs=1, ax=ax )

	fram=235
	ax = fig.add_subplot(fram, projection=projection)
	ztitle='MIMOC September MLS'
	zoutmap = BATHY_MAP( ztype='isol1000', ax=ax )
	PROJ_PLOT( lon_obs, lat_obs, mlS_obs[8,:,:].squeeze(), contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar='MLTSS', zplot_obs=1, ax=ax )

	plt.tight_layout()

	zfile_ext='_MTSClimS_'
	plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zc_year)+'.png',dpi=300)

	if zncout:
		ds_out = xr.Dataset()
		
		# ML T/S mean field
		#######################
		ds_out['S_mldM_march']= (('y','x'), S_mldM.values.astype('float32')) 
		ds_out['S_mldM_march'].attrs['long_name']='Model ML mean S in March'
		ds_out['S_mldM_march'].attrs['units']='PSU'
                
		ds_out['S_mldM_septe']= (('y','x'), S_mldS.values.astype('float32')) 
		ds_out['S_mldM_septe'].attrs['long_name']='Model ML mean S in September'
		ds_out['S_mldM_septe'].attrs['units']='PSU'
                
		ds_out['T_mldM_march']= (('y','x'), T_mldM.values.astype('float32')) 
		ds_out['T_mldM_march'].attrs['long_name']='Model ML mean T in March'
		ds_out['T_mldM_march'].attrs['units']='DegC'
                
		ds_out['T_mldM_septe']= (('y','x'), T_mldS.values.astype('float32')) 
		ds_out['T_mldM_septe'].attrs['long_name']='Model ML mean S in September'
		ds_out['T_mldM_septe'].attrs['units']='DegC'
                
		ds_out['mlS_obs_march']= (('yobs','xobs'), mlS_obs[2,:,:].values.squeeze().astype('float32')) 
		ds_out['mlS_obs_march'].attrs['long_name']='MIMOC ML mean S in March'
		ds_out['mlS_obs_march'].attrs['units']='PSU'
                
		ds_out['mlS_obs_septe']= (('yobs','xobs'), mlS_obs[8,:,:].values.squeeze().astype('float32')) 
		ds_out['mlS_obs_septe'].attrs['long_name']='MIMOC ML mean S in September'
		ds_out['mlS_obs_septe'].attrs['units']='PSU'
                
		ds_out['mlT_obs_march']= (('yobs','xobs'), mlT_obs[2,:,:].values.squeeze().astype('float32')) 
		ds_out['mlT_obs_march'].attrs['long_name']='MIMOC ML mean T in March'
		ds_out['mlT_obs_march'].attrs['units']='DegC'
                
		ds_out['mlT_obs_septe']= (('yobs','xobs'), mlT_obs[8,:,:].values.squeeze().astype('float32')) 
		ds_out['mlT_obs_septe'].attrs['long_name']='MIMOC ML mean T in September'
		ds_out['mlT_obs_septe'].attrs['units']='DegC'
                
		ds_out['lat_obs']= (('yobs','xobs'), lat_obs.values.astype('float32')) 
		ds_out['lat_obs'].attrs['long_name']='Degrees north'
		ds_out['lat_obs'].attrs['units']='Deg'
		
		ds_out['lon_obs']= (('yobs','xobs'), lon_obs.values.astype('float32')) 
		ds_out['lon_obs'].attrs['long_name']='Degrees east'
		ds_out['lon_obs'].attrs['units']='Deg'
		
		ds_out['lat_mod']= (('y','x'), zlat.values.astype('float32')) 
		ds_out['lat_mod'].attrs['long_name']='Degrees north'
		ds_out['lat_mod'].attrs['units']='Deg'
		
		ds_out['lon_mod']= (('y','x'), zlon.values.astype('float32')) 
		ds_out['lon_mod'].attrs['long_name']='Degrees east'
		ds_out['lon_mod'].attrs['units']='Deg'

		ds_out = ds_out.set_coords(['lat_obs','lon_obs','lat_mod','lon_mod'])

		# Write the NetCDF file 
		ds_out.attrs['History'] = 'Diagnostics have been calculated using the Arctic monitoring tool '
		ds_out.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
		nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_MTSClim_'+'y'+str(zc_year)+'.nc'
		ds_out.to_netcdf(nc_f,engine='netcdf4')

	return

################################################################################################################################
def AWT_MAPSF( zlon, zlat, zTemp, zSali, zTempIS, zSaliIS, zCONF, zCASE, zc_year, zncout) :
################################################################################################################################
	"""
	Function dedicated to plot the Atlantic Water maximum temperature as its associated depth
	
	Input:
	    zlon     : longitude 2D  
	    zlat     : latitude 2D  
	    zTemp    : Temperature 3D 
	    zSali    : Salinity 3D 
	    zTempIS  : Initial state Temperature 3D 
	    zSaliIS  : Initial state Salinity 3D 
	    zCONF    : configuration name 
	    zCASE    : experiment name associated to the configuration
	    zc_year  : current year 
	    zncout   : logical to outputs (or not) results into a netcdf file 
	
	Output:
	    None
	"""
	# ----------------------------------------------------------------------

	######################################################################################
	# Start with model output first 

        # To keep Temp where S > 33.5 only (away from the surface)
	Smask335 = (zSali > 33.5)
	Temp = zTemp.where(Smask335)
	Temp_filled = (Temp.fillna(-10)).compute()

	# Find the indices over z where T is max.
	depth_map = Temp_filled.argmax(dim='z',skipna=True)
        # Now get the effective depth of T max 
	deptht = zTemp['z']
	zAWTmax_depth1 = deptht[depth_map.compute()]
	zAWTmax_depth1 = xr.where( npy.isnan(zSali[0,:,:]), npy.nan, zAWTmax_depth1 )

	# Find the T max value over z 
	temp_map = Temp_filled.max(dim='z',skipna=True)
	mask = (temp_map.compute()>-10)
	zAWTmax1 = temp_map.where(mask,npy.nan)

	######################################################################################
	# Initial state 

        # To keep Temp where S > 33.5 only (away from the surface)
	SmaskIS335 = (zSaliIS > 33.5)
	TempIS = zTempIS.where(SmaskIS335)
	TempIS_filled = (TempIS.fillna(-10)).compute()

	# Find the indices over z where T is max.
	depth_mapIS = TempIS_filled.argmax(dim='z',skipna=True)
        # Now get the effective depth of T max 
	depthtIS = zTempIS['z']
	zAWTmax_depthIS = depthtIS[depth_mapIS.compute()]
	zAWTmax_depthIS = xr.where( npy.isnan(zSaliIS[0,:,:]), npy.nan, zAWTmax_depthIS )
	zAWTmax_depthIS = xr.where( zSaliIS[0,:,:] < 1., npy.nan, zAWTmax_depthIS )

	# Find the T max value over z 
	temp_mapIS = TempIS_filled.max(dim='z',skipna=True)
	maskIS = (temp_mapIS.compute()>-10)
	zAWTmaxIS = temp_mapIS.where(maskIS,npy.nan)
	zAWTmaxIS = xr.where( zSaliIS[0,:,:] <= 0., npy.nan, zAWTmaxIS )

	# Make the plot for the AW Max Temp 
	#############################################################################################
	vmin=0. ; vmax=3. ; vint=0.2
	contours=npy.arange(vmin,vmax+vint,vint)  
	limits=[vmin,vmax,vint]			 
	myticks=npy.arange(vmin,vmax+vint,vint) 
	
	my_cblab=r'($^\circ$C)'
	my_cmap= plt.get_cmap('Spectral_r')

	plt.clf()
	fig = plt.figure()
	projection = ccrs.NorthPolarStereo(central_longitude=-60, true_scale_latitude=65)

	ax1 = fig.add_subplot(221, projection=projection)
	ztitle=zCASE +' AW Max Temp \n'+str(zc_year)
	zoutmap = BATHY_MAP( ztype='isol1000', ax=ax1 )
	zMyvar = 'votemper'	
	PROJ_PLOT( zlon, zlat, zAWTmax1, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, ax=ax1 )

	ax2 = fig.add_subplot(222, projection=projection)
	ztitle=' AW Max Temp from \n'+' WOA2009'
	zoutmap = BATHY_MAP( ztype='isol1000', ax=ax2 )
	zMyvar = 'votemper'	
	PROJ_PLOT( zlon, zlat, zAWTmaxIS, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zplot_obs=0, ax=ax2 ) 
	
	# Make the plot for the AW Max Temp depth
	#############################################################################################
	vmin=0. ; vmax=800. ; vint=50.
	contours=npy.arange(vmin,vmax+vint,vint)  
	limits=[vmin,vmax,vint]			 
	myticks=npy.arange(vmin,vmax+vint,vint) 
	
	ztitle=zCASE+' AW Max Temp depth '
	my_cblab=r'(m)'
	my_cmap= plt.get_cmap('jet')

	ax3 = fig.add_subplot(223, projection=projection)
	zoutmap = BATHY_MAP( ztype='isol1000', ax=ax3 )
	PROJ_PLOT( zlon, zlat, zAWTmax_depth1, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, ax=ax3 )

	ax4 = fig.add_subplot(224, projection=projection)
	ztitle=' AW Max Temp depth from \n'+' WOA2009'
	zoutmap = BATHY_MAP( ztype='isol1000', ax=ax4 )
	PROJ_PLOT( zlon, zlat, zAWTmax_depthIS, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zplot_obs=0, ax=ax4 )
	plt.tight_layout()
	
	zfile_ext='_AWTmaxDepth_'
	plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zc_year)+'.png',dpi=300)

	if zncout:
		ds_out = xr.Dataset()
		
		# AWT field 
		#######################
		ds_out['AWTmax_mod']= (('y','x'), zAWTmax1.values.astype('float32')) 
		ds_out['AWTmax_mod'].attrs['long_name']='Model AWT max calculated using a salinity criteria value Sref=33.5'
		ds_out['AWTmax_mod'].attrs['units']='DegC'
                
		ds_out['AWTmax_init']= (('y','x'), zAWTmaxIS.values.astype('float32')) 
		ds_out['AWTmax_init'].attrs['long_name']='Initial state AWT max calculated using a salinity criteria value Sref=33.5'
		ds_out['AWTmax_init'].attrs['units']='DegC'
                
		ds_out['AWTmaxDepth_mod']= (('y','x'), zAWTmax_depth1.values.astype('float32')) 
		ds_out['AWTmaxDepth_mod'].attrs['long_name']='Model AWT max depth calculated using a salinity criteria value Sref=33.5'
		ds_out['AWTmaxDepth_mod'].attrs['units']='m'
                
		ds_out['AWTmaxDepth_init']= (('y','x'), zAWTmax_depthIS.values.astype('float32')) 
		ds_out['AWTmaxDepth_init'].attrs['long_name']='Initial state AWT max depth calculated using a salinity criteria value Sref=33.5'
		ds_out['AWTmaxDepth_init'].attrs['units']='m'
                
		ds_out['lat_mod']= (('y','x'), zlat.values.astype('float32')) 
		ds_out['lat_mod'].attrs['long_name']='Degrees north'
		ds_out['lat_mod'].attrs['units']='Deg'
		
		ds_out['lon_mod']= (('y','x'), zlon.values.astype('float32')) 
		ds_out['lon_mod'].attrs['long_name']='Degrees east'
		ds_out['lon_mod'].attrs['units']='Deg'
		
		ds_out['lat_mod']= (('y','x'), zlat.values.astype('float32')) 
		ds_out['lat_mod'].attrs['long_name']='Degrees north'
		ds_out['lat_mod'].attrs['units']='Deg'
		
		ds_out['lon_mod']= (('y','x'), zlon.values.astype('float32')) 
		ds_out['lon_mod'].attrs['long_name']='Degrees east'
		ds_out['lon_mod'].attrs['units']='Deg'

		ds_out = ds_out.set_coords(['lat_mod','lon_mod'])
		
		# Write the NetCDF file 
		ds_out.attrs['History'] = 'Diagnostics have been calculated using the Arctic monitoring tool '
		ds_out.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
		nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_AWTClim_'+'y'+str(zc_year)+'.nc'
		ds_out.to_netcdf(nc_f,engine='netcdf4')

	return

################################################################################################################################
def FWC_MAPSF( zlon, zlat, zSali, zSali_IS, zSSH, zCONF, zCASE, zc_year, ze3, ztmask, zncout, teos10 ) :
################################################################################################################################
	"""
	Function dedicated to plot the SSH and FWC 

	Input:
	    zlon     : longitude 2D  
	    zlat     : latitude 2D  
	    zSali    : Salinity 3D  
	    zSali_IS : Initial salinity 3D  
	    zSSH     : SSH  
	    zCONF    : configuration name 
	    zCASE    : experiment name associated to the configuration
	    zc_year  : current year 
  	    ze3      : 3D vertical scale factor 
  	    ztmask   : ocean/land mask at T-point 
	    zncout   : logical to outputs (or not) results into a netcdf file 
	    teos10   : logical to converts model temperature/salinity from CT/SA to Tpot/PS units
	
	Output:
	    None
	"""
	# ----------------------------------------------------------------------

	# FWC calculation over the year
	###########################################
	Sref = 34.80 

 	# Conversion from SA to PS
	if teos10 : 
		zSali = CONV_SA2PS( zSali, zlon.values, zlat.values ) 
		zSali_IS = CONV_SA2PS( zSali_IS, zlon.values, zlat.values ) 
	 
	print('				FWC calculation & plot ')

	# Freshwater content from model outputs 
	######################################################################################
	# Vales at depth are set to zero which leads to a problem 
	zSali = xr.where( npy.isnan(zSali), 9999., zSali )
	zSali = xr.where( zSali == 0., 9999., zSali )
	fwcmask = xr.where( zSali > Sref, 0., ze3 ) 
	FW4D = (Sref - zSali) / Sref * fwcmask
	# Sum over depth 
	fwc2D = FW4D.sum(dim="z")
	# Mask land area 
	fwc2D = xr.where( ztmask[0,:,:] < 1, npy.nan, fwc2D )

	# Freshwater content from initial state 
	######################################################################################
	# Vales at depth are set to zero which leads to a problem 
	zSali_IS = xr.where( npy.isnan(zSali_IS), 9999., zSali_IS )
	zSali_IS = xr.where( zSali_IS == 0., 9999., zSali_IS )
	fwcmask_init = xr.where( zSali_IS > Sref, 0., ze3 )
	FW4D_init = (Sref - zSali_IS) / Sref * fwcmask_init
	# Sum over depth 
	fwc2D_init = FW4D_init.sum(dim="z")
	# Mask land area 
	fwc2D_init = xr.where( ztmask[0,:,:] < 1, npy.nan, fwc2D_init )

	dbg=False
	if zCONF == 'CREG025.L75' :
		idbg=185 ; jdbg=515 ; kdbg=50	# CREG025.L75
	elif zCONF == 'CREG12.L75' :
		idbg=562 ; jdbg=1558 ; kdbg=50	 # CREG12.L75
	
	if dbg:
		print()
		print(' Print a specific point to debug within the Beaufort Gyre')
		print('   		zSali[0:kdbg,jdbg,idbg]: ',zSali[0:kdbg,jdbg,idbg])
		print('   		fwcmask[0:kdbg,jdbg,idbg]: ',fwcmask[0:kdbg,jdbg,idbg])
		print('   		ze3[0:kdbg,jdbg,idbg]: ',ze3[0:kdbg,jdbg,idbg])
		print('   		fwc2D[jdbg,idbg]: ',fwc2D[jdbg,idbg])


	# Read FWC inferred from Obs. 
	mean_FWCObs, lon2D_obs, lat2D_obs, obsper = FWC_OBS( t_year=int(zc_year) )

	# Read DOT observations 
	obs_ssh, lon_obs, lat_obs, obs_ssh_per = SSH_OBS( t_year=int(zc_year) )

	plt.clf()
	fig, axes = plt.subplots( 2, 3, figsize=(11.7,8.3),  subplot_kw={'projection': ccrs.NorthPolarStereo(central_longitude=-60, true_scale_latitude=65)} )
	axes = axes.flatten()

	# Plot the FWC map mean over the year
	#####################################
	
	zMyvar='FWC'
	seas=''
	contours, limits, myticks, ztitle, zfile_ext, my_cblab, my_cmap, m_alpha = SET_ARC_CNT( zCASE, zc_year, seas, zMyvar )
	zoutmap = BATHY_MAP( ztype='isol1000', ax=axes[0] )
	PROJ_PLOT( zlon.values, zlat.values, fwc2D, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, ax=axes[0] )

	maxval=npy.round(npy.nanmax(mean_FWCObs), decimals=2)
	fig.text(0.45,0.80,'Max: '+str(maxval)+' m',fontsize=9,color='r')
	ztitle=' FWC (m) from \n'+' BG Obs Sys. (Proshutinsky et al. GRL2018) \n '+ ' year ' + obsper
	zoutmap = BATHY_MAP( ztype='isol1000', ax=axes[1] )
	PROJ_PLOT( lon2D_obs, lat2D_obs, mean_FWCObs, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, ax=axes[1], zplot_obs=1 )

	ztitle=' FWC (m) from \n'+' Init State '
	zoutmap = BATHY_MAP( ztype='isol1000', ax=axes[2] )
	PROJ_PLOT( zlon, zlat, fwc2D_init, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, ax=axes[2] )
	
	zMyvar='ssh'
	seas=''
	contours, limits, myticks, ztitle, zfile_ext, my_cblab, my_cmap, m_alpha = SET_ARC_CNT( zCASE, zc_year, seas, zMyvar )
	zoutmap = BATHY_MAP( ztype='isol1000', ax=axes[3] )
	PROJ_PLOT( zlon, zlat, zSSH*m_alpha, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, ax=axes[3] )

	zMyvar='ssh'
	seas=''
	contours, limits, myticks, ztitle, zfile_ext, my_cblab, my_cmap, m_alpha = SET_ARC_CNT( zCASE, zc_year, seas, zMyvar, zplot_obs=1 )
	zoutmap = BATHY_MAP( ztype='isol1000', ax=axes[4] )
	ztitle = ztitle+obs_ssh_per
	PROJ_PLOT( lon_obs, lat_obs, obs_ssh*m_alpha, contours, limits, zmy_ticks=myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, ax=axes[4], zplot_obs=1 )

	fig.delaxes(axes[5])
	plt.tight_layout()

	zfile_ext='_FWCSSHClim_'
	plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zc_year)+'.png',dpi=300)

	if zncout:
		ds_out = xr.Dataset()
		
		# FWC field 
		#######################
		fwc2D = xr.where( ztmask[0,:,:] < 1, npy.nan, fwc2D )
		ds_out['fwc_mod']= (('y','x'), fwc2D.values.astype('float32')) 
		ds_out['fwc_mod'].attrs['long_name']='Model FWC calculated using a salinity reference value Sref='+str(Sref)
		ds_out['fwc_mod'].attrs['units']='m'
                
		fwc2D_init = xr.where( ztmask[0,:,:] < 1, npy.nan, fwc2D_init )
		ds_out['fwc_init']= (('y','x'), fwc2D_init.values.astype('float32')) 
		ds_out['fwc_init'].attrs['long_name']='Model Initial state FWC calculated using a salinity reference value Sref='+str(Sref)
		ds_out['fwc_init'].attrs['units']='m'
                
		ds_out['fwc_obs']= (('yobs','xobs'), mean_FWCObs.values.astype('float32')) 
		ds_out['fwc_obs'].attrs['long_name']='Obs. FWC calculated using a salinity reference value Sref=34.8 '+ \
						     ' from Proshutinsky et al. (GRL2018). Considered period: '+obsper
		ds_out['fwc_obs'].attrs['units']='m'

		ds_out['lat_obs']= (('yobs','xobs'), lat2D_obs.astype('float32')) 
		ds_out['lat_obs'].attrs['long_name']='Degrees north'
		ds_out['lat_obs'].attrs['units']='Deg'
		
		ds_out['lon_obs']= (('yobs','xobs'), lon2D_obs.astype('float32')) 
		ds_out['lon_obs'].attrs['long_name']='Degrees east'
		ds_out['lon_obs'].attrs['units']='Deg'
		
		ds_out['lat_mod']= (('y','x'), zlat.values.astype('float32')) 
		ds_out['lat_mod'].attrs['long_name']='Degrees north'
		ds_out['lat_mod'].attrs['units']='Deg'
		
		ds_out['lon_mod']= (('y','x'), zlon.values.astype('float32')) 
		ds_out['lon_mod'].attrs['long_name']='Degrees east'
		ds_out['lon_mod'].attrs['units']='Deg'

		ds_out = ds_out.set_coords(['lat_obs','lon_obs','lat_mod','lon_mod'])

		# Write the NetCDF file 
		ds_out.attrs['History'] = 'Diagnostics have been calculated using the Arctic monitoring tool '
		ds_out.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
		nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_FWCClim_'+'y'+str(zc_year)+'.nc'
		ds_out.to_netcdf(nc_f,engine='netcdf4')

		ds_out = xr.Dataset()

		# SSH field 
		#######################
		zSSH = xr.where( ztmask[0,:,:] < 1, npy.nan, zSSH )
		ds_out['ssh_mod']= (('y','x'), zSSH.values.astype('float32')) 
		ds_out['ssh_mod'].attrs['long_name']='Model SSH '
		ds_out['ssh_mod'].attrs['units']='m'
                
		ds_out['ssh_obs']= (('yobs','xobs'), obs_ssh.values.astype('float32')) 
		ds_out['ssh_obs'].attrs['long_name']='Observed DOT from Armitage et al. 2017 for the period: '+obs_ssh_per
		ds_out['ssh_obs'].attrs['units']='m'
                
		ds_out['lat_obs']= (('yobs','xobs'), obs_ssh.lat.values.astype('float32')) 
		ds_out['lat_obs'].attrs['long_name']='Degrees north'
		ds_out['lat_obs'].attrs['units']='Deg'
		
		ds_out['lon_obs']= (('yobs','xobs'), obs_ssh.lon.values.astype('float32')) 
		ds_out['lon_obs'].attrs['long_name']='Degrees east'
		ds_out['lon_obs'].attrs['units']='Deg'
		
		ds_out['lat_mod']= (('y','x'), zlat.values.astype('float32')) 
		ds_out['lat_mod'].attrs['long_name']='Degrees north'
		ds_out['lat_mod'].attrs['units']='Deg'
		
		ds_out['lon_mod']= (('y','x'), zlon.values.astype('float32')) 
		ds_out['lon_mod'].attrs['long_name']='Degrees east'
		ds_out['lon_mod'].attrs['units']='Deg'

		ds_out = ds_out.set_coords(['lat_obs','lon_obs','lat_mod','lon_mod'])

		# Write the NetCDF file 
		ds_out.attrs['History'] = 'Diagnostics have been calculated using the Arctic monitoring tool '
		ds_out.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
		nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_SSHClim_'+'y'+str(zc_year)+'.nc'
		ds_out.to_netcdf(nc_f,engine='netcdf4')

	return

################################################################################################################################
def CREG_MSK( zCONF, zCASE ) :
################################################################################################################################
        """
        Function dedicated to read the ocean/land mask at T-point as the geographical coordinates
        Filename to read in the directory './'+zCONF+'/GRID/' is zCONF-zCASE_mask.nc
        
        Input:
            zCONF : the configuration name 
            zCASE : the experiment name associated to the configuration
        
        Output:
            mask : ocean/land at T-point 
            lon  : 2D longitude  
            lat  : 2D latitude  
        """
        # ----------------------------------------------------------------------

        # Read the CREG025.L75 mask 
        locpath='./'+zCONF+'/GRID/'
        locfile=zCONF+'-'+zCASE+'_mask.nc'
        if chkfile(locpath+locfile) :
                ds_msk = xr.open_dataset( locpath+locfile )[['glamt','gphit','tmask']]
                lon  = ds_msk['glamt'].squeeze()
                lat  = ds_msk['gphit'].squeeze()
                mask = ds_msk['tmask'].squeeze()

        return mask, lon, lat

################################################################################################################################
def CREG_INIT( zCONF, zCASE ) :
################################################################################################################################
	"""
	Function dedicated to read the initial state file used for the concerned experiment 
	Filename to read in the directory zCONF+'/'+zCONF+'-'+zCASE+'-MEAN/ is zCONF-zCASE_init_gridT.nc
	
	Input:
	    zCONF : the configuration name 
	    zCASE : the experiment name associated to the configuration
	
	Output:
	    ds_TSinit : Initial temperature/salinity 3D  
	    zTemp_IS  : Initial temperature 3D  
	    zSali_IS  : Initial salinity 3D  
	"""
	# ----------------------------------------------------------------------
	# Read initial state to compare with
	print('                      Read initial state  ')
	locpath=zCONF+'/'+zCONF+'-'+zCASE+'-MEAN/'
	locfile=zCONF+'-'+zCASE+'_init_gridT.nc'
	if chkfile(locpath+locfile) : 
		ds_TSinit = xr.open_dataset(locpath+locfile, engine="netcdf4")[['votemper','vosaline']]
		ds_TSinit = ds_TSinit.rename({'nav_lev':'z'})
		zTemp_IS = ds_TSinit['votemper'].squeeze()
		zSali_IS = ds_TSinit['vosaline'].squeeze()
	
	return ds_TSinit, zTemp_IS, zSali_IS

################################################################################################################################
def WOA09_INIT( ) :
################################################################################################################################
	"""
	Function dedicated to read the WOA09 climatological temperature/salinity (TEOS10)
	Filename to read in the directory './DATA/':  
		- woa09_SalAbs_monthly_1deg_SA_CMA_drowned_Ex_L75_SM5.nc 
		- woa09_ConTem_monthly_1deg_CT_CMA_drowned_Ex_L75_SM5.nc 
	
	Input: 
	      None
	
	Output:
	    zTS_WOA   : Initial temperature/salinity 3D  
	    zTemp_WOA : Initial temperature 3D  
	    zSali_WOA : Initial salinity 3D  
	"""
	# ----------------------------------------------------------------------
	# Read initial state to compare with
	print('                      Read WOA09 climatology  ')
	locpath='./DATA/'
	locfileT='woa09_ConTem_monthly_1deg_CT_CMA_drowned_Ex_L75_SM5.nc'
	locfileS='woa09_SalAbs_monthly_1deg_SA_CMA_drowned_Ex_L75_SM5.nc'
	if chkfile(locpath+locfileS) and chkfile(locpath+locfileT) : 
		ds_T = xr.open_dataset(locpath+locfileT, engine="netcdf4")
		ds_S = xr.open_dataset(locpath+locfileS, engine="netcdf4")
		zTS_WOA = xr.Dataset()
		zTS_WOA['votemper'] = (('time_counter','z','lon0','lat0'),ds_T.CT.values )
		zTS_WOA['vosaline'] = (('time_counter','z','lon0','lat0'),ds_S.SA.values )
		zTS_WOA['lon'] = (('lon0','lat0'),ds_T.lon.values )
		zTS_WOA['lat'] = (('lon0','lat0'),ds_T.lat.values )

		zTemp_WOA = ds_T.CT.values
		zSali_WOA = ds_S.SA.values
	
	return zTS_WOA, zTemp_WOA, zSali_WOA

################################################################################################################################
def BATHY_MAP( ztype='isol1000', zarea='arctic', ax=None ) :
################################################################################################################################
	"""
	Function dedicated to plot bathymetry isolines  
	Filename to read is Bathymetry.nc
	
	Input:
	    ztype : (optionnal) to set bathymetry isolines contours to plot (default = 'isol1000')
	    zarea : (optionnal) specify the plot area (default = 'arctic')
	    ax    : (optionnal) considered plot axes 
	
	Output:
	    None
	"""
	# ----------------------------------------------------------------------

	locpath='./'
	locfile='Bathymetry.nc'
	if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
		ds_bat = xr.open_dataset(locpath+locfile)
		lon = ds_bat['nav_lon'].squeeze()
		lat = ds_bat['nav_lat'].squeeze()
		zBathy = ds_bat['bathy_meter'].squeeze()
	
	spval = 0.
	zBathy = xr.where( zBathy <= spval, npy.nan, zBathy )
	
	if ztype == 'isol1000' :
		vmin=1000. ; vmax=2000. 
		contours=[1000.]
		limits=[vmin,vmax]  
		myticks=[1000.]
	elif ztype == 'isol1500' :
		vmin=1500. ; vmax=2000. 
		contours=[1500.]
		limits=[vmin,vmax]  
		myticks=[1500.]
	elif ztype == 'isomonarc' :
		vmin=500. ; vmax=4000. 
		contours=[500.,2000.,4000.]
		limits=[vmin,vmax]  
		myticks=[500.,2000.,4000.]
	elif ztype == 'isol500' :
		vmin=500. ; vmax=500. 
		contours=[500.]
		limits=[vmin,vmax]  
		myticks=[500.]
	else:
		vmin=0. ; vmax=8000. 
		contours=[100.,500.,1000.,2000.,3000.,3500.,4000.]
		limits=[vmin,vmax] 
		myticks=[100.,500.,1000.,2000.,3000.,3500.,4000.]
	
	#
	plt.rcParams['text.usetex']=False
	plt.rcParams['font.family']='serif'
	plt.rcParams['axes.unicode_minus'] = False
	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	zcolorbat='grey'  ;  zalpha=0.4

	if zarea == 'arctic': # Focus on Arctic
		#m = Basemap(projection='npstere',boundinglat=65,lon_0=-60, resolution='i')
		ax.set_aspect('equal')  # Important pour éviter les distorsions
		ax.set_extent([-180, 180, 65, 90], crs=ccrs.PlateCarree())  # Optionnel : limite la vue
	elif zarea == 'labsea': # Focus on Labrador Sea
		m = Basemap(width=1400000,height=1600000,lat_1=50.,lat_2=65,lon_0=-50,lat_0=59.5,projection='aea',resolution='i')
	elif zarea == 'GulfS': # Focus on Gulf Stream area
		my_area = {'lonmin':-80., 'lonmax':-40.,'latmin':30.,'latmax':50.}
		m = Basemap(projection='cyl',llcrnrlat=my_area['latmin'],urcrnrlat=my_area['latmax'],\
					     llcrnrlon=my_area['lonmin'],urcrnrlon=my_area['lonmax'],resolution='i')
	elif zarea == 'irmsea': # Focus on Irminger Sea
		m = Basemap(width=1800000,height=1600000,lat_1=50.,lat_2=65,lon_0=-30,lat_0=59.5,projection='aea',resolution='i')
		############################################################################################################
		bx_ISB={'name':'ISB'  ,'lon_min':-37.,'lon_max':-37.,'lat_min':61.,'lat_max':61.}
		All_box=[bx_ISB]
		for box in All_box:
			lats = [box['lat_min'],box['lat_max']]
			lons = [box['lon_min'],box['lon_max']]
			x,y = m(lons,lats)
			m.scatter(x,y,1,marker='o', color='r')
		############################################################################################################
	elif zarea == 'ginsea': # Focus on GIN Seas
		m = Basemap(width=1400000,height=1600000,lat_1=50.,lat_2=65,lon_0=0,lat_0=74.,projection='aea',resolution='i')
	############################################################################################################
	elif zarea == 'cassis_BGZoom' :
		#m = Basemap(llcrnrlon=-180,llcrnrlat=66,urcrnrlon=-80,urcrnrlat=80, resolution='i',\
		#            projection='cass',lon_0=-140,lat_0=60)    
		ax.set_extent([-180, -100, 66, 80], crs=ccrs.PlateCarree())

		ax.add_feature(cartopy.feature.LAND, facecolor='dimgray')
		gridlines = ax.gridlines(draw_labels=True, linestyle=':', linewidth=0.4, alpha=0.7, xlocs=npy.arange(-180, -80, 10), y_inline=True, rotate_labels=False )
		gridlines.xlabel_style={'fontsize': 4}
		gridlines.ylabel_style={'fontsize': 4}
		gridlines.top_labels=False
		gridlines.left_labels=True
		gridlines.right_labels=False
	############################################################################################################
	elif zarea == 'cassis_BGZoom_HR' :
		#m = Basemap(llcrnrlon=-80,llcrnrlat=80,urcrnrlon=-180,urcrnrlat=60, resolution='i',\
		#            projection='cass',lon_0=0,lat_0=80)    
		ax.set_extent([-180, -80, 66, 80], crs=ccrs.PlateCarree())
	############################################################################################################
	else: # Focus on North Atlantic sector
		m = Basemap(width=6100000,height=5000000,lat_1=30.,lat_2=70,lon_0=-45,lat_0=45,projection='aea',resolution='i')
		zcolorbat='grey'   ;  zalpha=0.7

	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	pal = plt.get_cmap('binary')

	# contour (optional)
	if zarea == 'cassis_BGZoom' or zarea == 'cassis_BGZoom_HR' or zarea == 'arctic' :
		CS2 = ax.contour( lon, lat, zBathy.values, linewidths=0.5, levels=contours, colors=zcolorbat, alpha=zalpha, transform=ccrs.PlateCarree() )
	else :
		X,Y = m(lon.values,lat.values)
		CS2 = m.contour( X, Y, zBathy.values, linewidths=0.5,levels=contours, colors=zcolorbat, alpha=zalpha )
	plt.clabel(CS2, CS2.levels, inline=True, fmt='%.0f', fontsize=3)

	return

################################################################################################################################
def PROJ_PLOT( zlon, zlat, tab, contours, limits, name=None, zmy_ticks=None, zmy_cblab=None, zmy_cmap=None, zvar=None, zarea='arctic', data_ref=False, ax=None, zplot_obs=0, zmy_year=None ) :
################################################################################################################################
	"""
	Function dedicated to prepare the map projection depending the variable and geographical location to plot 
	
	Input:
	    zlon      : longitude 2D  
	    zlat      : latitude 2D  
	    tab       : the variable to plot 
	    contours  : contours to plot 
	    limits    : minimum & maximum values of the variable 
	    name      : (optionnal) the title name of the plot 
	    zmy_ticks : (optionnal) the colorbar ticks 
	    zmy_cblab : (optionnal) the colorbar labels 
	    zmy_cmap  : (optionnal) the colormap to use 
	    zvar      : (optionnal) the variable name to plot
	    zarea     : (optionnal) specify the plot area (default = 'arctic')
	    data_ref  : (optionnal) logical to specify if data is a reference one or not (default = False)
	    ax        : (optionnal) considered plot axes
	    zplot_obs : (optionnal) to speficy if the data to plot is from the model or Obs. (default = 0)
	
	Output:
	    None
	"""
	# ----------------------------------------------------------------------
	#
	plt.rcParams['text.usetex']=False
	plt.rcParams['font.family']='serif'
	plt.rcParams['axes.unicode_minus'] = False
	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	if zvar == 'sivolu' or zvar == 'siconc'  or zvar == 'ssh' or zvar == 'FWC' or zvar == 'voeke' or zvar == 'mldr10_1' :
		zfontsize=11.
	elif zvar == 'sobarstf' or zvar == 'topos' : 
		zfontsize=8.
	else:
		zfontsize=6.
	
	############################################################################################################
	if zarea == 'arctic': # Focus on Arctic basin
		ax.set_aspect('equal')  
		ax.set_extent([-180, 180, 65, 90], crs=ccrs.PlateCarree())
		
		ax.add_feature(cartopy.feature.LAND, facecolor='dimgray')
		gridlines = ax.gridlines(draw_labels=True, linestyle=':', linewidth=0.4, alpha=0.7, xlocs=npy.arange(-180, 181, 20), y_inline=False, rotate_labels=False )
		gridlines.xlabel_style={'fontsize': zfontsize}
		gridlines.ylabel_style={'fontsize': zfontsize}
		gridlines.top_labels=False
		gridlines.left_labels=False
		gridlines.right_labels=False

	############################################################################################################
	elif zarea == 'labsea': # Focus on Gulf Stream area
		m = Basemap(width=1400000,height=1600000,lat_1=50.,lat_2=65,lon_0=-50,lat_0=59.5,projection='aea',resolution='i')
		bx_LABK1={'name':'K1'  ,'lon_min':-52.4,'lon_max':-52.4,'lat_min':56.3,'lat_max':56.3}
		All_box=[bx_LABK1]
		for box in All_box:
			lats = [box['lat_min'],box['lat_max']]
			lons = [box['lon_min'],box['lon_max']]
			x,y = m(lons,lats)
			m.scatter(x,y,1,marker='o', color='r')
	############################################################################################################
	elif zarea == 'GulfS': # Focus on Gulf Stream area
		my_area = {'lonmin':-80., 'lonmax':-40.,'latmin':30.,'latmax':50.}
		m = Basemap(projection='cyl',llcrnrlat=my_area['latmin'],urcrnrlat=my_area['latmax'],\
					     llcrnrlon=my_area['lonmin'],urcrnrlon=my_area['lonmax'],resolution='i')
	############################################################################################################
	elif zarea == 'irmsea': # Focus on Irminger Sea
		m = Basemap(width=1800000,height=1600000,lat_1=50.,lat_2=65,lon_0=-30,lat_0=59.5,projection='aea',resolution='i')
		bx_ISB={'name':'ISB'  ,'lon_min':-37.,'lon_max':-37.,'lat_min':61.,'lat_max':61.}
		All_box=[bx_ISB]
		for box in All_box:
			lats = [box['lat_min'],box['lat_max']]
			lons = [box['lon_min'],box['lon_max']]
			x,y = m(lons,lats)
			m.scatter(x,y,1,marker='o', color='r')
	############################################################################################################
	elif zarea == 'cassis_BGZoom' :
		#m = Basemap(llcrnrlon=-180,llcrnrlat=66,urcrnrlon=-80,urcrnrlat=80, resolution='i',\
		#            projection='cass',lon_0=-140,lat_0=60)    
		#ax.set_extent([-180, -80, 66, 80], crs=ccrs.PlateCarree())
		#ax.set_extent([-180, -120, 66, 80], crs=ccrs.PlateCarree())
		ax.set_extent([-180, -80, 66, 85], crs=ccrs.PlateCarree())

		ax.add_feature(cartopy.feature.LAND, facecolor='dimgray')
		gridlines = ax.gridlines(draw_labels=False, linestyle=':', linewidth=0.4, alpha=0.7, xlocs=npy.arange(-180, -80, 10), y_inline=True, rotate_labels=False )
		gridlines.xlabel_style={'fontsize': 4}
		gridlines.ylabel_style={'fontsize': 4}
		gridlines.top_labels=False
		gridlines.left_labels=True
		gridlines.right_labels=False
	############################################################################################################
	elif zarea == 'cassis_BGZoom_HR' :
		m = Basemap(llcrnrlon=-80,llcrnrlat=80,urcrnrlon=-180,urcrnrlat=60, resolution='i',\
		            projection='cass',lon_0=0,lat_0=80)    
	############################################################################################################
	elif zarea == 'ginsea': # Focus on GIN Seas
		m = Basemap(width=1400000,height=1600000,lat_1=50.,lat_2=65,lon_0=0,lat_0=74.,projection='aea',resolution='i')
	############################################################################################################
	elif zarea == 'natl': # Focus on North Atlantic sector
		 m = Basemap(width=6100000,height=5000000,lat_1=30.,lat_2=70,lon_0=-45,lat_0=45,projection='aea',resolution='i')
	
	# Need to kepp the following lines for compatibility with Basemap
	if zarea == 'GulfS' or zarea == 'labsea' or zarea == 'irmsea' or zarea == 'ginsea' or zarea == 'natl':
		m.drawparallels(npy.arange(-90.,91.,2.),labels=[True,False,False,False], size=zfontsize, linewidth=0.3, color='grey',alpha=0.70 )
		m.drawmeridians(npy.arange(-180.,181.,5.),labels=[False,False,False,True], size=zfontsize, latmax=90.,linewidth=0.3, color='grey',alpha=0.70 )
		m.fillcontinents(color='grey',lake_color='white')
	
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	
	if zmy_cmap != None :
		pal = zmy_cmap
	else:
		pal = plt.get_cmap('coolwarm')
	
	if zarea == 'GulfS' and zvar == 'votemper' :
		if data_ref :
			zlinewidths=1.1   ; zcolor='g'
		else:	
			zlinewidths=0.8   ; zcolor='r'
		X,Y = m(zlon,zlat)
		C = m.contour( X, Y, tab, linewidths=zlinewidths, levels=[17.], colors=zcolor )

	elif zarea == 'labsea' or zarea == 'irmsea' or zarea == 'ginsea' or zarea == 'natl' :
		X,Y = m(zlon,zlat)
		C = m.contourf( X,Y,tab,contours,cmap=pal,norm=norm,extend='both' )

	else:
		if zplot_obs == 0 : 
			C = tab.plot.pcolormesh( ax=ax, levels=contours, cmap=pal, transform = ccrs.PlateCarree(), add_colorbar=False, add_labels=True )
		else :
			if zvar == 'ssh' or zvar == 'mldr10_1' or zvar == 'siconc' :
				C = ax.pcolormesh( zlon, zlat, tab, vmin=limits[0], vmax=limits[1], cmap=zmy_cmap, transform = ccrs.PlateCarree() )
			elif zvar == 'voeke' :
				if zmy_year >= 2015 : # Plot Cryosat data 
					C = ax.pcolormesh( zlon, zlat, tab, vmin=limits[0], vmax=limits[1], cmap=zmy_cmap, transform = ccrs.PlateCarree() )
				else :
					C = ax.contourf( zlon, zlat, tab, levels=contours, cmap=pal, norm=norm, extend='both', transform = ccrs.PlateCarree() )
			else :
				C = ax.contourf( zlon, zlat, tab, levels=contours, cmap=pal, norm=norm, extend='both', transform = ccrs.PlateCarree() )

		############################################################################################################
		############################################################################################################
		moorplot=1
		if moorplot == 1 :
				bx_ARCB={'name':'B'  ,'lon_min':-150.,'lon_max':-150.,'lat_min':78.,'lat_max':78.}
				bx_ARCM={'name':'M1' ,'lon_min': 125.,'lon_max': 125.,'lat_min':78.,'lat_max':78.}
				bx_EURA={'name':'EUR','lon_min':  60.,'lon_max':  60.,'lat_min':85.,'lat_max':85.}

				All_box=[bx_ARCB,bx_EURA]
				for box in All_box:
					lats = [box['lat_min'],box['lat_max']]
					lons = [box['lon_min'],box['lon_max']]
					ax.scatter(lons,lats,1,marker='o', color='r', transform=ccrs.PlateCarree())
		############################################################################################################
		############################################################################################################
	
		# colorbar	
		if zmy_ticks is None:
			cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
		else:
			if zvar == 'votemper' or zvar == 'vosaline' or zvar == 'sivolu' :
				cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8,drawedges=True)
			elif zvar == 'sobarstf' or zvar == 'topos' :
				cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.6,drawedges=True)
			elif zvar == 'MLTSS' :
				cbar = plt.colorbar(C,ticks=zmy_ticks,format='%.0f',orientation='vertical',shrink=0.8,drawedges=True)
			elif zvar == 'voeke' :
				if zplot_obs == 1 :
					if zmy_year >= 2015 : # Plot Cryosat data 
						cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8,extend='both')
					else :
						cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8,drawedges=True)
				else :
					cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8,drawedges=True)
			elif zvar == 'mldr10_1' :
				cbar = plt.colorbar(C,ticks=zmy_ticks,format='%.0f',orientation='vertical',shrink=0.6,extend='both')
			elif zvar == 'ssh' or zvar == 'FWC' :
				cbar = plt.colorbar(C,ticks=zmy_ticks,format='%.0f',orientation='vertical',shrink=0.6,extend='both')
			elif zvar == 'siconc' :
				cbar = plt.colorbar(C,format='%.0f',orientation='vertical',shrink=0.8,extend='both')
			else:
				cbar = plt.colorbar(C,format='%.0f',orientation='vertical',shrink=0.8,drawedges=True)

			cbar.set_label(zmy_cblab,fontsize=zfontsize)
			cl = plt.getp(cbar.ax, 'ymajorticklabels')
			plt.setp(cl, fontsize=zfontsize)
			#if ztickslabels != None and zvar == 'voeke' : 
			#	zticks = npy.linspace(1e-6, 1e-2, 5)
			#	cbar.set_ticks(10**zticks)
			#	cbar.ax.set_yticklabels(ztickslabels)
	
	ax.set_title(name,fontsize=zfontsize)
	
	return

################################################################################################################################
def DO_MAPS( zlon, zlat, zCONF, zCASE, zVar, zVarname, zc_year, slev=None, seas='', plot_obs=0, ano=0, ax=None, over_title=None ) :
################################################################################################################################
	"""
	Function dedicated to prepare the data projection on a map
	
	Input:
	    zlon     : longitude 2D  
	    zlat     : latitude 2D  
	    zCONF    : configuration name 
	    zCASE    : experiment name associated to the configuration
	    zVar     : variable to plot 
	    zVarname : variable name to plot 
	    zc_year  : the current year 
	    slev     : (optionnal) vertical index level to plot 
	    seas     : (optionnal) specify the month to plot either 'm03' or 'm09' (default = '')
	    plot_obs : (optionnal) to speficy if the data to plot is from the model or Obs. (default = 0)
	    ano      : (optionnal) specify if the plot is an anomaly or the variable itself (default = 0)
	    ax       : (optionnal) considered plot axes
	
	Output:
	    m : the plot projection done 
	"""
	# ----------------------------------------------------------------------

	# Do the plot 
	print() 
	print('                    plot '+zVarname+' field')
	print() 
	
	m_alpha = 1.
	contours, limits, myticks, ztitle, zfile_ext, my_cblab, my_cmap, m_alpha = SET_ARC_CNT( zCASE, zc_year, seas, zVarname, zslev=slev, zplot_obs=plot_obs, zdiff=ano )
	if over_title != None : ztitle = over_title 
	BATHY_MAP( ztype='isol1000', ax=ax )
	m = PROJ_PLOT( zlon, zlat, zVar[:,:]*m_alpha, contours, limits, name=ztitle, zmy_ticks=myticks, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zVarname, zplot_obs=plot_obs, ax=ax, zmy_year=zc_year )

	return m 

################################################################################################################################
def BFG_COMPUTE( lon, lat, ssh_raw, depth, var_type, increment, grid_area, rm_landbarrier=0 ) :
################################################################################################################################
	"""
	Function dedicated to computes the largest closed contour in the Western Arctic basin

	History: This code has been developed by Heather Regan and slightly adapted to be included into the MONARC
		 See Regan et al. JPO2020 ; https://doi.org/10.1175/JPO-D-19-0234.1
	
	Input:
	    zlon           : longitude 2D  
	    zlat           : latitude 2D  
	    ssh_raw        : the SSH we're examining
	    depth          : the bathymetry (depth > 0) to check off-shelf regions
	    var_type       : usually set to "SSH": identifies the variable type (has also been used for MSL in past)
	    increment      : the increment with which to iterate out from the maximum. Usually start at 10cm.
	                     Higher resolution needs a smaller increment than lower resolution, because field varies 
	                     more and so larger increment may miss small features. But smaller increment takes longer
	    rm_landbarrier : use (or not) coastline as a valid edge of contour (e.g. set to 1 if using MSL, as atmospheric variable)
	
	Output:
	    mask_full  : the identified closed contour
	    BG_max_val : Max. SSH value 
	    BGcalcmin  : Min. SSH value
	    lat        : latitude
	    lon        : longitude
	    BG_area    : Surface area of the closed contour 
	"""
	# ----------------------------------------------------------------------

	## Step 1: set up coastline to determine when the contour is no longer closed
	# ---------------------------------------------------------------------------
	if rm_landbarrier==0: 
	  land_arr = npy.nan*npy.ones([npy.shape(lon)[0],npy.shape(lon)[1]]); land_arr[npy.isnan(ssh_raw)] = 1;
	  #First, reduce the array so that it takes less time
	  lonmask = npy.nan*npy.ones(lon.shape);
	  lonmask[lon<-80] = 1; lonmask[lon>140] = 1; lonmask[lat<68] = npy.nan;
	  lonmask[(lon<20) & (lon > -130) & (lat < 70)] = npy.nan;        
	  lonmask[(lon<20) & (lon > -120) & (lat < 75.5)] = npy.nan;
	  lonmask[(lon<20) & (lon > -110) & (lat < 73)] = npy.nan;
	  lonmask[(lon<20) & (lon > -100) & (lat < 80)] = npy.nan;
	  lonmask[(lon<20) & (lon > -90) & (lat < 80) & (depth < 1000)] = npy.nan;
	  ## extra
	  lonmask[lat>80]  = npy.nan
	  
	  ssh_full = ssh_raw*lonmask;
	  landmask = npy.zeros(lon.shape); landmask[ssh_raw==0] = 1; landmask[npy.isnan(ssh_raw)] = 1;
	
	else:
	  land_arr = npy.zeros([npy.shape(lon)[0],npy.shape(lon)[1]]); land_arr[npy.isnan(ssh_raw)] = 1;
	  lonmask = npy.ones(lon.shape)
	  ssh_full = ssh_raw*lonmask
	
	## Step 2: Identify off-shelf maximum    
	# ---------------------------------------------------------------------------
	# This ensures that the maximum nonzero value off the shelf is found. Otherwise can get a high maxima near the coast
	# Artificially force by depth field if rm_landbarrier isn't there 
	shelfmask = npy.ones(ssh_full.shape);
	if var_type == 'MSL':
	  shelfmask[depth<0] = npy.nan; shelfmask[depth>=0] = 1;
	else:
	  shelfmask[depth<3000] = npy.nan; shelfmask[depth>=3000] = 1; ## Here, make sure that depth array is in form of depths > 0!
	
	masked_shelf = ssh_full*shelfmask;
	masked_shelf[masked_shelf==0] = npy.nan;
	masked_shelf[npy.isnan(masked_shelf)] = npy.nan;
	maxarr = npy.nanmax(masked_shelf);
	print('			max '+str(maxarr))
	maxarr_whole = maxarr.copy();
	
	if var_type == 'MSL':
	  inc_min = increment/10000;
	else:
	  inc_min = increment/1000; ## Was previously 100, for monthly values. Changed to 1000 for yearly. RERUN FOR MONTHLY
	
	#Here the array is reduced to a more manageable size
	ssh_xsum = npy.nansum(ssh_full*lonmask,axis=1);
	ssh_ysum = npy.nansum(ssh_full*lonmask,axis=0);
	xfind = npy.nonzero(ssh_xsum); x1 = npy.nanmin(xfind); x2 = npy.nanmax(xfind);
	yfind = npy.nonzero(ssh_ysum); y1 = npy.nanmin(yfind); y2 = npy.nanmax(yfind);
	
	ssh = ssh_full[x1:x2,y1:y2].copy(); dsmall = depth[x1:x2,y1:y2].copy();
	ssh_3000 = ssh.copy(); 
	if var_type != 'MSL':
	  ssh_3000[dsmall<3000] = npy.nan
	
	mask_full = npy.zeros(ssh_full.shape);
	maskarr = npy.zeros(ssh.shape);
	land_arr = npy.zeros(ssh.shape); land_arr[npy.isnan(ssh)] = 1;
	 
	if var_type=='FW': 
	  land_arr[ssh==0] = 1; 
	  land_arr[npy.isnan(ssh)] = 1; 
	
	## Step 3: now loop over increments to find largest contour
	# ---------------------------------------------------------------------------
	# For each new increment, check that all cells in this new contour are a) not adjacent to land, and b) not higher than previous maximum
	  
	## LOOP #########################################
	#################################################
	if abs(maxarr) == 0:
	    print('		BG not found')
	    all_met = 1;
	else:
	    all_met = 0;
	    for x in range (0,ssh.shape[0]):
	        for y in range (0,ssh.shape[1]):
	            if ssh_3000[x,y] >= maxarr_whole:
	                maskarr[x,y] = 1;
	
	size_of_old_mask = 0;
	size_of_new_mask = 0;
	while_loop = 0;
	cond = 1;
	  
	## We have the maximum value. Basically store checkarr coordinates and loop over it
	while all_met == 0:
	    maskarr_new = maskarr.copy();
	    reloop = 1;
	    while_loop = while_loop + 1;
	
	    ###################
	    #Here the new maximum contour is found
	    #Define a new edge array based on mask
	    checkarr = FND_CEDGE (maskarr_new);
	    #Generate list of coordinates of new edge
	    [cx,cy] = npy.where(checkarr==1);
	    near_ocean = 1;
	
	    length_of_carr = cx.shape[0];
	    no_in_mask = npy.sum(maskarr_new+checkarr);
	    inc_land_mask = maskarr_new+checkarr; inc_land_mask[land_arr==cond] = 0;
	    if no_in_mask == npy.sum(inc_land_mask):
	        looping = 1; end_of_loop = length_of_carr;
	    else:
	        looping = 1; end_of_loop = 1;
	        near_ocean = 0; reloop = 0;
	
	    while looping < end_of_loop:
	        thisx = cx[looping];
	        thisy = cy[looping];
	        #First check that it's not land. If it is, exit the loop        
	        if land_arr[thisx,thisy] != cond:
	            #Need to check that this is next to the mask containing the maximum
	            if ssh[thisx,thisy]>= maxarr:
	                maskarr_new[thisx,thisy] = 1;
	                #Now loop over surrounding cells
	                for yval in range (-1,2):
	                    for xval in range (-1,2):
	                        nexty = cy[looping]+yval; nextx = cx[looping]+xval;
	                        #Check edges of domain
	                        if nexty>=1: 
	                            if nexty<ssh.shape[1]: 
	                                if nextx>=1: 
	                                    if nextx<ssh.shape[0]:
	                                        maxdim = npy.nanmax(ssh.shape);
	                                        coords_1D = cy + 100*maxdim*maxdim*cx;
	                                        maskarr_new[nextx,nexty] = 1;
	                                        if nexty + 100*maxdim*maxdim*nextx not in coords_1D:
	                                            cy=npy.append(cy,nexty);
	                                            cx=npy.append(cx,nextx);
	                                            end_of_loop = end_of_loop + 1;
	        else:
	            near_ocean = 0; looping = end_of_loop;
	
	        looping = looping + 1;
	        size_of_new_mask = npy.nansum(maskarr_new);
	        maskarr_new_withland = maskarr_new.copy();
	        maskarr_new_withland[land_arr==cond] = 0;
	        size_of_new_mask_withland = npy.nansum(maskarr_new_withland);
	    
	    new_edges = FND_CEDGE (maskarr_new)+maskarr_new;
	    new_edges_withland = new_edges.copy(); new_edges_withland[land_arr==cond] = 0;
	    #########################
	    #check if its reached land or not
	    if maxarr < npy.nanmin(ssh):
	        all_met = 1; maskarr_out = maskarr.copy();
	    
	    if abs(increment) > abs(inc_min):        
	        all_met = 0;
	        if near_ocean == 0:
	            maxarr = maxarr + increment;
	            increment = increment/2;
	            maxarr = maxarr - increment;
	            maskarr_new = maskarr.copy();
	            print('		non-ocean cells. Try lower increment')
	        elif npy.nansum(new_edges_withland) < npy.nansum(new_edges):
	            maxarr = maxarr + increment;
	            increment = increment/2;
	            maxarr = maxarr - increment;
	            maskarr_new = maskarr.copy();
	            print('		met a wall? reversing')
	        elif size_of_new_mask > size_of_old_mask:
	            print('		continuing to increment out')
	            maxarr = maxarr - increment;
	            maskarr = maskarr_new.copy();
	            size_of_old_mask = size_of_new_mask.copy();
	        else:
	            maxarr = maxarr + increment;
	            increment = increment/2;
	            maxarr = maxarr - increment;
	            maskarr_new = maskarr.copy();
	            print('		other')
	    else:
	        all_met = 1; maskarr_out = maskarr.copy();
	    
	    print('		next '+str(while_loop)+' max '+str(maxarr)+' no vals '+str(npy.sum(maskarr_new))+' inc '+str(increment))
	
	############################
	mask_full[x1:x2,y1:y2] = maskarr_out.copy();
	
	BGcalcarr = mask_full*ssh_full;
	BGcalcarr[BGcalcarr==0] = npy.nan
	BGcalcmin = npy.nanmin(BGcalcarr)
	
	## final metrics for netcdf
	mask_nan = mask_full*1
	mask_nan[mask_full==0] = npy.nan
	BG_max_val = npy.nanmax(mask_nan*ssh_full)
	[r,c] = npy.where(mask_nan*ssh_full == BG_max_val)
	print('			',r,c,len(r))
	BG_max_lat = lat[r[0],c[0]]
	BG_max_lon = lon[r[0],c[0]]
	BG_area = npy.nansum((grid_area*mask_nan)[:]) 
	  
	## netcdf file output
	mask_nan = mask_full*1; mask_full[mask_full!=1] = npy.nan
	
	return mask_full, BG_max_val, BGcalcmin, BG_max_lat, BG_max_lon, BG_area  

################################################################################################################################
def FND_CEDGE ( oldarr ) :
################################################################################################################################
	"""
	Function identifies the edge of a contour by looking at the four adjacent cells

	History: This algorithm has been originaly developed by Heather Regan and slightly adapeted to be included into the MONARC
		 See Regan et al. JPO2020 ; https://doi.org/10.1175/JPO-D-19-0234.1
	
	Input:
	    oldarr : mask array 
	
	Output:
	    newarr : a new mask for contour 
	"""
	# ----------------------------------------------------------------------

	#Oldarr is a mask. Newarr finds coordinates next to ones
	newarr = npy.zeros(oldarr.shape);    
	
	dx = npy.diff(oldarr,axis=0);
	dy = npy.diff(oldarr,axis=1);
	y_offset_top = npy.zeros(dy.shape);
	y_offset_bottom = npy.zeros(dy.shape);
	x_offset_left = npy.zeros(dx.shape);
	x_offset_right = npy.zeros(dx.shape);
	    
	y_offset_top[dy==1] = 1; 
	y_offset_bottom[dy==-1] = 1;
	x_offset_left[dx==1] = 1;
	x_offset_right[dx==-1] = 1;
	    
	## putting it into new mask
	newmask_y = npy.zeros(newarr.shape);
	newmask_x = npy.zeros(newarr.shape);
	newmask_y[:,0] = y_offset_top[:,0];
	newmask_y[:,-1] = y_offset_bottom[:,-1];
	newmask_y[:,1:-2] = y_offset_top[:,1:-1] + y_offset_bottom[:,0:-2];
	newmask_x[0,:] = x_offset_left[0,:];
	newmask_x[-1,:] = x_offset_right[-1,:];
	newmask_x[1:-2,:] = x_offset_left[1:-1,:] + x_offset_right[0:-2,:];
	   
	newarr[newmask_x+newmask_y > 0] = 1; 
	    
	return newarr

################################################################################################################################
def EKE_CALC( zlon, zlat, zCONF, zCASE, xiosfreq, zc_year, zdatadir, zncout ) :
################################################################################################################################
	"""
	Function dedicated to compute EKE using monthly a yearly velocities 
	
	Input:
	    zlon     : longitude 2D  
	    zlat     : latitude 2D  
	    zCONF    : the configuration name 
	    zCASE    : the experiment name associated to the configuration
	    xiosfreq : model outputs frequency 
	    zc_year  : the current year 
	    zdatadir : path to access model outputs 
	    zncout   : logical to outputs (or not) results into a netcdf file 
	
	Output:
	    ds_eke : annual mean 3D EKE in a Dataset
	"""
	# ----------------------------------------------------------------------

	# Prepare all metrics 
	datadir = Path('./'+zCONF+'/GRID/')
	domcfg = open_domain_cfg(datadir=datadir, files=[zCONF+'_domain_cfg.nc'])
	
	metrics = { #define the name of the scaling factors
	    ('X',): ['e1t', 'e1u', 'e1v', 'e1f'], # X distances
	    ('Y',): ['e2t', 'e2u', 'e2v', 'e2f'], # Y distances
	    ('Z',): ['e3t_0', 'e3u_0', 'e3v_0', 'e3f_0', 'e3w_0'], # Z distances
	}   
	grid = xgcm.Grid(domcfg, metrics=metrics, periodic=False) #create the grid
	
	# Read one file with T-point variables 
	#locpath=Path(zdatadir+'/'+str(zc_year)+'/'+xiosfreq+'/')
	zpath=zdatadir+'/'+str(zc_year)+'/'+xiosfreq+'/'
	locfile=zCONF+'-'+zCASE+'_y'+str(zc_year)+'.'+xiosfreq+'_gridT.nc'
	nemoT = xr.open_dataset(zpath+locfile, engine="netcdf4")[['nav_lon','nav_lat','deptht','time_counter']]
	nemoT = nemoT.rename({'deptht':'z'})

	# Read yearly [UV] files 
	nemo_yyUV = process_nemo(
	    positions=[
	        (xr.open_mfdataset(zpath+zCONF+'-'+zCASE+'_y'+str(zc_year)+'.'+xiosfreq+'_gridU.nc',concat_dim=["time_counter"], combine='nested', parallel=True)[['vozocrtx']], 'U'),
	        (xr.open_mfdataset(zpath+zCONF+'-'+zCASE+'_y'+str(zc_year)+'.'+xiosfreq+'_gridV.nc',concat_dim=["time_counter"], combine='nested', parallel=True)[['vomecrty']], 'V'),
	    ],
	    domcfg=domcfg
	)
	
	# Read all monthly [UV] files 
	#locpath=Path(zdatadir+'/'+str(c_year)+'/'+xiosfreq+'/')
	zpath=zdatadir+'/'+str(zc_year)+'/'+xiosfreq+'/'
	#locfile=zCONF+'-'+zCASE+'_y'+str(zc_year)+'m*.'+xiosfreq+'_grid[UV].nc'
	#if chkfile(locpath+locfile) : 
	#	nemo_mmUV = open_nemo(domcfg=domcfg, files=locpath.glob(locfile))
	#	nemo_mmUV
	nemo_mmUV = process_nemo(
	    positions=[
	        (xr.open_mfdataset(zpath+zCONF+'-'+zCASE+'_y'+str(zc_year)+'m*.'+xiosfreq+'_gridU.nc',concat_dim=["time_counter"], combine='nested', parallel=True)[['vozocrtx']], 'U'),
	        (xr.open_mfdataset(zpath+zCONF+'-'+zCASE+'_y'+str(zc_year)+'m*.'+xiosfreq+'_gridV.nc',concat_dim=["time_counter"], combine='nested', parallel=True)[['vomecrty']], 'V'),
	    ],
	    domcfg=domcfg
	)
	
	# Compute the annual mean EKE using monthly mean velocities
	###############################################################
	# Velocities anomalies against annual mean (The following syntax is for keeping the dimensions order as [t,z,y,x] )
	Up_mmU = -1.*(nemo_mmUV['vozocrtx'] - nemo_yyUV['vozocrtx'].squeeze())
	Vp_mmV = -1.*(nemo_mmUV['vomecrty'] - nemo_yyUV['vomecrty'].squeeze())
	
	# Velocity squared at [UV]-point
	Up2_mmU = Up_mmU**2 
	Vp2_mmV = Vp_mmV**2
	
	# Velocity squared at T-point
	Up2_mmT = grid.interp(Up2_mmU,axis='X')
	Vp2_mmT = grid.interp(Vp2_mmV,axis='Y')
	
	# Monthly EKE at T-point 
	EKE_mmT = 0.5 * ( Up2_mmT + Vp2_mmT )
	
	# Save EKE annual mean in a dataset with appropriate coordinates 
	ds_eke = xr.Dataset()
	ds_eke = ds_eke.assign_coords( z=nemoT.z, longitude=zlon, latitude=zlat )
	ds_eke['voeke'] = (('z','y','x'), EKE_mmT.mean(dim='t').values)

	return ds_eke 

################################################################################################################################
def CONV_SA2PS( zSAL, zlon, zlat ) :
################################################################################################################################
	"""
	Function to convert CT/SA to Tpot/PS units using the GSW package 
	
	Input:
	    zSAL : 3D temperature/Salinity data array 
	    zlon : 2D longitude 
	    zlat : 2D latitude
	
	Output:
	    z_SP : 3D practical salinity 
	"""

	# Compute the pressure at each depth
	pressure = gsw.p_from_z( -zSAL.z.values.squeeze(), 77. )
	pressure3D = P3Dzyx( pressure, zlat ) 
	
	# Apply the conversion
	z_SP = gsw.conversions.SP_from_SA( zSAL, pressure3D, zlon, zlat )
	
	return z_SP.astype('float32')

################################################################################################################################
def P3Dzyx( zvector, zlat ):
################################################################################################################################
	"""
	Function to extend a 1D pressure to 3D 
	
	Input:
	    zvector : 1D pressure 
	    zlat    : 2D latitude 
	
	Output:
	    zpressure3D : 3D pressure 
	"""

	# Prepare this 1D field to be duplicated in 3D : z,y,x
	z2dt = npy.reshape( zvector, (len(zvector), 1, 1) )
	
	# Horizontaly to fit the T/S on a global grid
	zpressure3D = npy.tile( z2dt,( 1, zlat.shape[0], zlat.shape[1] ) )
	
	return zpressure3D

################################################################################################################################
def MLD_CALC( zTS, zlon, zlat, depthW, teos10, dtype='None' ) :
################################################################################################################################
	"""
	Function to compute the MLD using a density criteria of 0.1 kg/m3 
	
	Input:
	    zTS    : 3D temperature/salinity in March 
	    zlon   : 2D longitude 
	    zlat   : 2D latitude
	    depthW : W-points depth
	    teos10 : logical to specify the T/S nature, either EOS10 or EOS80
	    dtype  : a name to set a specific treatment 
	
	Output:
	    z_MLD : MLD 
	"""

	# Set the densitity criteria value [kg m3]
	rn_crit = 0.1 

	# Compute the Sigma_0 potential density from the surface 
	if teos10 : 
	#if not teos10 : 
		Sigma_0 = gsw.sigma0( zTS['vosaline'].values, zTS['votemper'].values ) 
	else :
		# Need to compute SA & CT fields first
		# Compute the pressure at each depth
		pressure = gsw.p_from_z( -zTS.z.values.squeeze(), 77. )
		pressure3D = P3Dzyx( pressure, zlat ) 

		SA = gsw.SA_from_SP(zTS['vosaline'].values, pressure3D, zlon.values, zlat.values)
		CT = gsw.CT_from_pt(SA, zTS['votemper'].values)

		Sigma_0 = gsw.sigma0( SA, CT )

	# Compute the density anomaly against the surface one
	Sig0_ano = Sigma_0 - Sigma_0[0,:,:]

	# Detect the MLD using a 0.1 kg/m3 density criteria
	Sig0_ano_msk = ( Sig0_ano <= rn_crit )                        # Create a mask to identify T-points in the MLD 
	if dtype == 'woa09' : Sig0_ano_msk[55:74] = False             # Specific tretment for woa09 data because data at depth are not masked 
	Sig0_ano_msk_roll = npy.roll( Sig0_ano_msk, shift=1, axis=0 ) # Shift by 1 level over depth to get the effective W-point 
	DeptW_mld = depthW.where( Sig0_ano_msk_roll )                 # Select only W-point depth in the MLD  
	z_MLD = DeptW_mld.max(dim='z',skipna=True).compute()          # Finally, found the greater depth over depth
	
	return z_MLD.astype('float32')

################################################################################################################################
def TOPOS_CALC( zlon, zlat, zCONF, zCASE, zc_year, zdatadir ) :
################################################################################################################################
	"""
	Function dedicated to compute topostrophy using monthly velocities 
	
	Input:
	    zlon     : longitude 2D  
	    zlat     : latitude 2D  
	    zCONF    : the configuration name 
	    zCASE    : the experiment name associated to the configuration
	    zc_year  : the current year 
	    zdatadir : path to access model outputs 
	
	Output:
	    ds_tsphy : annual mean 2D Topostrophy in a Dataset
	"""
	# ----------------------------------------------------------------------

	# Prepare all metrics 
	datadir = Path('./'+zCONF+'/GRID/')
	domcfg = open_domain_cfg(datadir=datadir, files=[zCONF+'_domain_cfg.nc'])
	
	metrics = { #define the name of the scaling factors
	    ('X',): ['e1t', 'e1u', 'e1v', 'e1f'], # X distances
	    ('Y',): ['e2t', 'e2u', 'e2v', 'e2f'], # Y distances
	    ('Z',): ['e3t_0', 'e3u_0', 'e3v_0', 'e3f_0', 'e3w_0'], # Z distances
	}   
	grid = xgcm.Grid(domcfg, metrics=metrics, periodic=False) #create the grid

	# Need to get the bathymetry
	fileid = zCONF+'_domain_cfg.nc'
	ds_bathy = xr.open_dataset(datadir / fileid )[['bathy_meter','time_counter']]
	ds_bathy['time_counter'].attrs['bounds'] = ""  
	ds_bathy = ds_bathy.swap_dims({'t':'time_counter'})
	ds_bathy
		
	nemoBat = process_nemo( positions=[(ds_bathy, 'T')],domcfg=domcfg)
	nemoBat
	
	# Set the proper coordinates for each mask
	fileid = zCONF+'-'+zCASE+'_mask.nc'
	ds_mesh=xr.open_dataset(datadir / fileid )[['tmask']]
	ds_mesh = ds_mesh.rename({'nav_lev': 'z_c'}).assign_coords(z_c=domcfg.z_c.values)
	ds_mesh['time_counter'].attrs['bounds'] = ""  
	
	ds_Tmsk = process_nemo( positions=[(ds_mesh.tmask, 'T')],domcfg=domcfg).squeeze()
	ds_Tmsk = ds_Tmsk.reset_coords('t', drop=True)
	
	ds_mesh=xr.open_dataset(datadir / fileid )[['umask']]
	ds_mesh = ds_mesh.rename({'nav_lev': 'z_c'}).assign_coords(z_c=domcfg.z_c.values)
	ds_mesh['time_counter'].attrs['bounds'] = ""  
	
	ds_Umsk = process_nemo( positions=[(ds_mesh.umask, 'U')],domcfg=domcfg).squeeze()
	ds_Umsk = ds_Umsk.reset_coords('t', drop=True)
	
	ds_mesh=xr.open_dataset(datadir / fileid )[['vmask']]
	ds_mesh = ds_mesh.rename({'nav_lev': 'z_c'}).assign_coords(z_c=domcfg.z_c.values)
	ds_mesh['time_counter'].attrs['bounds'] = ""  
	
	ds_Vmsk = process_nemo( positions=[(ds_mesh.vmask, 'V')],domcfg=domcfg).squeeze()
	ds_Vmsk = ds_Vmsk.reset_coords('t', drop=True)
	
	# Read monthly velocities components 
	nemo_mmUV = open_nemo(domcfg=domcfg, files=fs.glob(zdatadir+zCONF+'-'+zCASE+'_y'+str(zc_year)+'m*.1m_grid[UV].nc'))[['vozocrtx','vomecrty']]
	
	#1 Identify the last ocean T-point 
	mask_bot_ano = (ds_Tmsk.tmask.roll(z_c=-1) - ds_Tmsk.tmask ).squeeze()
	mask_bot = ( mask_bot_ano == -1)  
	
	#2 Bathymetry horizontal gradient 
	# bathymetry interpolated at U-point and compute the zonal gradient at T-point
	dbathy_dx = grid.derivative(grid.interp(nemoBat.bathy_meter,axis='X',to='right'),axis='X').squeeze()
	# bathymetry interpolated at V-point and compute the meridional gradient at T-point
	dbathy_dy = grid.derivative(grid.interp(nemoBat.bathy_meter,axis='Y',to='right'),axis='Y').squeeze()
	
	#3 Set velocities at T-point taking into account land point 
	umsk_T = 1./ (ds_Umsk.umask.roll(x_f=1) + ds_Umsk.umask )
	vmsk_T = 1./ (ds_Vmsk.vmask.roll(y_f=1) + ds_Vmsk.vmask )
	
	Ut_ztmp = (nemo_mmUV.vozocrtx.roll(x_f=1) + nemo_mmUV.vozocrtx ) * umsk_T
	Vt_ztmp = (nemo_mmUV.vomecrty.roll(y_f=1) + nemo_mmUV.vomecrty ) * vmsk_T
	
	# Set the proper dimensions names as coordinates as well i.e. x_c etc ..
	U_T= xr.DataArray( Ut_ztmp , dims=['t','z_c','y_c','x_c'], 
	                  coords={'t':nemo_mmUV['t'].values,
	                          'z_c':domcfg['z_c'].values, 
	                          'y_c':domcfg['y_c'].values, 
	                          'x_c':domcfg['x_c'].values})
	V_T= xr.DataArray( Vt_ztmp , dims=['t','z_c','y_c','x_c'], 
	                  coords={'t':nemo_mmUV['t'].values,
	                          'z_c':domcfg['z_c'].values, 
	                          'y_c':domcfg['y_c'].values, 
	                          'x_c':domcfg['x_c'].values})
	
	#4 Select bottom velocoties only 
	Ubot_zT = U_T.where(mask_bot)
	Vbot_zT = V_T.where(mask_bot)
	mask = ~npy.isnan(Ubot_zT)
	# idz 3D (t,y,x) of bottom indices 
	idz = mask.argmax(dim='z_c')  # Trouve l'index de la première valeur non-NaN (ici 1)
	
	# .compute() is required to get a nd.array field to fid .isel (.persist would also work)
	Ubot_T = Ubot_zT.isel(z_c=idz.compute())
	Vbot_T = Vbot_zT.isel(z_c=idz.compute())
	
	#5 Final calculation of fxU . grad H
	topos_T= -1.*domcfg.ff_t*Vbot_T*dbathy_dx + domcfg.ff_t*Ubot_T*dbathy_dy

	# Save annual mean topostrophy in a dataset with appropriate coordinates 
	ds_tsphy = xr.Dataset()
	ds_tsphy = ds_tsphy.assign_coords( t=nemo_mmUV.t, longitude=zlon, latitude=zlat )
	ds_tsphy['topos'] = (('y','x'), topos_T.mean(dim='t').values)

	return ds_tsphy 
