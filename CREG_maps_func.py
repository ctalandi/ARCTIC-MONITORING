#!/usr/bin/env python

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
import pandas as pd
import calendar
import time
from pathlib import Path 
from xnemogcm import open_domain_cfg, open_nemo, process_nemo, open_namelist, open_nemo_and_domain_cfg
from xnemogcm import __version__ as xnemogcm_version
import xgcm


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
def BFG_mapsf( zlon, zlat, zvar_ssh, zbathy, zarea, zCONF, zCASE, zs_year, ze_year, zlgTS_ys, zlgTS_ye, zncout ) :
################################################################################################################################

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
	msk_ym, BG_maxval_ym, BG_minval_ym, BG_maxlat_ym, BG_maxlon_ym, BG_area_ym = BFG_compute( zlon, zlat, zssh_ym, zbathy, 'SSH', zincr, zarea )
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
		msk[zmm,:,:], BG_maxval[zmm], BG_minval[zmm], BG_maxlat[zmm], BG_maxlon[zmm], BG_area[zmm] = BFG_compute( zlon, zlat, zssh, zbathy, 'SSH', zincr, zarea )
		print('			Computing wall time (s) : ',time.time() - start_time)

	#------------------------------------------------------------------------------------------------------------------------

        # Plot the yearly closed contours as the BFG center as well 
	plt.figure()
	plt.subplot(211)
	zoutmap, X, Y = Iso_Bat( ztype='isol1000', zarea='cassis_BGZoom' )	
	zoutmap.drawparallels(npy.arange(-90.,91.,2.),labels=[True,False,False,False], size=5, linewidth=0.3)
	zoutmap.drawmeridians(npy.arange(-180.,181.,10.),labels=[False,False,False,True], size=5, latmax=90.,linewidth=0.3)
	zoutmap.fillcontinents(color='grey',lake_color='white')

	if npy.nansum(msk_ym) > 0:
		msk_plot = xr.where( npy.isnan(msk_ym), 0., msk_ym*1 )
		CS2 = zoutmap.contour( X, Y, msk_plot, linewidths=0.5, colors='k' )
	# Get indices of the BFG center 
	[r,c] = npy.nonzero( msk_ym*zssh_ym == npy.nanmax(msk_ym*zssh_ym) )
	
	clat = [zlat[r.item(),c.item()],]
	clon = [zlon[r.item(),c.item()],]
	cx,cy = zoutmap(clon,clat)
	zoutmap.scatter(cx,cy, s=10, marker='o', color='k')
	plt.title( zCASE+' BFG SSH contours \n yearly mean SSH '+str(zs_year), fontsize=6 )


        # Plot the monthly closed contours as the BFG center as well 
	cmap = plt.get_cmap('Spectral_r')
	colors = [cmap(i) for i in npy.linspace(0, 1, 12)]

	plt.subplot(212)
	zoutmap, X, Y = Iso_Bat( ztype='isol1000', zarea='cassis_BGZoom' )	
	zoutmap.drawparallels(npy.arange(-90.,91.,2.),labels=[True,False,False,False], size=5, linewidth=0.3)
	zoutmap.drawmeridians(npy.arange(-180.,181.,10.),labels=[False,False,False,True], size=5, latmax=90.,linewidth=0.3)
	zoutmap.fillcontinents(color='grey',lake_color='white')

	for zmm in range(0,12):
		if npy.nansum(msk[zmm,:,:]) > 0:
			msk_plot = xr.where( npy.isnan(msk[zmm,:,:]), 0., msk[zmm,:,:]*1 )
			CS2 = zoutmap.contour( X, Y, msk_plot, linewidths=0.5, colors=colors[zmm] )
		# Get indices of the BFG center 
		[r,c] = npy.nonzero( msk[zmm,:,:]*zvar_ssh.isel(time_counter=zmm) == npy.nanmax(msk[zmm,:,:]*zvar_ssh.isel(time_counter=zmm)) )
		
		clat = [zlat[r.values.item(),c.values.item()],]
		clon = [zlon[r.values.item(),c.values.item()],]
		cx,cy = zoutmap(clon,clat)
		zoutmap.scatter(cx,cy, s=10, marker='o', color=colors[zmm])

	legend_elements = [ Line2D([0], [0], color=colors[i], lw=1, label=calendar.month_name[i+1]) for i in range(len(colors)) ]
	plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 0.85), ncol=1, fontsize=6, handlelength=1.5, handletextpad=0.5, borderpad=0.2, labelspacing=0.2)
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

		#ds_out = ds_out.set_coords(['BGlat','BGlon','time'])
		
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
		mod_max = ds_rssh["BGmax"].values.squeeze()
		mod_area = ds_rssh["BGarea"].values.squeeze()
		mod_maxhlat = ds_rssh["BGmaxhlat"].values.squeeze()
		mod_maxhlon = ds_rssh["BGmaxhlon"].values.squeeze()
	
		# Make plots 
		#######################################################################################################
		plt.clf()
		xwind=410
		
		# Max SSH 
		################
		ax=plt.subplot(xwind+1)
		ax.set_title(zCASE,size=7)
		ds_rssh['BGmax'].plot(color='k', linewidth=0.6, label='model')
		ds_obs['maxheight'].plot(color='g', linewidth=0.6, label='obs')
		plt.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
		plt.xlabel('years' ,size=5)
		plt.ylabel(' max ssh at gyre centre \n (metres)',size=6)
		plt.setp(ax.get_xticklabels(),visible=False)
		plt.yticks(size=6)
		plt.legend(fontsize='small', ncol=2)

		# Gyre area 
		################
		ax=plt.subplot(xwind+2)
		(ds_rssh['BGarea']*1e-12).plot(color='k',  linewidth=0.6)
		(ds_obs['area_m2']*1e-12).plot(color='g', linewidth=0.6)
		plt.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
		plt.xlabel('years' ,size=5)
		plt.ylabel('gyre area \n'+'x10^6 (km^2)',size=6)
		plt.setp(ax.get_xticklabels(),visible=False)
		plt.yticks(size=6)
	
		# Latitude of Max SSH 
		######################
		ax=plt.subplot(xwind+3)
		ds_rssh['BGmaxhlat'].plot(color='k', linewidth=0.6)
		ds_obs['maxh_lat'].plot(color='g',linewidth=0.6)
		plt.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
		plt.xlabel('years' ,size=5)
		plt.ylabel('latitude of max ssh',size=6)
		plt.setp(ax.get_xticklabels(),visible=False)
		plt.yticks(size=6)
	
		# Longitude of Max SSH 
		######################
		ax=plt.subplot(xwind+4)
		ds_rssh['BGmaxhlon'].plot(color='k', linewidth=0.6)
		ds_obs['maxh_lon'].plot(color='g', linewidth=0.6)
		plt.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
		plt.xlabel('years' ,size=5)
		plt.ylabel('longitude of max ssh',size=6)
		plt.xticks(size=5)
		plt.yticks(size=6)
	
		plt.tight_layout()

		zfile_ext='_BFG_metrics_LGTS_y'+str(zlgTS_ys)+'LASTy'
		plt.savefig(zCONF+'-'+zCASE+zfile_ext+'.png',dpi=300)

	return

################################################################################################################################
def ICE_maps( zlon, zlat, zMy_var1, zMy_var1frld_SeasM, zMy_var1frld_SeasS, zCONF, zCASE, zclimyear, zs_year, zc_year, zncout ) :
################################################################################################################################

        num_fram=320
        # Annual mean Ice thickness
        zMyvar='sivolu'   ; fram=num_fram+1
        zMy_var1 = xr.where( zMy_var1 == 0., npy.nan, zMy_var1 )

        simple_maps( zlon, zlat, zCONF, zCASE, zMy_var1, zMyvar, zclimyear, zfram=fram )
        # March mean Ice fraction 
        zMyvar='siconc'   ; fram=num_fram+3
        simple_maps( zlon, zlat, zCONF, zCASE, zMy_var1frld_SeasM, zMyvar, zclimyear, seas='m03', zfram=fram )
        # September mean Ice fraction 
        zMyvar='siconc'   ; fram=num_fram+5
        simple_maps( zlon, zlat, zCONF, zCASE, zMy_var1frld_SeasS, zMyvar, zclimyear, seas='m09', zfram=fram )

        # Annual mean Ice thickness PIOMASS observations
        # WARNING this part has been interpolated on CREG025 grid directly 
        # Read lat, lon from CREG025 grid since the PIOMASS data is not available on CREG12.L75 grid
        zmask, plon, plat = CREG_MSK( zCASE )

        obs_thick = ICE_THICK_OBS( zconfig=zCONF, t_year=zs_year )
        obs_thick = xr.where( npy.squeeze(zmask[0,:,:]) < 1., npy.nan, obs_thick )
        obs_thick = xr.where( obs_thick == 0., npy.nan, obs_thick )

        zMyvar='sivolu'   ; fram=num_fram+2
        simple_maps( plon, plat, zCONF, zCASE, obs_thick, zMyvar, zs_year, zfram=fram, plot_obs=1 )

        # Read NSIDC obs. data 
        obs_conc_m03, obs_conc_m09, obs_lon, obs_lat = ICE_CONCE_OBS( t_year=zc_year )

        # March mean Ice fraction 
        zMyvar='siconc'   ; fram=num_fram+4
        simple_maps( obs_lon, obs_lat, zCONF, zCASE, obs_conc_m03, zMyvar, zs_year, seas='m03', zfram=fram, plot_obs=1 )
        # September mean Ice fraction 
        zMyvar='siconc'   ; fram=num_fram+6
        simple_maps( obs_lon, obs_lat, zCONF, zCASE, obs_conc_m09, zMyvar, zs_year, seas='m09', zfram=fram, plot_obs=1 )
        plt.tight_layout()

        zfile_ext='_ICEClim_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

        if zncout:
                ds_out = xr.Dataset()
                
                # ICE fields
                #######################
                ds_out['IceThick_mod']= (('y','x'), zMy_var1.values.astype('float32')) 
                ds_out['IceThick_mod'].attrs['long_name']='Model annual mean Ice thickness'
                ds_out['IceThick_mod'].attrs['units']='m'
                
                ds_out['IceConceM03_mod']= (('y','x'), zMy_var1frld_SeasM.values.astype('float32')) 
                ds_out['IceConceM03_mod'].attrs['long_name']='Model monthly mean Ice concentration in march'
                ds_out['IceConceM03_mod'].attrs['units']='-'
                
                ds_out['IceConceM09_mod']= (('y','x'), zMy_var1frld_SeasS.values.astype('float32')) 
                ds_out['IceConceM09_mod'].attrs['long_name']='Model monthly mean Ice concentration in september'
                ds_out['IceConceM09_mod'].attrs['units']='-'
                
                ds_out['IceThick_obs']= (('y','x'), obs_thick.values.astype('float32')) 
                if zs_year >= 1979 and zs_year <= 2013 :
                        ds_out['IceThick_obs'].attrs['long_name']='PIOMAS annual mean Ice thickness over '+str(zs_year)
                else :
                        ds_out['IceThick_obs'].attrs['long_name']='PIOMAS climatological mean Ice thickness over 1979-2013'
                ds_out['IceThick_obs'].attrs['units']='m'
                
                ds_out['IceConceM03_obs']= (('yobs','xobs'), obs_conc_m03.values.astype('float32')) 
                if zs_year >= 1979 and zs_year <= 2015 :
                        ds_out['IceConceM03_obs'].attrs['long_name']='NSDIC monthly mean Ice concentration in march '+str(zs_year)
                else :
                        ds_out['IceConceM03_obs'].attrs['long_name']='NSDIC climatological mean Ice concentration in march over 1979-2015'
                ds_out['IceConceM03_obs'].attrs['units']='-'
                
                ds_out['IceConceM09_obs']= (('yobs','xobs'), obs_conc_m09.values.astype('float32')) 
                if zs_year >= 1979 and zs_year <= 2015 :
                        ds_out['IceConceM09_obs'].attrs['long_name']='NSDIC monthly mean Ice concentration in september '+str(zs_year)
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
                nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_ICEClim_'+'y'+zclimyear+'.nc'
                ds_out.to_netcdf(nc_f,engine='netcdf4')

        return

################################################################################################################################
def MLD_maps( zlon, zlat, zMy_var1SeasM, zMy_var1SeasS, zCONF, zCASE, zclimyear, zncout ) :
################################################################################################################################

        num_fram=220
        # March mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+1
        simple_maps( zlon, zlat, zCONF, zCASE, zMy_var1SeasM, zMyvar, zclimyear, seas='m03', zfram=fram )
        # September mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+3
        simple_maps( zlon, zlat, zCONF, zCASE, zMy_var1SeasS, zMyvar, zclimyear, seas='m09', zfram=fram )

        # MLD from observation
        mld_obs_m03, mld_obs_m09, lon_obs, lat_obs = MLD_OBS()
        # March mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+2
        simple_maps( lon_obs, lat_obs, zCONF, zCASE, mld_obs_m03, zMyvar, zclimyear, seas='m03', zfram=fram, plot_obs=1 )
        # September mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+4
        simple_maps( lon_obs, lat_obs, zCONF, zCASE, mld_obs_m09, zMyvar, zclimyear, seas='m09', zfram=fram, plot_obs=1 )
        plt.tight_layout()

        zfile_ext='_MLDClim_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

        if zncout:
                ds_out = xr.Dataset()
                
                # MLD fields
                #######################
                ds_out['MLDd001M03_mod']= (('y','x'), zMy_var1SeasM.values.astype('float32')) 
                ds_out['MLDd001M03_mod'].attrs['long_name']='Model monthly mean MLD in march based on a density criteria of 0.01 kg/m^3'
                ds_out['MLDd001M03_mod'].attrs['units']='m'

                ds_out['MLDd001M09_mod']= (('y','x'), zMy_var1SeasS.values.astype('float32')) 
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
                nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_MLDClim_'+'y'+zclimyear+'.nc'
                ds_out.to_netcdf(nc_f,engine='netcdf4')

        return

################################################################################################################################
def ENE_maps( zlon, zlat, zMy_var1, ds_eke, zdepth, zCONF, zCASE, zclimyear, zs_year, zncout ) :
################################################################################################################################

        num_fram=220
        # mean PSI
        zMyvar='sobarstf'   ; fram=num_fram+4
        simple_maps( zlon, zlat, zCONF, zCASE, zMy_var1, zMyvar, zclimyear, zfram=fram )

        # Get the model depth at the surface and at ~100m 
        lev0 =  0 # ~  0m in the model
        lev1 = 20 # ~ 69m in the model; corresponds to the halocline depth 
        lev2 = 39 # ~503m in the model; corresponds to the AW depth 
        zd0 = int(zdepth.isel(nav_lev=lev0).values.item()) ;  zd1 = int(zdepth.isel(nav_lev=lev1).values.item())  ; zd2 = int(zdepth.isel(nav_lev=lev2).values.item())

        # EKE inferred from Von Appen et al. VONAPPEN_EKE_OBS()
        obs_VonAppeneke = VONAPPEN_EKE_OBS( )
        #obs_vonapeneke = xr.where( obs_eke >= 9e20, npy.nan, obs_eke )

        # Mean EKE at the surface 
        zMyvar='voeke'   ; fram=num_fram+1
        m = simple_maps( zlon, zlat, zCONF, zCASE, npy.log10(ds_eke[zMyvar].isel(z=lev0)), zMyvar, zclimyear, slev=str(zd0) , zfram=fram )

        # EKE from DOT observations (Armitage et al. 2017)
        obs_eke, lon_obs, lat_obs = EKE_OBS( t_year=zs_year )
        obs_eke = xr.where( obs_eke >= 9e20, npy.nan, obs_eke )

        zMyvar='voeke'   ; fram=num_fram+2
        simple_maps( lon_obs, lat_obs, zCONF, zCASE, npy.log10(obs_eke), zMyvar, zclimyear, slev=str(zd0), zfram=fram, plot_obs=1 )
        #plt.tight_layout()
        #plt.subplots_adjust(wspace=0.05)

        zfile_ext='_DYNPSIClim_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)


        plt.clf()
        num_fram=120
        # Mean EKE at 69m 
        zMyvar='voeke'   ; fram=num_fram+1
        m = simple_maps( zlon, zlat, zCONF, zCASE, npy.log10(ds_eke[zMyvar].isel(z=lev1)), zMyvar, zclimyear, slev=str(zd1) , zfram=fram )
	# Select data in the halocline depth
        mask_halo = npy.logical_and(obs_VonAppeneke['Mean depth'] >= 50, obs_VonAppeneke['Mean depth'] <= 100)
        lons_halo = obs_VonAppeneke.Longitude.where(mask_halo)
        lats_halo = obs_VonAppeneke.Latitude.where(mask_halo)
        X_halo, Y_halo = m(lons_halo, lats_halo)
        eke_halo = obs_VonAppeneke.EKE.where(mask_halo)

        vmin = -6; vmax = -2
        cmap = plt.get_cmap('RdYlBu_r')
        plt.scatter(X_halo, Y_halo, s = 10, c = npy.log10(eke_halo), cmap = cmap, vmin = vmin, vmax = vmax, edgecolors = 'k', linewidths = 0.5)

        # Mean EKE at 508m 
        zMyvar='voeke'   ; fram=num_fram+2
        m = simple_maps( zlon, zlat, zCONF, zCASE, npy.log10(ds_eke[zMyvar].isel(z=lev2)), zMyvar, zclimyear, slev=str(zd2), zfram=fram )
        # Select data in the AW layer 
        mask_aw = ~npy.isnan(obs_VonAppeneke['EKE at depth'])
        lons_aw = obs_VonAppeneke.Longitude.where(mask_aw)
        lats_aw = obs_VonAppeneke.Latitude.where(mask_aw)
        X_aw, Y_aw = m(lons_aw, lats_aw)
        eke_aw = obs_VonAppeneke['EKE at depth'].where(mask_aw)

        vmin = -6; vmax = -2
        cmap = plt.get_cmap('RdYlBu_r')
        plt.scatter(X_aw, Y_aw, s = 10, c = npy.log10(eke_aw), cmap = cmap, vmin = vmin, vmax = vmax, edgecolors = 'k', linewidths = 0.5)

        plt.tight_layout()

        zfile_ext='_DYNEKEClim_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

        if zncout:
                ds_out = xr.Dataset()
                
                # DYN fields
                #######################
                ds_out['PSI_mod']= (('y','x'), zMy_var1.values.astype('float32')) 
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

                ds_out['EKESurf_obs']= (('yobs','xobs'), obs_eke.values.astype('float32')) 
                if zs_year >= 2003 and zs_year <= 2014 :
                        ds_out['EKESurf_obs'].attrs['long_name']='EKE annual mean derived from DOT field (Armitage et al. 2017) in '+str(zs_year)
                else :
                        ds_out['EKESurf_obs'].attrs['long_name']='EKE climatological mean derived from DOT field (Armitage et al. 2017) over 2003-2014'
                ds_out['EKESurf_obs'].attrs['units']='m2/s2'

                ds_out['lat_obs']= (('yobs','xobs'), lat_obs.astype('float32')) 
                ds_out['lat_obs'].attrs['long_name']='Degrees north'
                ds_out['lat_obs'].attrs['units']='Deg'
                
                ds_out['lon_obs']= (('yobs','xobs'), lon_obs.astype('float32')) 
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
                nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_DYNClim_'+'y'+zclimyear+'.nc'
                ds_out.to_netcdf(nc_f,engine='netcdf4')
        return

################################################################################################################################
def TSD_maps( zlon, zlat, zMy_var1T, zMy_var1S, zMy_varTinit, zMy_varSinit, zdepth, zCONF, zCASE, zclimyear, zncout ) :
################################################################################################################################

        # Get the model depth at the surface and at ~100m 
        zd1 = int(zdepth.isel(nav_lev=0).values.item())  ; zd2 = int(zdepth.isel(nav_lev=23).values.item())
        #zd1 = npy.round(zdepth.isel(nav_lev=0).values).item()  ; zd2 = npy.round(zdepth.isel(nav_lev=23).values).item()

        num_fram=220
        # Surface temperature
        zMyvar='votemper'   ; fram=num_fram+1
        simple_maps( zlon, zlat, zCONF, zCASE, npy.squeeze(zMy_var1T[0,:,:]-zMy_varTinit[0,:,:]), zMyvar, zclimyear, slev=str(zd1) , zfram=fram, ano=1 )
        # ~100m temperature
        zMyvar='votemper'   ; fram=num_fram+2
        simple_maps( zlon, zlat, zCONF, zCASE, npy.squeeze(zMy_var1T[23,:,:]-zMy_varTinit[23,:,:]), zMyvar, zclimyear, slev=str(zd2), zfram=fram, ano=1 )
        # Surface salinity
        zMyvar='vosaline'   ; fram=num_fram+3
        simple_maps( zlon, zlat, zCONF, zCASE, npy.squeeze(zMy_var1S[0,:,:]-zMy_varSinit[0,:,:]), zMyvar, zclimyear, slev=str(zd1), zfram=fram, ano=1 )
        # ~100m  salinity
        zMyvar='vosaline'   ; fram=num_fram+4
        simple_maps( zlon, zlat, zCONF, zCASE, npy.squeeze(zMy_var1S[23,:,:]-zMy_varSinit[23,:,:]), zMyvar, zclimyear, slev=str(zd2), zfram=fram, ano=1 )
        plt.tight_layout()

        zfile_ext='_TSDIffClim_@'+str(zd1)+'m@'+str(zd2)+'m_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

        plt.clf()

        # Get the model depth at ~200m and at ~300m 
        zd3 = int(zdepth.isel(nav_lev=30).values.item())  ; zd4 = int(zdepth.isel(nav_lev=34).values.item())

        num_fram=220
        # ~200m temperature
        zMyvar='votemper'   ; fram=num_fram+1
        simple_maps( zlon, zlat, zCONF, zCASE, npy.squeeze(zMy_var1T[30,:,:]-zMy_varTinit[30,:,:]), zMyvar, zclimyear, slev=str(zd3) , zfram=fram, ano=1 )
        # ~300m temperature
        zMyvar='votemper'   ; fram=num_fram+2
        simple_maps( zlon, zlat, zCONF, zCASE, npy.squeeze(zMy_var1T[34,:,:]-zMy_varTinit[34,:,:]), zMyvar, zclimyear, slev=str(zd4), zfram=fram, ano=1 )
        # ~200m salinity
        zMyvar='vosaline'   ; fram=num_fram+3
        simple_maps( zlon, zlat, zCONF, zCASE, npy.squeeze(zMy_var1S[30,:,:]-zMy_varSinit[30,:,:]), zMyvar, zclimyear, slev=str(zd3), zfram=fram, ano=1 )
        # ~300m  salinity
        zMyvar='vosaline'   ; fram=num_fram+4
        simple_maps( zlon, zlat, zCONF, zCASE, npy.squeeze(zMy_var1S[34,:,:]-zMy_varSinit[34,:,:]), zMyvar, zclimyear, slev=str(zd4), zfram=fram, ano=1 )
        plt.tight_layout()

        zfile_ext='_TSDIffClim_@'+str(zd3)+'m@'+str(zd4)+'m_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

        return

################################################################################################################################
def ATL_maps( zlon, zlat, zMy_var1SeasM, zMdata_read, zMy_var1T, zMy_varTinit, zMy_var1ssh, zdepth, zCONF, zCASE, zclimyear, zs_year, ze_year, zncout ) :
################################################################################################################################

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
        ztitle=zCASE +' mean MLD01 over \n'+zclimyear+'  m03'
        vmin=0. ; vmax=2400. ; vint=100.   ;   contours=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]
        limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]

        plt.subplot(fram)
        zoutmap = Iso_Bat( ztype='isol1000',zarea='labsea' )
        Proj_plot( zlon, zlat, zMy_var1SeasM, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='labsea' )
        plt.tight_layout()

        zfile_ext='_LAB_MLDClimm03_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

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
        plt.plot( time_axis, -1.*npy.squeeze(zMdata_read[:,j_K1:j_K1+1,i_K1:i_K1+1]), linewidth=0.8, color='k', label='Lab Sea K1' )
        # Plot obs. MLD in March
        year_obs = npy.arange(1995,2006,1)+0.20547945   ; mld_obs=[-2300.,-1300.,-1400.,-1000.,-1000.,-1100.,-1100.,-1200.,-1400.,-700.,-1300.]
        plt.scatter(year_obs,mld_obs)

        # In Irm. Sea
        if zCONF == 'CREG12.L75' : 
                i_K1=697   ;   j_K1=577   # CREG12.L75 C-type indices geo loc   60.88N  -36.99W
        else :
                i_K1=232   ;   j_K1=192   # CREG025.L75 C-type indices geo loc   60.88N  -36.99W
        plt.plot( time_axis, -1.*npy.squeeze(zMdata_read[:,j_K1:j_K1+1,i_K1:i_K1+1]), linewidth=0.8, color='g', label='Irm Sea ')
        plt.title( zCASE+' MLD 0.01 in Lab. & Irm. Seas \n '+str(zs_year)+str(ze_year), size=9 )
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
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

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
        plt.plot( time_axis,-1.*npy.squeeze(zMdata_read[:,j_K1:j_K1+1,i_K1:i_K1+1]),linewidth=0.8, color='k', label='Lab Sea DeepConv' )
        # Plot obs. MLD in March
        year_obs = npy.arange(1995,2006,1)+0.20547945   ; mld_obs=[-2300.,-1300.,-1400.,-1000.,-1000.,-1100.,-1100.,-1200.,-1400.,-700.,-1300.]
        plt.scatter(year_obs,mld_obs)
        plt.title(zCASE+' MLD 0.01 in Lab. @ -54W,58N \n '+str(zs_year)+str(ze_year),size=9)
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
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

        plt.clf()
        plt.figure()
        # MLD IN THE IRMINGER SEA IN MARCH
        ###################################
        num_fram=110
        # March mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+1
        my_cblab=r'MLD (m)'   ;   my_cmap=plt.get_cmap('Blues')
        ztitle=zCASE +' mean MLD01 over \n'+zclimyear+'  m03'
        vmin=0. ; vmax=1600. ; vint=100.   ;   contours=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
        limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
        plt.subplot(fram)
        zoutmap = Iso_Bat( ztype='isol1000',zarea='irmsea' )
        Proj_plot( zlon, zlat, zMy_var1SeasM, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='irmsea' )
        plt.tight_layout()

        zfile_ext='_IRM_MLDClimm03_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

        plt.clf()
        plt.figure()
        # MLD IN THE GIN SEAS IN MARCH
        ###################################
        num_fram=110
        # March mean MLD
        zMyvar='mldr10_1'   ; fram=num_fram+1
        my_cblab=r'MLD (m)'   ;   my_cmap=plt.get_cmap('Blues')
        ztitle=zCASE +' mean MLD01 over \n'+zclimyear+'  m03'
        vmin=0. ; vmax=1600. ; vint=100.   ;   contours=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
        limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
        plt.subplot(fram)
        zoutmap = Iso_Bat( ztype='isol1000',zarea='ginsea' )
        Proj_plot( zlon, zlat, zMy_var1SeasM, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='ginsea' )
        plt.tight_layout()

        zfile_ext='_GIN_MLDClimm03_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

        plt.clf()
        plt.figure()
        # PLOT ISOTHERM 17 Deg OFF CAPE HATTERAS
        ########################################
        num_fram=110
        zMyvar='votemper'   ; fram=num_fram+1
        my_cblab=r'ISO 17 (DegC)'   ;   my_cmap=plt.get_cmap('jet')
        ztitle=zCASE +' mean Iso 17 DegC over \n'+zclimyear
        vmin=0. ; vmax=2400. ; vint=100.   ;   contours=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]
        limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]

        plt.subplot(fram)
        klev=29
        zoutmap = Iso_Bat( ztype='isol1000', zarea='GulfS' )
        zzlon = zlon.values   ; zzlat = zlat.values 
        Proj_plot( zzlon, zzlat, npy.squeeze(zMy_var1T[klev,:,:])   , contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='GulfS' )
        Proj_plot( zzlon, zzlat, npy.squeeze(zMy_varTinit[klev,:,:]), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='GulfS', data_ref=True )
        plt.tight_layout()

        zfile_ext='_ATL_ISO17Clim_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

        plt.clf()
        plt.figure()
        # PLOT SSH OVER THE ATLANTIC AREA
        #################################
        num_fram=110
        zMyvar='sossheig'   ; fram=num_fram+1
        my_cblab=r'SSH (cm)'   ;   my_cmap=plt.get_cmap('coolwarm')
        ztitle=zCASE +' mean SSH over \n'+zclimyear
        vmin=-100. ; vmax=100. ; vint=5.  ;   contours=npy.arange(vmin,vmax+vint,vint)
        limits=[vmin,vmax,vint]           ;   myticks=npy.arange(vmin,vmax+vint,vint)

        #plt.subplot(fram)
        zoutmap = Iso_Bat( ztype='isol1000' )
        Proj_plot( zzlon,  zzlat, zMy_var1ssh*100. , contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar )
        #Proj_plot( zzlon,  zzlat, npy.squeeze(zMy_varTinit[klev,:,:]), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='GulfS', data_ref=True )
        #plt.tight_layout()

        zfile_ext='_ATL_SSHClim_'
        plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

        return

################################################################################################################################
def MOC_maps( zlon, zlat, zMy_MOC, zdepth, zCONF, zCASE, zclimyear, zs_year, ze_year, zncout ) :
################################################################################################################################

        plt_AMOCTS=True

        if  plt_AMOCTS: 
        #if not plt_AMOCTS: 
                # AMOC 
                #######

                # Prepare 2 dimnsional (lat,depth) array for ploting   
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
                ztitle=zCASE +' mean AMOC over \n'+zclimyear
                vmin=-15. ; vmax=15. ; vint=1.    ;   contours=npy.arange(vmin,vmax+vint,vint)
                limits=[vmin,vmax,vint]  ;             myticks=npy.arange(vmin,vmax+vint,vint)
                norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])

                num_fram=210
                fram=num_fram+1
                plt.subplot(fram)
                plt.title( ztitle, fontsize=6 )
                plt.contourf( ypltz, zplt*(-1.e-3), npy.squeeze(zMy_MOC), contours, cmap=my_cmap, norm=norm, extend='both' )
                contours = npy.arange(vmin,vmax+vint,2.*vint)
                C = plt.contour( ypltz, zplt*(-1.e-3), npy.squeeze(zMy_MOC), linewidths=0.5, levels=contours, colors='k' )
                plt.clabel( C, C.levels, inline=True, fmt='%3.0f', fontsize=6 )
                #plt.tight_layout()

                zfile_ext='_AMOCClim_'
                plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)


        if plt_AMOCTS :
                # Max. AAMOC Time-series @ 40N or 43N
                #######################################

                time_grid = npy.arange(zs_year,ze_year+1,1.,dtype=int)
                newlocsx  = npy.array(time_grid,'f')
                newlabelsx = npy.array(time_grid,'i')

                jloc= 205 # Lat 40N
                #jloc= 252 # Lat 43N
                if jloc == 205   : zlat='L40N'
                elif jloc == 252 : zlat='L43N'

                plt.clf()
                ax=plt.subplot(211)
                plt.plot( time_grid, zMy_MOC.isel(y=jloc).max(dim='depthw'), 'g',linewidth=1., label=zCONF )
                plt.title( zCASE+' Max AAMOCz @ '+zlat+'\n '+str(zs_year)+str(ze_year), size=9 )
                plt.ylabel( 'Max AMOCz  \n'+r'(Sv)', size=7 )
                plt.ylim([10.,15.])
                plt.xticks(newlocsx,newlabelsx,size=5)
                plt.setp(ax.get_xticklabels(),rotation=90)
                plt.yticks(size=6)
                plt.grid(True, linewidth=0.7,linestyle='--',alpha=0.7,color='grey')

                ax=plt.subplot(212)
                plt.plot( time_grid, zMy_MOC.max().item(), 'k-*',linewidth=1., label=zCONF )
                #plt.plot( time_grid, npy.max(zMy_MOC[:,:,:],axis=(1,2)), 'k', linewidth=1., label=zCONF )
                #plt.title( zCASE+' Max AMOCz '+'\n '+str(zs_year)+str(ze_year), size=9 )
                plt.ylabel( 'Global Max AMOCz  \n'+r'(Sv)', size=7 )
                plt.ylim([10.,20.])
                plt.xticks(newlocsx,newlabelsx,size=5)
                plt.setp(ax.get_xticklabels(),rotation=90)
                plt.yticks(size=6)
                plt.grid(True, linewidth=0.7,linestyle='--',alpha=0.7,color='grey')
                #plt.tight_layout()

                zfile_ext='_MaxAMOCz_'+zlat+'_TiSe_'
                plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

        return

################################################################################################################################
def MTS_maps( zlon, zlat, zCONF, zCASE, zMLD_M, zMLD_S, zMy_varM, zMy_varS, zgdepw_0, ze3t_0, zclimyear, indir, grdir, zncout ) :
################################################################################################################################

	# Mask all levels below the MLD 
	e3t_0msk_SeasM= xr.where( zgdepw_0 < zMLD_M, ze3t_0, 0.)
	e3t_0msk_SeasS= xr.where( zgdepw_0 < zMLD_S, ze3t_0, 0.)

	# Sum all e3t scale factors over the vertcal axis
	e3t_0sum_SeasM = e3t_0msk_SeasM.sum(dim='z').squeeze()
	e3t_0sum_SeasS = e3t_0msk_SeasS.sum(dim='z').squeeze()

	# Compute the T/S mean within the ML in March/September
	T_mldM = (e3t_0msk_SeasM * zMy_varM['votemper']).sum(dim='z').squeeze()/e3t_0sum_SeasM
	S_mldM = (e3t_0msk_SeasM * zMy_varS['vosaline']).sum(dim='z').squeeze()/e3t_0sum_SeasM
                                                       
	T_mldS = (e3t_0msk_SeasS * zMy_varM['votemper']).sum(dim='z').squeeze()/e3t_0sum_SeasS
	S_mldS = (e3t_0msk_SeasS * zMy_varS['vosaline']).sum(dim='z').squeeze()/e3t_0sum_SeasS

	# MLTS from MIMOC observations
	mlT_obs, mlS_obs, lon_obs, lat_obs = MLTS_OBS()

	# Plots Temperature maps 
	########################
	plt.clf()
	plt.subplot(221)
	vmin=-2. ; vmax=6. ; vint=0.2
	contours=npy.arange(vmin,vmax+vint,vint)  
	limits=[vmin,vmax,vint]			 
	myticks=npy.arange(vmin,vmax+vint,vint) 
	
	ztitle=zCASE +' March ML mean T over \n'+str(zclimyear)
	my_cblab=r'($^\circ$C)'
	my_cmap= plt.get_cmap('Spectral_r')
	zoutmap = Iso_Bat( ztype='isol1000' )
	Proj_plot( zlon, zlat, T_mldM, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap )

	plt.subplot(222)
	ztitle=zCASE +' September ML mean T over \n'+str(zclimyear)
	zoutmap = Iso_Bat( ztype='isol1000' )
	Proj_plot(zlon, zlat, T_mldS, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap)

	plt.subplot(223)
	ztitle='MIMOC March ML mean T'
	zoutmap = Iso_Bat( ztype='isol1000' )
	Proj_plot( lon_obs, lat_obs, mlT_obs[2,:,:].squeeze(), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap )

	plt.subplot(224)
	ztitle='MIMOC September ML mean T'
	zoutmap = Iso_Bat( ztype='isol1000' )
	Proj_plot( lon_obs, lat_obs, mlT_obs[8,:,:].squeeze(), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap )

	plt.tight_layout()

	zfile_ext='_MTSClimT_'
	plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zclimyear)+'.png',dpi=300)

	# Plots Salinity maps 
	##################### 
	plt.clf()
	plt.subplot(221)
	vmin=26. ; vmax=36. ; vint=0.5
	contours=npy.arange(vmin,vmax+vint,vint)  
	limits=[vmin,vmax,vint]			 
	myticks=npy.arange(vmin,vmax+2.*vint,2.*vint) 

	ztitle=zCASE +' March ML mean S over \n'+str(zclimyear)
	my_cblab=r'(PSU)'
	my_cmap= plt.get_cmap('Spectral_r')
	zoutmap = Iso_Bat( ztype='isol1000' )
	Proj_plot( zlon, zlat, S_mldM, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar='MLTSS' )

	plt.subplot(222)
	ztitle=zCASE +' September ML mean S over \n'+str(zclimyear)
	zoutmap = Iso_Bat( ztype='isol1000' )
	Proj_plot( zlon, zlat, S_mldS, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar='MLTSS' )

	plt.subplot(223)
	ztitle='MIMOC March ML mean S'
	zoutmap  =Iso_Bat( ztype='isol1000' )
	Proj_plot( lon_obs, lat_obs, mlS_obs[2,:,:].squeeze(), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar='MLTSS' )

	plt.subplot(224)
	ztitle='MIMOC September ML mean S'
	zoutmap = Iso_Bat( ztype='isol1000' )
	Proj_plot( lon_obs, lat_obs, mlS_obs[8,:,:].squeeze(), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar='MLTSS' )

	plt.tight_layout()

	zfile_ext='_MTSClimS_'
	plt.savefig(zCONF+'-'+zCASE+zfile_ext+'y'+str(zclimyear)+'.png',dpi=300)

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
		nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_MTSClim_'+'y'+zclimyear+'.nc'
		ds_out.to_netcdf(nc_f,engine='netcdf4')

	return

################################################################################################################################
def AWT_maps( zlon, zlat, zMy_var1T, zMy_var1S, zdepth, zCONFIG, zCASE, zclimyear, zncout) :
################################################################################################################################

	zCONFCASE=zCONFIG+'-'+zCASE

	######################################################################################
	# Start with model output first 

        # To keep Temp where S > 33.5 only (away from the surface)
	Smask335 = (zMy_var1S > 33.5)
	Temp = zMy_var1T.where(Smask335)
	Temp_filled = (Temp.fillna(-10)).compute()

	# Find the indices over z where T is max.
	depth_map = Temp_filled.argmax(dim='z',skipna=True)
        # Now get the effective depth of T max 
	deptht = zMy_var1T['z']
	zAWTmax_depth1 = deptht[depth_map.compute()]

	# Find the T max value over z 
	temp_map = Temp_filled.max(dim='z',skipna=True)
	mask = (temp_map.compute()>-10)
	zAWTmax1 = temp_map.where(mask,npy.nan)

	######################################################################################
	# Now switch to PHC 3.0 Obs. data

	# Read the associated data
	zMy_varTinit, zMy_varSinit ,lon_obs, lat_obs = PHC3_OBS()
	## Set depths manualy 
	#zdepth_phc3=npy.array([0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, \
	#	     700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2500, \
	#	     3000, 3500, 4000, 4500, 5000, 5500 ])

        # To keep Temp where S > 33.5 only (away from the surface)
	Smask335 = (zMy_varSinit > 33.5)
	TempInit = zMy_varTinit.where(Smask335)
	TempInit_filled = (TempInit.fillna(-10)).compute()

	# Find the indices over z where T is max.
	depthInit_map = TempInit_filled.argmax(dim='depth',skipna=True)
        # Now get the effective depth of T max 
	zAWTmax_depthI = zMy_varTinit.depth[depthInit_map.compute()]
	#zAWTmax_depthI = zdepth_phc3[depthInit_map].compute()]

	# Find the T max value over z 
	TempInit_map = TempInit_filled.max(dim='depth',skipna=True)
	maskInit = (TempInit_map.compute()>-10)
	zAWTmaxI = TempInit_map.where(maskInit,npy.nan)

	# Make the plot for the AW Max Temp 
	#############################################################################################
	vmin=0. ; vmax=7. ; vint=0.5
	contours=npy.arange(vmin,vmax+vint,vint)  
	limits=[vmin,vmax,vint]			 
	myticks=npy.arange(vmin,vmax+vint,vint) 
	
	my_cblab=r'($^\circ$C)'
	my_cmap= plt.get_cmap('jet')

	plt.clf()
	plt.figure()
	plt.subplot(221)
	ztitle=zCASE +' AW Max Temp over \n'+zclimyear
	zoutmap = Iso_Bat( ztype='isol1000' )
	zMyvar = 'votemper'	
	Proj_plot( zlon, zlat, zAWTmax1, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar )

	plt.subplot(222)
	ztitle=' AW Max Temp from \n'+' PHC 3.0'
	zoutmap = Iso_Bat( ztype='isol1000' )
	zMyvar = 'votemper'	
	Proj_plot( lon_obs, lat_obs, zAWTmaxI, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar ) 
	
	# Make the plot for the AW Max Temp depth
	#############################################################################################
	vmin=0. ; vmax=800. ; vint=100.
	contours=npy.arange(vmin,vmax+vint,vint)  
	limits=[vmin,vmax,vint]			 
	myticks=npy.arange(vmin,vmax+vint,vint) 
	
	ztitle=zCASE+' AW Max Temp depth '
	my_cblab=r'(m)'
	my_cmap= plt.get_cmap('jet')

	plt.subplot(223) 
	zoutmap = Iso_Bat( ztype='isol1000' )
	Proj_plot( zlon, zlat, zAWTmax_depth1, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap )

	plt.subplot(224) 
	ztitle=' AW Max Temp depth from \n'+' PHC 3.0'
	zoutmap = Iso_Bat( ztype='isol1000' )
	Proj_plot( lon_obs, lat_obs, zAWTmax_depthI, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap )
	plt.tight_layout()
	
	zfile_ext='_AWTmaxDepth_'
	plt.savefig(zCONFCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

	if zncout:
		ds_out = xr.Dataset()
		
		# FWC field 
		#######################
		ds_out['AWTmax_mod']= (('y','x'), zAWTmax1.values.astype('float32')) 
		ds_out['AWTmax_mod'].attrs['long_name']='Model AWT max calculated using a salinity criteria value Sref=33.5'
		ds_out['AWTmax_mod'].attrs['units']='DegC'
                
		ds_out['AWTmax_init']= (('yobs','xobs'), zAWTmaxI.values.astype('float32')) 
		ds_out['AWTmax_init'].attrs['long_name']='Model Initial AWT max calculated using a salinity criteria value Sref=33.5'
		ds_out['AWTmax_init'].attrs['units']='DegC'
                
		ds_out['AWTmaxDepth_mod']= (('y','x'), zAWTmax_depth1.values.astype('float32')) 
		ds_out['AWTmaxDepth_mod'].attrs['long_name']='Model AWT max depth calculated using a salinity criteria value Sref=33.5'
		ds_out['AWTmaxDepth_mod'].attrs['units']='m'
                
		ds_out['AWTmaxDepth_init']= (('yobs','xobs'), zAWTmax_depthI.values.astype('float32')) 
		ds_out['AWTmaxDepth_init'].attrs['long_name']='Model AWT max depth calculated using a salinity criteria value Sref=33.5'
		ds_out['AWTmaxDepth_init'].attrs['units']='m'
                
		ds_out['lat_obs']= (('yobs','xobs'), lat_obs.astype('float32')) 
		ds_out['lat_obs'].attrs['long_name']='Degrees north'
		ds_out['lat_obs'].attrs['units']='Deg'
		
		ds_out['lon_obs']= (('yobs','xobs'), lon_obs.astype('float32')) 
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
		nc_f = './NETCDF/'+zCONFIG+'-'+zCASE+'_AWTClim_'+'y'+zclimyear+'.nc'
		ds_out.to_netcdf(nc_f,engine='netcdf4')

	return

################################################################################################################################
def FWC_maps( zlon, zlat, zMy_var1S, zMy_varSinit, zMy_var1ssh, zCONFIG, zCASE, zclimyear, ze3, ztmask, zncout ) :
################################################################################################################################

	# FWC calculation over the year
	###########################################
	Sref=34.80*1.004715
	#Sref=34.80 
	 
	print('				FWC calculation & plot ')

	# Freshwater content from model outputs 
	######################################################################################
	fwcmask = xr.where( zMy_var1S > Sref, 0., ze3 ) 
	FW4D = (Sref - zMy_var1S) / Sref * fwcmask
	# Sum over depth 
	fwc2D = FW4D.sum(dim="z")
	# Mask land area 
	fwc2D = xr.where( ztmask[0,:,:] < 1, npy.nan, fwc2D )

	# Freshwater content from initial state 
	######################################################################################
	fwcmask_init = xr.where( zMy_varSinit > Sref, 0., ze3 )
	FW4D_init = (Sref - zMy_varSinit) / Sref * fwcmask_init
	# Sum over depth 
	fwc2D_init = FW4D_init.sum(dim="z")
	# Mask land area 
	fwc2D_init = xr.where( ztmask[0,:,:] < 1, npy.nan, fwc2D_init )

	dbg=False
	if zCONFIG == 'CREG025.L75' :
		idbg=185 ; jdbg=515 ; kdbg=50	# CREG025.L75
	elif zCONFIG == 'CREG12.L75' :
		idbg=562 ; jdbg=1558 ; kdbg=50	 # CREG12.L75
	
	if dbg:
		print()
		print(' Print a specific point to debug within the Beaufort Gyre')
		print('   		ze33Dtime_msk[ti,0:kdbg,jdbg,idbg]: '+str(ze33Dtime_msk[0,0:kdbg,jdbg,idbg])+' '+str(zMy_var1S[0,0:kdbg,jdbg,idbg]))
		print('   		fwc2D[jdbg,idbg]: '+fwc2D[jdbg,idbg])

	# Use obs. from Proshutinsky GRL2018
	locpath='./DATA/'
	locfile='BeaufortGyreFWC-Obs-Proshutinsky_GRL2018_y2003-2017.nc'
	if chkfile(locpath+locfile,zstop=True) :
		ds_obsfwc = xr.open_dataset(locpath+locfile)
		fwc2D_obs = ds_obsfwc['BGFWC']
		lat2D_obs = ds_obsfwc['lat']
		lon2D_obs = ds_obsfwc['lon']

		lat2D_obs = xr.where( lat2D_obs == 9999., npy.nan, lat2D_obs ) 
		lon2D_obs = xr.where( lon2D_obs == 9999., npy.nan, lon2D_obs ) 

	if int(zclimyear) >= 2003 and int(zclimyear) <= 2017 :
		s_ind=int(zclimyear)-2003  ; e_ind=s_ind+1
		mean_FWCObs = fwc2D_obs[s_ind:e_ind,:,:].squeeze()
		mean_FWCObs = xr.where( mean_FWCObs == 9999., npy.nan, mean_FWCObs )
		obsper=zclimyear
	else:
		mean_FWCObs = fwc2D_obs.mean(dim='t').squeeze()
		mean_FWCObs = fwc2D_obs.mean(dim='t').squeeze()
		mean_FWCObs = xr.where( mean_FWCObs == 9999., npy.nan, mean_FWCObs )
		obsper='2003-2017'

	# Get observations SSH
	obs_ssh, lon_obs, lat_obs, obs_ssh_per = SSH_OBS( t_year=int(zclimyear[0:4]) )
	obs_ssh = xr.where( obs_ssh >= 9e20, npy.nan, obs_ssh )

	plt.clf()
	fig=plt.figure()

	# Plot the FWC map mean over the year
	#####################################
	my_cblab=r'(m)'
	ztitle=zCASE+' mean FWC (m) over \n'+zclimyear
	vmin=0. ; vmax=25. ; vint=2.
	contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	limits=[vmin,vmax,vint]			       # limits for eke
	myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)
	my_cmap=plt.get_cmap('Spectral_r')
	
	plt.subplot(231)
	zMyvar='FWC'
	zoutmap = Iso_Bat( ztype='isol1000' )
	Proj_plot( zlon, zlat, fwc2D, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar )

	plt.subplot(232)
	maxval=npy.nanmax(mean_FWCObs)
	fig.text(0.45,0.78,'Max: '+str(maxval)+' m',fontsize=5,color='r')
	ztitle=' Mean FWC (m) from \n'+' BG Obs Sys. (Proshutinsky et al. GRL2018) \n '+ ' over year ' + obsper
	zoutmap = Iso_Bat( ztype='isol1000' )
	Proj_plot( lon2D_obs, lat2D_obs, mean_FWCObs, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar )

	plt.subplot(233)
	ztitle=' Mean FWC (m) from \n'+' Init State '
	zoutmap = Iso_Bat( ztype='isol1000' )
	Proj_plot( zlon, zlat, fwc2D_init, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar )
	
	plt.subplot(234)
	zMyvar='ssh'
	seas=''
	contours, limits, myticks, ztitle, zfile_ext, my_cblab, my_cmap, m_alpha = SET_ARC_CNT( zCASE, zclimyear, seas, zMyvar )
	zoutmap = Iso_Bat( ztype='isol1000' )
	Proj_plot( zlon, zlat, zMy_var1ssh*m_alpha, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar )

	plt.subplot(235)
	zMyvar='ssh'
	seas=''
	contours, limits, myticks, ztitle, zfile_ext, my_cblab, my_cmap, m_alpha = SET_ARC_CNT( zCASE, zclimyear, seas, zMyvar, zplot_obs=1 )
	zoutmap = Iso_Bat( ztype='isol1000' )
	Proj_plot( lon_obs, lat_obs, obs_ssh*m_alpha, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar )
	#plt.tight_layout()

	zfile_ext='_FWCSSHClim_'
	plt.savefig(zCONFIG+'-'+zCASE+zfile_ext+'y'+zclimyear+'.png',dpi=300)

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

		ds_out['lat_obs']= (('yobs','xobs'), lat2D_obs.values.astype('float32')) 
		ds_out['lat_obs'].attrs['long_name']='Degrees north'
		ds_out['lat_obs'].attrs['units']='Deg'
		
		ds_out['lon_obs']= (('yobs','xobs'), lon2D_obs.values.astype('float32')) 
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
		nc_f = './NETCDF/'+zCONFIG+'-'+zCASE+'_FWCClim_'+'y'+zclimyear+'.nc'
		ds_out.to_netcdf(nc_f,engine='netcdf4')

		ds_out = xr.Dataset()

		# SSH field 
		#######################
		zMy_var1ssh = xr.where( ztmask[0,:,:] < 1, npy.nan, zMy_var1ssh )
		ds_out['ssh_mod']= (('y','x'), zMy_var1ssh.values.astype('float32')) 
		ds_out['ssh_mod'].attrs['long_name']='Model SSH '
		ds_out['ssh_mod'].attrs['units']='m'
                
		ds_out['ssh_obs']= (('yobs','xobs'), obs_ssh.values.astype('float32')) 
		ds_out['ssh_obs'].attrs['long_name']='Observed DOT from Armitage et al. 2017 for the period: '+obs_ssh_per
		ds_out['ssh_obs'].attrs['units']='m'
                
		ds_out['lat_obs']= (('yobs','xobs'), lat_obs.astype('float32')) 
		ds_out['lat_obs'].attrs['long_name']='Degrees north'
		ds_out['lat_obs'].attrs['units']='Deg'
		
		ds_out['lon_obs']= (('yobs','xobs'), lon_obs.astype('float32')) 
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
		nc_f = './NETCDF/'+zCONFIG+'-'+zCASE+'_SSHClim_'+'y'+zclimyear+'.nc'
		ds_out.to_netcdf(nc_f,engine='netcdf4')

	return

################################################################################################################################
def CREG_MSK( zCASE ) :
################################################################################################################################

        # Read the CREG025.L75 mask 
        locpath='./CREG025.L75/GRID/'
        locfile='CREG025.L75-'+zCASE+'_mask.nc'
        if chkfile(locpath+locfile) :
                ds_msk=xr.open_dataset(locpath+locfile)[['glamt','gphit','tmask']]
                lon = ds_msk['glamt'].squeeze()
                lat = ds_msk['gphit'].squeeze()
                mask= ds_msk['tmask'].squeeze()

        return mask, lon, lat

################################################################################################################################
def CREG_INIT( zCONFIG, zCASE ) :
################################################################################################################################
	# Read initial state to compare with
	print('                      Read initial state  ')
	locpath=zCONFIG+'/'+zCONFIG+'-'+zCASE+'-MEAN/'
	locfile=zCONFIG+'-'+zCASE+'_init_gridT.nc'
	if chkfile(locpath+locfile) : 
		ds_TSinit = xr.open_dataset(locpath+locfile, engine="netcdf4")[['votemper','vosaline']]
		ds_TSinit = ds_TSinit.rename({'nav_lev':'z'})
		My_varTinit = ds_TSinit['votemper'].squeeze()
		My_varSinit = ds_TSinit['vosaline'].squeeze()
	
	return My_varTinit, My_varSinit

################################################################################################################################
def Iso_Bat( ztype='isol1000', zarea='arctic' ) :
################################################################################################################################

	locpath='./'
	locfile='Bathymetry.nc'
	if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
		ds_bat = xr.open_dataset(locpath+locfile)
		lon = ds_bat['nav_lon'].squeeze()
		lat = ds_bat['nav_lat'].squeeze()
		My_var = ds_bat['bathy_meter'].squeeze()
	
	spval = 0.
	My_var = xr.where( My_var <= spval, npy.nan, My_var )
	
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
		m = Basemap(projection='npstere',boundinglat=65,lon_0=-60, resolution='i')
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
		m = Basemap(llcrnrlon=-180,llcrnrlat=66,urcrnrlon=-80,urcrnrlat=80, resolution='i',\
		            projection='cass',lon_0=-140,lat_0=60)    
	############################################################################################################
	elif zarea == 'cassis_BGZoom_HR' :
		m = Basemap(llcrnrlon=-80,llcrnrlat=80,urcrnrlon=-180,urcrnrlat=60, resolution='i',\
		            projection='cass',lon_0=0,lat_0=80)    
	############################################################################################################
	else: # Focus on North Atlantic sector
		m = Basemap(width=6100000,height=5000000,lat_1=30.,lat_2=70,lon_0=-45,lat_0=45,projection='aea',resolution='i')
		zcolorbat='grey'   ;  zalpha=0.7

	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	pal = plt.get_cmap('binary')
	X,Y = m(lon.values,lat.values)

	# contour (optional)
	CS2 = m.contour( X, Y, My_var.values, linewidths=0.5,levels=contours, colors=zcolorbat, alpha=zalpha )
	plt.clabel(CS2, CS2.levels, inline=True, fmt='%.0f', fontsize=3)

	return m, X, Y

################################################################################################################################
def Proj_plot( lon, lat, tab, contours, limits, myticks=None, name=None, zmy_cblab=None, zmy_cmap=None, filename='test.png', zvar=None, ztickslabels=None, zarea='arctic', data_ref=False ) :
################################################################################################################################
	#
	plt.rcParams['text.usetex']=False
	plt.rcParams['font.family']='serif'
	plt.rcParams['axes.unicode_minus'] = False
	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	############################################################################################################
	if zarea == 'arctic': # Focus on Arctic basin
		m = Basemap(projection='npstere',boundinglat=65,lon_0=-60, resolution='i')
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
		m = Basemap(llcrnrlon=-180,llcrnrlat=66,urcrnrlon=-80,urcrnrlat=80, resolution='i',\
		            projection='cass',lon_0=-140,lat_0=60)    
	############################################################################################################
	elif zarea == 'cassis_BGZoom_HR' :
		m = Basemap(llcrnrlon=-80,llcrnrlat=80,urcrnrlon=-180,urcrnrlat=60, resolution='i',\
		            projection='cass',lon_0=0,lat_0=80)    
	############################################################################################################
	elif zarea == 'ginsea': # Focus on GIN Seas
		m = Basemap(width=1400000,height=1600000,lat_1=50.,lat_2=65,lon_0=0,lat_0=74.,projection='aea',resolution='i')
	############################################################################################################
	else: # Focus on North Atlantic sector
		 m = Basemap(width=6100000,height=5000000,lat_1=30.,lat_2=70,lon_0=-45,lat_0=45,projection='aea',resolution='i')
	
	if zvar == 'sivolu' or zvar == 'siconc'  or zvar == 'ssh' or zvar == 'FWC' or zvar == 'voeke' or zvar == 'sobarstf' :
		zfontsize=4.
	else:
		zfontsize=6.
	
	if zvar != 'Bathy' :
		if zarea == 'arctic': 
			m.drawparallels(npy.arange(-90.,91.,5.),labels=[False,False,False,False], size=zfontsize, linewidth=0.3)
			m.drawmeridians(npy.arange(-180.,181.,20.),labels=[True,False,False,True], size=zfontsize, latmax=90.,linewidth=0.3)
		else :
			m.drawparallels(npy.arange(-90.,91.,2.),labels=[True,False,False,False], size=zfontsize, linewidth=0.3, color='grey',alpha=0.70 )
			m.drawmeridians(npy.arange(-180.,181.,5.),labels=[False,False,False,True], size=zfontsize, latmax=90.,linewidth=0.3, color='grey',alpha=0.70 )
		m.fillcontinents(color='grey',lake_color='white')
	
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	
	if zmy_cmap != None :
		pal = zmy_cmap
	else:
		pal = plt.get_cmap('coolwarm')
	
	X,Y = m(lon,lat)
	
	if zarea == 'GulfS' and zvar == 'votemper' :
		if data_ref :
			zlinewidths=1.1   ; zcolor='g'
		else:	
			zlinewidths=0.8   ; zcolor='r'
		C = m.contour( X,Y,tab,linewidths=zlinewidths,levels=[17.], colors=zcolor )
	else:
		C = m.contourf( X,Y,tab,contours,cmap=pal,norm=norm,extend='both' )
		if zvar == 'ssh' :
			CS=m.contour(X, Y, tab, linewidths=0.5, levels=npy.arange(limits[0],limits[1],5.), colors='k', alpha=0.4)

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
				x,y = m(lons,lats)
				m.scatter(x,y,1,marker='o', color='r')
				#m.plot(x,y,linewidth=2, color='g')
		############################################################################################################
		############################################################################################################
	
		# colorbar	
		if myticks is None:
			cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
		else:
			if zvar == 'votemper' or zvar == 'vosaline' or zvar == 'sivolu' or zvar == 'sobarstf' : 
				cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8,drawedges=True)
			elif zvar == 'voeke' :
				cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.6,drawedges=True)
			elif zvar == 'MLTSS' :
				cbar = plt.colorbar(C,ticks=myticks,format='%.0f',orientation='vertical',shrink=0.8,drawedges=True)
			else:
				cbar = plt.colorbar(C,format='%.0f',orientation='vertical',shrink=0.8,drawedges=True)
	
			cbar.set_label(zmy_cblab,fontsize=zfontsize)
			cl = plt.getp(cbar.ax, 'ymajorticklabels')
			plt.setp(cl, fontsize=zfontsize)
			#if ztickslabels != None and zvar == 'voeke' : 
			#	zticks = npy.linspace(1e-6, 1e-2, 5)
			#	cbar.set_ticks(10**zticks)
			#	cbar.ax.set_yticklabels(ztickslabels)
	
	plt.title(name,fontsize=zfontsize)
	
	return m

################################################################################################################################
def simple_maps( zlon, zlat, zCONF, zCASE, zMy_var1, zMyvar, zclimyear, slev=None, seas='', zfram=111, plot_obs=0, ano=0 ) :
################################################################################################################################

	m_alpha=1.
	
	# Do the plot 
	print() 
	print('                    plot '+zMyvar+' field')
	print() 
	
	plt.subplot(zfram)
	contours, limits, myticks, ztitle, zfile_ext, my_cblab, my_cmap, m_alpha, mytickslabels = SET_ARC_CNT( zCASE, zclimyear, seas, zMyvar, zslev=slev, zplot_obs=plot_obs, zdiff=ano )
	zoutmap = Iso_Bat( ztype='isol1000' )
	m = Proj_plot( zlon, zlat, zMy_var1[:,:]*m_alpha, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, ztickslabels=mytickslabels )

	return m 

################################################################################################################################
def BFG_compute( lon, lat, ssh_raw, depth, var_type, increment, grid_area, rm_landbarrier=0 ) :
################################################################################################################################

	# History: This algorithme has been originaly developed by Heather Regan and slightly adapeted to be included into the MONARC
	# 	   See Regan et al. JPO2020 ; https://doi.org/10.1175/JPO-D-19-0234.1

	## This function computes the largest closed contour in the Western Arctic basin
	###### INPUTS:
	##        lon:            longitude
	##        lat:            latitude
	##        ssh_raw:        the ssh field we're examining
	##        depth:          the bathymetry (depth > 0) to check off-shelf regions 
	##        var_type:       usually set to "SSH": identifies the variable type (has also been used for MSL in past)
	##        increment:      the increment with which to iterate out from the maximum. Usually start at 10cm. 
	##                        Higher resolution needs a smaller increment than lower resolution, because field varies more and so 
	##                        larger increment may miss small features. But smaller increment takes longer
	##        rm_landbarrier  whether or not to use coastline as a valid edge of contour (e.g. set to 1 if using MSL, as atmospheric variable)
	###### OUTPUTS
	##        mask_full:      the identified closed contour
	##        BG_max_val:     Max. SSH value 
	##        BGcalcmin:      Min. SSH value
	##        lat:            latitude
	##        lon:            longitude
	##        BG_area:        Surface area of the closed contour 
	
	## Step 1: set up coastline to determine when the contour is no longer closed
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
	#This ensures that the maximum nonzero value off the shelf is found. Otherwise can get a high maxima near the coast
	## Artificially force by depth field if rm_landbarrier isn't there 
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
	## For each new increment, check that all cells in this new contour are a) not adjacent to land, and b) not higher than previous maximum
	  
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
	    checkarr = fn_getEdge(maskarr_new);
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
	    
	    new_edges = fn_getEdge(maskarr_new)+maskarr_new;
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
def fn_getEdge( oldarr ) :
################################################################################################################################

	# History: This algorithm has been originaly developed by Heather Regan and slightly adapeted to be included into the MONARC
	# 	   See Regan et al. JPO2020 ; https://doi.org/10.1175/JPO-D-19-0234.1

	## This function identifies the edge of a contour by looking at the four adjacent cells
	
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
def EKE_compute( zlon, zlat, zCONF, zCASE, xiosfreq, zc_year, zdatadir, zncout ) :
################################################################################################################################

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
	ds_eke = ds_eke.assign_coords( z=nemoT.z, y=nemoT.y, x=nemoT.x, glamt=zlon, gphit=zlat )
	ds_eke['voeke'] = (('z','y','x'), EKE_mmT.mean(dim='t').values)

	return ds_eke 
