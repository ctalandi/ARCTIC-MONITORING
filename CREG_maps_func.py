#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import sys
import numpy as npy
from CREG_maps_cont import *
from checkfile import *
import subprocess
import xarray as xr 
from datetime import datetime

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

def MTS_maps( zlon, zlat, zCONF, zCASE, zMLD_M, zMLD_S, zMy_varM, zMy_varS, zgdepw_0, ze3t_0, zclimyear, indir, grdir, zncout ) :

	#e3t_0msk_SeasM=npy.zeros((ze3t_0.shape[0],zlon.shape[0],zlon.shape[1]))
	#e3t_0msk_SeasS=npy.zeros((ze3t_0.shape[0],zlon.shape[0],zlon.shape[1]))
	#e3t_0sum_SeasM=npy.zeros((zlon.shape[0],zlon.shape[1]))
	#e3t_0sum_SeasS=npy.zeros((zlon.shape[0],zlon.shape[1]))
	#T_mldM=npy.zeros((2,zlon.shape[0],zlon.shape[1]))
	#S_mldM=npy.zeros((2,zlon.shape[0],zlon.shape[1]))

	# Mask all levels below the MLD 
	e3t_0msk_SeasM= xr.where( zgdepw_0 < zMLD_M, ze3t_0, 0.)
	e3t_0msk_SeasS= xr.where( zgdepw_0 < zMLD_S, ze3t_0, 0.)
	#e3t_0msk_SeasM= npy.squeeze(npy.where(zgdept_0[:,:,:] < zMy_var1SeasM[:,:], ze3t_0[:,:,:], 0.))
	#e3t_0msk_SeasS= npy.squeeze(npy.where(zgdept_0[:,:,:] < zMy_var1SeasS[:,:], ze3t_0[:,:,:], 0.))

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
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, T_mldM, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap)

	plt.subplot(222)
	ztitle=zCASE +' September ML mean T over \n'+str(zclimyear)
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, T_mldS, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap)

	plt.subplot(223)
	ztitle='MIMOC March ML mean T'
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(lon_obs, lat_obs, mlT_obs[2,:,:].squeeze(), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap)

	plt.subplot(224)
	ztitle='MIMOC September ML mean T'
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(lon_obs, lat_obs, mlT_obs[8,:,:].squeeze(), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap)

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
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, S_mldM, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar='MLTSS')

	plt.subplot(222)
	ztitle=zCASE +' September ML mean S over \n'+str(zclimyear)
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, S_mldS, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar='MLTSS')

	plt.subplot(223)
	ztitle='MIMOC March ML mean S'
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(lon_obs, lat_obs, mlS_obs[2,:,:].squeeze(), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar='MLTSS')

	plt.subplot(224)
	ztitle='MIMOC September ML mean S'
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(lon_obs, lat_obs, mlS_obs[8,:,:].squeeze(), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar='MLTSS')

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

def simple_maps( zlon, zlat, zCONF, zCASE, zMy_var1, zMyvar, zclimyear, slev=None, seas='', zfram=111, plot_obs=0, ano=0) :

	m_alpha=1.
	
	# Do the plot 
	print(	zMyvar+' plot')
	print() 
	
	plt.subplot(zfram)
	contours,limits,myticks,ztitle,zfile_ext,my_cblab,my_cmap,m_alpha = SET_ARC_CNT(zCASE,zclimyear,seas,zMyvar,zslev=slev,zplot_obs=plot_obs,zdiff=ano)
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, zMy_var1[:,:]*m_alpha, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar)

	return

def AWTmax_maps( zlon, zlat, zMy_var1T, zMy_var1S, zdepth, zMyvar, zCONFIG, zCASE, zclimyear, zncout) :

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

# 	# Start to mask field where Salinity is lower than 33.5 PSU
# 	# Necessary to remove the surface temperature maxima than can arise
# 	zMy_var1=npy.ma.masked_where((zMy_var1S <= 33.5),zMy_var1)
# 	zMy_varTinit=npy.ma.masked_where((zMy_varSinit <= 33.5),zMy_varTinit)
# 
# 	# Find the Max temp. over depth
# 	zAWTmax1=npy.amax(zMy_var1[:,:,:],axis=0).squeeze()
# 	zAWTmaxI=npy.amax(zMy_varTinit[:,:,:],axis=0).squeeze()
# 
# 	# Find the depth of the Max temp. 
# 	zAWTmax_depth1_ind=npy.argmax(zMy_var1[:,:,:],axis=0).squeeze()
# 	zAWTmax_depthI_ind=npy.argmax(zMy_varTinit[:,:,:],axis=0).squeeze()
# 
# 	zAWTmax_depth1=npy.zeros((zAWTmax_depth1_ind.shape[0],zAWTmax_depth1_ind.shape[1]))
# 	zAWTmax_depthI=npy.zeros((zAWTmax_depthI_ind.shape[0],zAWTmax_depthI_ind.shape[1]))
# 
# 	zdepth_phc3=npy.array([0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, \
# 		     700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2500, \
# 		     3000, 3500, 4000, 4500, 5000, 5500 ])
# 
# 	for jj in set(npy.arange(zAWTmax_depth1_ind.shape[0])) :
# 	    for ji in set(npy.arange(zAWTmax_depth1_ind.shape[1])) :
# 		    zAWTmax_depth1[jj,ji]=zdepth[zAWTmax_depth1_ind[jj,ji]]
# 
# 	for jj in set(npy.arange(zAWTmax_depthI.shape[0])) :
# 	    for ji in set(npy.arange(zAWTmax_depthI.shape[1])) :
# 		    zAWTmax_depthI[jj,ji]=zdepth_phc3[zAWTmax_depthI_ind[jj,ji]]
# 
# 	zAWTmax_depth1=npy.ma.masked_where((zMy_var1S[0,:,:].squeeze() == 0),zAWTmax_depth1)
# 	zAWTmax_depthI=npy.ma.masked_where((zMy_varSinit[0,:,:].squeeze() == 0),zAWTmax_depthI)


	# Make the plot for the AW Max Temp 
	#############################################################################################
	vmin=0. ; vmax=7. ; vint=0.5
	contours=npy.arange(vmin,vmax+vint,vint)  
	limits=[vmin,vmax,vint]			 
	myticks=npy.arange(vmin,vmax+vint,vint) 
	
	ztitle=zCASE +' AW Max Temp over \n'+zclimyear
	my_cblab=r'($^\circ$C)'
	my_cmap= plt.get_cmap('jet')


	plt.clf()
	plt.figure()
	plt.subplot(221)
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, zAWTmax1, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar)

	plt.subplot(222)
	ztitle=' AW Max Temp from \n'+' PHC 3.0'
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(lon_obs, lat_obs, zAWTmaxI, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar)
	
	# Make the plot for the AW Max Temp depth
	#############################################################################################
	zfile_ext='_AWTmaxDepth_'
	
	vmin=0. ; vmax=800. ; vint=100.
	contours=npy.arange(vmin,vmax+vint,vint)  
	limits=[vmin,vmax,vint]			 
	myticks=npy.arange(vmin,vmax+vint,vint) 
	
	ztitle=zCASE+' AW Max Temp depth '
	my_cblab=r'(m)'
	my_cmap= plt.get_cmap('jet')

	plt.subplot(223) 
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, zAWTmax_depth1, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap)

	plt.subplot(224) 
	ztitle=' AW Max Temp depth from \n'+' PHC 3.0'
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(lon_obs, lat_obs, zAWTmax_depthI, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap)
	
	plt.tight_layout()
	
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
		nc_f = './NETCDF/'+zCONFIG+'-'+zCASE+'_AWTmaxClim_'+'y'+zclimyear+'.nc'
		ds_out.to_netcdf(nc_f,engine='netcdf4')

	return

def FWC_maps( zlon, zlat, zMy_var1S, zMy_varSinit, zMy_var1ssh, zCONFIG, zCASE, zclimyear, ze3, ztmask, ztime_dim, zncout) :

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

#	# Use the initial state as reference
#	ze33Dtime=ze3[:,:,:]*ztmask[:,:,:]
#	ze33Dtime_msk = ze33Dtime.copy()
#	ze33Dtime_mskI = ze33Dtime.copy()
#	ze33Dtime_msk[npy.where(zMy_var1S > Sref )] = 0.e0
#	ze33Dtime_mskI[npy.where(zMy_varSinit > Sref )] = 0.e0
#	fwc2D = npy.zeros((zlon.shape[0],zlon.shape[1]))    ;	 fwc2D_init = npy.zeros((zlon.shape[0],zlon.shape[1]))
#	
#	zMy_varSinit=CREG_INIT(zCONFIG,zCASE,zlon,ze3)
#
	#print(' a specific point to debug within the Beaufoirt Gyre')
	dbg=False
	if zCONFIG == 'CREG025.L75' :
		idbg=185 ; jdbg=515 ; kdbg=50	# CREG025.L75
	elif zCONFIG == 'CREG12.L75' :
		idbg=562 ; jdbg=1558 ; kdbg=50	 # CREG12.L75
	
#	for jk in set(npy.arange(75)) :
#	    fwc2D[:,:] = fwc2D[:,:] + ( Sref - zMy_var1S[jk,:,:] ) / Sref * ze33Dtime_msk[jk,:,:] 
#	    fwc2D_init[:,:] = fwc2D_init[:,:] + ( Sref - zMy_varSinit[jk,:,:] ) / Sref * ze33Dtime_mskI[jk,:,:] 
	
	if dbg:
		print('  ze33Dtime_msk[ti,0:kdbg,jdbg,idbg]: '+str(ze33Dtime_msk[0,0:kdbg,jdbg,idbg])+' '+str(zMy_var1S[0,0:kdbg,jdbg,idbg]))
		print('  fwc2D[jdbg,idbg]: '+fwc2D[jdbg,idbg])
#	
#	tmask2D=npy.squeeze(ztmask[0,:,:].copy())
#	fwc2D=npy.ma.masked_where((tmask2D == 0),fwc2D)
#	fwc2D_init=npy.ma.masked_where((tmask2D == 0),fwc2D_init)

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

#		fwc2D_obs=npy.squeeze(fid.variables['BGFWC'])
#		lat2D_obs=npy.squeeze(fid.variables['lat'])
#		lon2D_obs=npy.squeeze(fid.variables['lon'])
#
#		lat2D_obs=npy.ma.masked_where(lat2D_obs == 9999.,lat2D_obs)
#		lon2D_obs=npy.ma.masked_where(lon2D_obs == 9999.,lon2D_obs)

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

#	if int(zclimyear) >= 2003 and int(zclimyear) <= 2017 :
#		s_ind=int(zclimyear)-2003  ; e_ind=s_ind+1
#		mean_FWCObs=npy.squeeze(npy.array(fwc2D_obs[s_ind:e_ind,:,:]).copy())
#		mean_FWCObs=npy.ma.masked_where(mean_FWCObs == 9999.,mean_FWCObs)
#		obsper=zclimyear
#	else:
#		mean_FWCObs=npy.squeeze(npy.nanmean(npy.array(fwc2D_obs),axis=0))
#		mean_FWCObs=npy.ma.masked_where(mean_FWCObs == 9999.,mean_FWCObs)
#		obsper='2003-2017'
	
	# Get observations SSH
	obs_ssh,lon_obs,lat_obs,obs_ssh_per=SSH_OBS(t_year=int(zclimyear[0:4]))
	obs_ssh = xr.where( obs_ssh >= 9e20, npy.nan, obs_ssh )
#	#obs_ssh=npy.ma.masked_where(obs_ssh >= 9e20, obs_ssh)

	plt.clf()
	fig=plt.figure()

	# Plot the FWC map mean over the year
	#####################################
	my_cblab=r'(m)'
	zfile_ext='_FWCClim_'
	ztitle=zCASE+' mean FWC (m) over \n'+zclimyear
	vmin=0. ; vmax=25. ; vint=2.
	contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	limits=[vmin,vmax,vint]			       # limits for eke
	myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)
	my_cmap=plt.get_cmap('coolwarm')
	my_cmap=plt.get_cmap('Spectral_r')
	zMyvar='FWC'
	
	plt.subplot(231)
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, fwc2D, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar)

	plt.subplot(232)
	maxval=npy.nanmax(mean_FWCObs)
	fig.text(0.45,0.78,'Max: '+str(maxval)+' m',fontsize=5,color='r')
	ztitle=' Mean FWC (m) from \n'+' BG Obs Sys. (Proshutinsky et al. GRL2018) \n '+ ' over year ' + obsper
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(lon2D_obs, lat2D_obs, mean_FWCObs, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar)

	plt.subplot(233)
	ztitle=' Mean FWC (m) from \n'+' Init State '
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, fwc2D_init, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar)
	
	plt.subplot(234)
	zMyvar='ssh'
	seas=''
	contours,limits,myticks,ztitle,zfile_ext,my_cblab,my_cmap,m_alpha = SET_ARC_CNT(zCASE,zclimyear,seas,zMyvar)
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, zMy_var1ssh*m_alpha, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar)

	plt.subplot(235)
	zMyvar='ssh'
	seas=''
	contours,limits,myticks,ztitle,zfile_ext,my_cblab,my_cmap,m_alpha = SET_ARC_CNT(zCASE,zclimyear,seas,zMyvar,zplot_obs=1)
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(lon_obs, lat_obs, obs_ssh*m_alpha, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar)

	zfile_ext='_FWCSSHClim_'
	#plt.tight_layout()
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


def Arc_plot(lon,lat,tab,contours,limits,myticks=None,name=None,zmy_cblab=None,zmy_cmap=None,filename='test.png',zvar=None):
	#
	plt.rcParams['text.usetex']=False
	plt.rcParams['font.family']='serif'
	plt.rcParams['axes.unicode_minus'] = False
	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	m = Basemap(projection='npstere',boundinglat=65,lon_0=-60, resolution='i')

	if zvar == 'sivolu' or zvar == 'siconc'  or zvar == 'ssh' or zvar == 'FWC' :
		zfontsize=4.
	else:
		zfontsize=6.

	if zvar != 'Bathy' :
		m.drawparallels(npy.arange(-90.,91.,5.),labels=[False,False,False,False], size=zfontsize, linewidth=0.3)
		m.drawmeridians(npy.arange(-180.,181.,20.),labels=[True,False,False,True], size=zfontsize, latmax=90.,linewidth=0.3)
		m.fillcontinents(color='grey',lake_color='white')

	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])

	if zmy_cmap != None :
		pal = zmy_cmap
	else:
		pal = plt.get_cmap('coolwarm')

	X,Y = m(lon,lat)
	C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm,extend='both')
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
		elif zvar == 'MLTSS' :
			cbar = plt.colorbar(C,ticks=myticks,format='%.0f',orientation='vertical',shrink=0.8,drawedges=True)
		else:
			cbar = plt.colorbar(C,format='%.0f',orientation='vertical',shrink=0.8,drawedges=True)
		
		cbar.set_label(zmy_cblab,fontsize=zfontsize)
		cl = plt.getp(cbar.ax, 'ymajorticklabels')
		plt.setp(cl, fontsize=zfontsize)

	plt.title(name,fontsize=zfontsize)

	return 

def Arc_Bat(ztype='isol1000') :

	locpath='./'
	locfile='Bathymetry.nc'
	if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
		ds_bat = xr.open_dataset(locpath+locfile)
		lon = ds_bat['nav_lon'].squeeze()
		lat = ds_bat['nav_lat'].squeeze()
		My_var = ds_bat['bathy_meter'].squeeze()
	
	spval = 0.
	My_var = xr.where( My_var <= spval, npy.nan, My_var )
	#My_var = npy.ma.masked_where(My_var <= spval,My_var)
	
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
	m = Basemap(projection='npstere',boundinglat=65,lon_0=-60, resolution='i')
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	pal = plt.get_cmap('binary')
	X,Y = m(lon,lat)

	# contour (optional)
	CS2 = m.contour(X, Y, My_var, linewidths=0.5,levels=contours, colors='grey', alpha=0.4)
	##plt.clabel(CS2, CS2.levels, inline=True, fmt='%.0f', fontsize=3)

	return m, X, Y

def SSH_OBS(t_year=1959):

	locpath='./DATA/'
	locfile='EKE_DOT_based_2003-2014.nc'
	if chkfile(locpath+locfile,zstop=True) :
		ds_sshobs = xr.open_dataset(locpath+locfile)
		lon = ds_sshobs['lon']
		lat = ds_sshobs['lat']
		ssh_init = ds_sshobs['DOT']

		SSH_lon2D=npy.tile(lon,(lat.size,1))
		SSH_lat2D=npy.tile(lat,(lon.size,1)).T

	if t_year >= 2003 and t_year <=2014 :
		# Get the specific year
		s_ind=(t_year-2003)*12	 
		out_ssh_OBS = ssh_init.isel(date=slice(s_ind,s_ind+12)).mean(dim='date').squeeze()
		#out_ssh_OBS = npy.mean(ssh_init[s_ind:s_ind+12,:,:],axis=0).squeeze()
		ssh_OBS_obsper=str(t_year)
	else:
		# Compute the mean over the obs. monthly period 2003-2014
		out_ssh_OBS = ssh_init.mean(dim='date').squeeze()
		ssh_OBS_obsper='2003-2014'

	# Remove the domain mean to get an anomaly
	out_ssh_OBS = out_ssh_OBS - ssh_init.mean()
	#out_ssh_OBS = out_ssh_OBS - npy.nanmean(ssh_init)

	return out_ssh_OBS, SSH_lon2D, SSH_lat2D, ssh_OBS_obsper

def MLD_OBS():

	locpath='./DATA/'
	locfile='MLD_MIMOC_based_monthlyClim_rhocrit0.01.nc'
	if chkfile(locpath+locfile,zstop=True) :
		ds_obsmld = xr.open_dataset(locpath+locfile)
		lon = ds_obsmld['lon'].squeeze()
		lat = ds_obsmld['lat'].squeeze()
		mld_init = ds_obsmld['MLD'].squeeze()

	MLD_lon2D = npy.tile(lon,(lat.size,1))
	MLD_lat2D = npy.tile(lat,(lon.size,1)).T

	mld_init = xr.where( mld_init > 1e9, npy.nan, mld_init )

	# Compute the mean over the obs. monthly period 2003-2014
	mld_m03 = npy.squeeze(mld_init[2,:,:])
	mld_m09 = npy.squeeze(mld_init[8,:,:])

	return mld_m03, mld_m09, MLD_lon2D, MLD_lat2D

def MLTS_OBS():

	locpath='./DATA/'
	locfile='MIMOC_ML_v2.2_PT_S_MLP_Clim.nc'
	if chkfile(locpath+locfile,zstop=True) :
		ds_tsmld = xr.open_dataset(locpath+locfile)
		lon2D = ds_tsmld['longitude'].squeeze()
		lat2D = ds_tsmld['latitude'].squeeze()
		mlT_init = ds_tsmld['POTENTIAL_TEMPERATURE_MIXED_LAYER'].squeeze()
		mlS_init = ds_tsmld['SALINITY_MIXED_LAYER'].squeeze()
		mld_init = ds_tsmld['DEPTH_MIXED_LAYER'].squeeze()

	mlT_init = xr.where( mlT_init > 1e9, npy.nan, mlT_init )
	mlS_init = xr.where( mlS_init > 1e9, npy.nan, mlS_init )
	mld_init = xr.where( mld_init > 1e9, npy.nan, mld_init )

	return mlT_init, mlS_init, lon2D, lat2D


def EKE_OBS(t_year=1959):

	locpath='./DATA/'
	locfile='EKE_DOT_based_2003-2014.nc'
	if chkfile(locpath+locfile,zstop=True) :
		ds_ekedot = xr.open_dataset(locpath+locfile)
		lon = ds_ekedot['lon'].squeeze()
		lat = ds_ekedot['lat'].squeeze()
		EKE_init = ds_ekedot['EKE_yearly'].squeeze()

	EKE_lon2D = npy.tile(lon,(lat.size,1))
	EKE_lat2D = npy.tile(lat,(lon.size,1)).T

	if t_year >= 2003 and t_year <=2014 :
		# Get the specific year
		s_ind = (t_year-2003)
		out_EKE_OBS = npy.squeeze(EKE_init[s_ind,:,:].copy())
	else:
		# Compute the mean over the obs. monthly period 2003-2014
		out_EKE_OBS = npy.mean(EKE_init,axis=0).squeeze()

	return out_EKE_OBS,EKE_lon2D,EKE_lat2D


def ICE_THICK_OBS(zconfig='CREG025.L75',t_year=1959):

	locpath='./DATA/'
	if zconfig == 'CREG025.L75' : 
		locfile='PIOMAS_icethic_interpCREG025.L75_1-12_1979-2020.nc'
	elif zconfig == 'CREG12.L75' :
		locfile='PIOMAS_icethic_interpCREG12.L75_1-12_1979-2020.nc'
	if chkfile(locpath+locfile,zstop=True) :
		ds_icet = xr.open_dataset(locpath+locfile)
		ICE_thick_init = ds_icet['icethic'].squeeze()

	if t_year >= 1979 and t_year <= 2020 :
		s_ind = (t_year-1979)*12
		mean_ICE_thick = ICE_thick_init[s_ind:s_ind+12,:,:].copy()
	else:
		mean_ICE_thick = ICE_thick_init.copy()
	
	# Annual or climatological mean
	out_ICE_thick = npy.mean(mean_ICE_thick,axis=0).squeeze()

	return out_ICE_thick


def ICE_CONCE_OBS(t_year=1959):

	locpath='./DATA/'
	locfile='NSIDC-0051_92585_monthly.nc'
	if chkfile(locpath+locfile,zstop=True) :
		ds_icec = xr.open_dataset(locpath+locfile)
		lon = ds_icec['longitude'].squeeze()
		lat = ds_icec['latitude'].squeeze()
		CONC_init = ds_icec['Average_Sea_Ice_Concentration_with_Final_Version'].squeeze()

	# Initial data are based on add_offset
	# 251 > missing pole
	# 252 > not used
	# 253 > coastline
	# 254 > land
	# 255 > missing value
	# data will be recovereed in dividing it by 250
	CONC_init = xr.where( CONC_init == 255, npy.nan, CONC_init )
	CONC_init = xr.where( CONC_init == 254, npy.nan, CONC_init )
	CONC_init = xr.where( CONC_init == 253, npy.nan, CONC_init )
	CONC_init = xr.where( CONC_init == 251, npy.nan, CONC_init )
	COR_CONC_init = CONC_init/250.

	CONC_init_land = npy.squeeze(CONC_init[0,:,:].copy())

	if t_year >= 1979 and t_year <= 2015 : 
		CONC_m03=npy.empty((1,CONC_init.shape[1],CONC_init.shape[2]))  ;    CONC_m09=npy.empty((1,CONC_init.shape[1],CONC_init.shape[2]))
		CONC_m03[:,:,:]=npy.nan   ;    CONC_m09[:,:,:]=npy.nan
		s_indm03=(t_year-1979)*12+4	;    s_indm09=(t_year-1979)*12+10
		CONC_m03[0,:,:] = COR_CONC_init[s_indm03,:,:].copy()
		CONC_m09[0,:,:] = COR_CONC_init[s_indm09,:,:].copy()
	else:
		CONC_m03=npy.empty((37,CONC_init.shape[1],CONC_init.shape[2]))	;    CONC_m09=npy.empty((37,CONC_init.shape[1],CONC_init.shape[2]))
		CONC_m03[:,:,:]=npy.nan   ;    CONC_m09[:,:,:]=npy.nan

		# Get all March month
		c_month=4  ; it=0
		while c_month <= 446 :
			CONC_m03[it,:,:] = COR_CONC_init[c_month,:,:]
			c_month+=12  ;	it+=1

		# Get all September month
		c_month=10  ; it=0
		while c_month <= 446 :
			CONC_m09[it,:,:] = COR_CONC_init[c_month,:,:]
			c_month+=12  ;	it+=1

	# Annual march/septemper mean
	mean_CONC_m03 = npy.nanmean(CONC_m03,axis=0).squeeze()
	mean_CONC_m09 = npy.nanmean(CONC_m09,axis=0).squeeze()

	mean_CONC_m03 = xr.where( CONC_init_land == 254 , npy.nan, mean_CONC_m03 )
	mean_CONC_m09 = xr.where( CONC_init_land == 254 , npy.nan, mean_CONC_m09 )

	return mean_CONC_m03, mean_CONC_m09, lon, lat


def PHC3_OBS():

	print('				Read PHC 3.0 Obs. state  ')
	locpath='./DATA/'
	locfile='phc3.0_annual.nc'
	if chkfile(locpath+locfile) : 
		ds_obs = xr.open_dataset(locpath+locfile)
		lon_obs = ds_obs['lon']
		lat_obs = ds_obs['lat']
		PHC_lon2D=npy.tile(lon_obs,(lat_obs.size,1))
		PHC_lat2D=npy.tile(lat_obs,(lon_obs.size,1)).T

		My_varTinit = ds_obs['temp']
		My_varSinit = ds_obs['salt']

		My_varTinit = xr.where(My_varTinit > 1.e9, npy.nan, My_varTinit)
		My_varSinit = xr.where(My_varSinit > 1.e9, npy.nan, My_varSinit)

	return My_varTinit, My_varSinit, PHC_lon2D, PHC_lat2D
	
def CREG_INIT(zCONFIG,zCASE,llon,lzd):
	# Read initial state to compare with
	Ar_size=(lzd.shape[0],llon.shape[0],llon.shape[1])
	My_varSinit=npy.zeros(Ar_size)
	print('				Read initial state  ')
	locpath=zCONFIG+'/'+zCONFIG+'-'+zCASE+'-MEAN/'
	locfile=zCONFIG+'-'+zCASE+'_init_gridT.nc'
	if chkfile(locpath+locfile) : 
		ds_sinit = xr.open_dataset(locpath+locfile)
		My_varSinit = ds_sinit['vosaline'].squeeze()
	
	return My_varSinit

def Atl_plot(lon,lat,tab,contours,limits,myticks=None,name=None,zmy_cblab=None,zmy_cmap=None,filename='test.png',zvar=None, zarea=None, data_ref=False):
	#
	plt.rcParams['text.usetex']=False
	plt.rcParams['font.family']='serif'
	plt.rcParams['axes.unicode_minus'] = False
	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	if zarea == 'labsea': # Focus on Labrador Sea
		m = Basemap(width=1400000,height=1600000,lat_1=50.,lat_2=65,lon_0=-50,lat_0=59.5,projection='aea',resolution='i')
		############################################################################################################
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
	else: # Focus on North Atlantic sector
		#m = Basemap(projection='cyl',llcrnrlat=24,urcrnrlat=65,llcrnrlon=-83,urcrnrlon=-15.,resolution='i')
		#m = Basemap(projection='ortho',lat_0=45.,lon_0=-45.,resolution='i', height=1000000.)
		 m = Basemap(width=6100000,height=5000000,lat_1=30.,lat_2=70,lon_0=-45,lat_0=45,projection='aea',resolution='i')
		 #m = Basemap(width=2700000,height=2000000,lat_1=50.,lat_2=70,lon_0=-40,lat_0=60,projection='aea',resolution='i')
	
	if zvar == 'sivolu' or zvar == 'siconc'  :
		zfontsize=4.
	else:
		zfontsize=6.
	
	if zvar != 'Bathy' :
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
		C = m.contour(X,Y,tab,linewidths=zlinewidths,levels=[17.], colors=zcolor )
	else:
		C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm,extend='both')
		#zlinewidths=0.5   ; zcolor='k'
		#C2 = m.contour(X,Y,tab,linewidths=zlinewidths,levels=[-75.,-50.,-25.,25.,50.,75.], colors=zcolor )
		#plt.clabel(C2, C2.levels, inline=True, fmt='%.0f', fontsize=6)
		#zlinewidths=0.8   ; zcolor='k'
		#C3 = m.contour(X,Y,tab,linewidths=zlinewidths,levels=[0.], colors=zcolor )
		#plt.clabel(C3, C3.levels, inline=True, fmt='%.0f', fontsize=6)
	
	
		# colorbar	
		if myticks is None:
			cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
		else:
			if zvar == 'votemper' or zvar == 'vosaline' or zvar == 'sivolu' :
				cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8,drawedges=True)
			else:
				cbar = plt.colorbar(C,format='%.0f',orientation='vertical',shrink=0.8,drawedges=True)
	
			cbar.set_label(zmy_cblab,fontsize=zfontsize)
			cl = plt.getp(cbar.ax, 'ymajorticklabels')
			plt.setp(cl, fontsize=zfontsize)
	
	plt.title(name,fontsize=zfontsize)
	
	return 


def Atl_Bat(ztype='isol1000', zarea=None) :

	locpath='./'
	locfile='Bathymetry.nc'
	if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
		ds_bat = xr.open_dataset(locpath+locfile)
		lon = ds_bat['nav_lon'].squeeze()
		lat = ds_bat['nav_lat'].squeeze()
		My_var = ds_bat['bathy_meter'].squeeze()
	
	spval = 0.
	My_var = xr.where( My_var <= spval, npy.nan, My_var )
	#My_var= npy.ma.masked_where(My_var <= spval,My_var)
	
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
	zcolorbat='grey'  ;  zalpha=0.7
	if zarea == 'labsea': # Focus on Labrador Sea
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
	else: # Focus on North Atlantic sector
		#m = Basemap(projection='cyl',llcrnrlat=20,urcrnrlat=60,llcrnrlon=-90,urcrnrlon=0,resolution='i')
		m = Basemap(width=6100000,height=5000000,lat_1=30.,lat_2=70,lon_0=-45,lat_0=45,projection='aea',resolution='i')
		zcolorbat='grey'   ;  zalpha=0.7

	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	pal = plt.get_cmap('binary')
	X,Y = m(lon.values,lat.values)

	# contour (optional)
	CS2 = m.contour(X, Y, My_var.values, linewidths=0.5,levels=contours, colors=zcolorbat, alpha=zalpha)
	plt.clabel(CS2, CS2.levels, inline=True, fmt='%.0f', fontsize=4)

	return m, X, Y
