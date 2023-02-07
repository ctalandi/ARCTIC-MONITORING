#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import sys
import numpy as npy
from CREG_maps_cont import *
from checkfile import *
import subprocess

# Matplotlib
try:
        import matplotlib.pylab as plt
        import matplotlib as mpl
        from matplotlib import rcParams
except:
        print 'matplotlib is not available on your machine'
        print 'check python path or install this package' ; exit()

# Basemap
try:
        from mpl_toolkits.basemap import Basemap
except:
        print 'Basemap is not available on your machine'
        print 'check python path or install this package' ; exit()

# NetCDF 
try:
	from netCDF4 import Dataset
except:
        print 'netCDF4 is not available on your machine'
        print 'check python path or install this package' ; exit()


def MTS_maps( zlon, zlat, zCONF, zCASE, zMy_var1SeasM, zMy_var1SeasS, zMy_varT, zMy_varS, zgdept_0, ze3t_0, zclimyear, indir, grdir ) :

	e3t_0msk_SeasM=npy.zeros((ze3t_0.shape[0],zlon.shape[0],zlon.shape[1]))
	e3t_0msk_SeasS=npy.zeros((ze3t_0.shape[0],zlon.shape[0],zlon.shape[1]))
	e3t_0sum_SeasM=npy.zeros((zlon.shape[0],zlon.shape[1]))
	e3t_0sum_SeasS=npy.zeros((zlon.shape[0],zlon.shape[1]))
	T_mldM=npy.zeros((2,zlon.shape[0],zlon.shape[1]))
	S_mldM=npy.zeros((2,zlon.shape[0],zlon.shape[1]))

	# Mask all levels below the MLD 
	e3t_0msk_SeasM= npy.squeeze(npy.where(zgdept_0[:,:,:] < zMy_var1SeasM[:,:], ze3t_0[:,:,:], 0.))
	e3t_0msk_SeasS= npy.squeeze(npy.where(zgdept_0[:,:,:] < zMy_var1SeasS[:,:], ze3t_0[:,:,:], 0.))

	# Sum all e3t scale factors over the vertcal axis
	e3t_0sum_SeasM= e3t_0msk_SeasM.sum(axis=0).squeeze()
	e3t_0sum_SeasS= e3t_0msk_SeasS.sum(axis=0).squeeze()

	# Compute the T/S mean within the ML in March/September
	T_mldM[0,:,:] = (npy.sum(e3t_0msk_SeasM[:,:,:]*zMy_varT[0,:,:,:].squeeze(),axis=0))/e3t_0sum_SeasM[:,:]
	S_mldM[0,:,:] = (npy.sum(e3t_0msk_SeasM[:,:,:]*zMy_varS[0,:,:,:].squeeze(),axis=0))/e3t_0sum_SeasM[:,:]

	T_mldM[1,:,:] = (npy.sum(e3t_0msk_SeasS[:,:,:]*zMy_varT[1,:,:,:].squeeze(),axis=0))/e3t_0sum_SeasS[:,:]
	S_mldM[1,:,:] = (npy.sum(e3t_0msk_SeasS[:,:,:]*zMy_varS[1,:,:,:].squeeze(),axis=0))/e3t_0sum_SeasS[:,:]

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
	my_cmap= plt.cm.get_cmap('Spectral_r')
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, T_mldM[0,:,:], contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap)

	plt.subplot(222)
        ztitle=zCASE +' September ML mean T over \n'+str(zclimyear)
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, T_mldM[1,:,:].squeeze(), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap)

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
	my_cmap= plt.cm.get_cmap('Spectral_r')
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, S_mldM[0,:,:].squeeze(), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar='MLTSS')

	plt.subplot(222)
        ztitle=zCASE +' September ML mean S over \n'+str(zclimyear)
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, S_mldM[1,:,:].squeeze(), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar='MLTSS')

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

	savefile=True
        if savefile:
		# ML T/S mean field 
		#######################
		cmd_ddate="date"  ;  get_output = subprocess.check_output(cmd_ddate)
                nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_MTSClim_'+'y'+zclimyear+'.nc'
                w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
                w_nc_fid.description = "Diagnostics has been calculated using the Arctic monitoring tool "
		w_nc_fid.date=get_output.decode("utf-8")
                w_nc_fid.createDimension('x', S_mldM.shape[2])
                w_nc_fid.createDimension('y', S_mldM.shape[1])
                w_nc_fid.createDimension('xobs', lat_obs.shape[1])
                w_nc_fid.createDimension('yobs', lat_obs.shape[0])
                #w_nc_fid.createDimension('time_counter', None)

                w_nc_var = w_nc_fid.createVariable('S_mldM_march', 'f4', ('y','x'))
                w_nc_var.long_name='Model ML mean S in March'
                w_nc_var.units="PSU"
                w_nc_fid.variables['S_mldM_march'][:,:] = S_mldM[0,:,:].squeeze()

                w_nc_var = w_nc_fid.createVariable('S_mldM_septe', 'f4', ('y','x'))
                w_nc_var.long_name='Model ML mean S in September'
                w_nc_var.units="PSU"
                w_nc_fid.variables['S_mldM_septe'][:,:] = S_mldM[1,:,:].squeeze()

                w_nc_var = w_nc_fid.createVariable('T_mldM_march', 'f4', ('y','x'))
                w_nc_var.long_name='Model ML mean T in March'
                w_nc_var.units="DegC"
                w_nc_fid.variables['T_mldM_march'][:,:] = T_mldM[0,:,:].squeeze()

                w_nc_var = w_nc_fid.createVariable('T_mldM_septe', 'f4', ('y','x'))
                w_nc_var.long_name='Model ML mean T in September'
                w_nc_var.units="DegC"
                w_nc_fid.variables['T_mldM_septe'][:,:] = T_mldM[1,:,:].squeeze()

                w_nc_var = w_nc_fid.createVariable('mlS_obs_march', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='MIMOC ML mean S in March'
                w_nc_var.units="PSU"
                w_nc_fid.variables['mlS_obs_march'][:,:] = mlS_obs[2,:,:].squeeze()

                w_nc_var = w_nc_fid.createVariable('mlS_obs_septe', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='MIMOC ML mean S in September'
                w_nc_var.units="PSU"
                w_nc_fid.variables['mlS_obs_septe'][:,:] = mlS_obs[8,:,:].squeeze()

                w_nc_var = w_nc_fid.createVariable('mlT_obs_march', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='MIMOC ML mean T in March'
                w_nc_var.units="DegC"
                w_nc_fid.variables['mlT_obs_march'][:,:] = mlT_obs[2,:,:].squeeze()

                w_nc_var = w_nc_fid.createVariable('mlT_obs_septe', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='MIMOC ML mean T in September'
                w_nc_var.units="DegC"
                w_nc_fid.variables['mlT_obs_septe'][:,:] = mlT_obs[8,:,:].squeeze()

                w_nc_var = w_nc_fid.createVariable('lat_obs', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='Degrees north'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lat_obs'][:,:] = lat_obs

                w_nc_var = w_nc_fid.createVariable('lon_obs', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='Degrees east'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lon_obs'][:,:] = lon_obs
                
                w_nc_var = w_nc_fid.createVariable('lat_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Degrees north'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lat_mod'][:,:] = zlat

                w_nc_var = w_nc_fid.createVariable('lon_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Degrees east'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lon_mod'][:,:] = zlon
                
                w_nc_fid.close()  # close the file

	return

def simple_maps( zlon, zlat, zCONF, zCASE, zMy_var1, zMyvar, zclimyear, slev=None, seas='', zfram=111, plot_obs=0, ano=0) :

        m_alpha=1.
	
	# Do the plot 
	print zMyvar+' plot'
	print 
	
	plt.subplot(zfram)
	contours,limits,myticks,ztitle,zfile_ext,my_cblab,my_cmap,m_alpha = SET_ARC_CNT(zCASE,zclimyear,seas,zMyvar,zslev=slev,zplot_obs=plot_obs,zdiff=ano)
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, zMy_var1[:,:]*m_alpha, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar)

	return

def AWTmax_maps( zlon, zlat, zMy_var1, zMy_var1S, zdepth, zMyvar, zCONFIG, zCASE, zclimyear) :

	zCONFCASE=zCONFIG+'-'+zCASE

	zMy_varTinit, zMy_varSinit ,lon_obs, lat_obs=PHC3_OBS()

	# Start to mask field where Salinity is lower than 33.5 PSU
	# Necessary to remove the surface temperature maxima than can arise
	zMy_var1=npy.ma.masked_where((zMy_var1S <= 33.5),zMy_var1)
	zMy_varTinit=npy.ma.masked_where((zMy_varSinit <= 33.5),zMy_varTinit)

	# Find the Max temp. over depth
        zAWTmax1=npy.amax(zMy_var1[:,:,:],axis=0).squeeze()
        zAWTmaxI=npy.amax(zMy_varTinit[:,:,:],axis=0).squeeze()

	# Find the depth of the Max temp. 
        zAWTmax_depth1_ind=npy.argmax(zMy_var1[:,:,:],axis=0).squeeze()
        zAWTmax_depthI_ind=npy.argmax(zMy_varTinit[:,:,:],axis=0).squeeze()

	zAWTmax_depth1=npy.zeros((zAWTmax_depth1_ind.shape[0],zAWTmax_depth1_ind.shape[1]))
	zAWTmax_depthI=npy.zeros((zAWTmax_depthI_ind.shape[0],zAWTmax_depthI_ind.shape[1]))

	zdepth_phc3=npy.array([0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, \
                     700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2500, \
                     3000, 3500, 4000, 4500, 5000, 5500 ])

	for jj in set(npy.arange(zAWTmax_depth1_ind.shape[0])) :
	    for ji in set(npy.arange(zAWTmax_depth1_ind.shape[1])) :
        	zAWTmax_depth1[jj,ji]=zdepth[zAWTmax_depth1_ind[jj,ji]]

	for jj in set(npy.arange(zAWTmax_depthI.shape[0])) :
	    for ji in set(npy.arange(zAWTmax_depthI.shape[1])) :
        	zAWTmax_depthI[jj,ji]=zdepth_phc3[zAWTmax_depthI_ind[jj,ji]]

	zAWTmax_depth1=npy.ma.masked_where((zMy_var1S[0,:,:].squeeze() == 0),zAWTmax_depth1)
	zAWTmax_depthI=npy.ma.masked_where((zMy_varSinit[0,:,:].squeeze() == 0),zAWTmax_depthI)


	# Make the plot for the AW Max Temp 
	#############################################################################################
        vmin=0. ; vmax=7. ; vint=0.5
        contours=npy.arange(vmin,vmax+vint,vint)  
        limits=[vmin,vmax,vint]                  
        myticks=npy.arange(vmin,vmax+vint,vint) 
	
        ztitle=zCASE +' AW Max Temp over \n'+zclimyear
        my_cblab=r'($^\circ$C)'
	my_cmap= plt.cm.get_cmap('jet')

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
        my_cmap= plt.cm.get_cmap('jet')

	plt.subplot(223) 
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(zlon, zlat, zAWTmax_depth1, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap)

	plt.subplot(224) 
        ztitle=' AW Max Temp depth from \n'+' PHC 3.0'
	zoutmap=Arc_Bat(ztype='isol1000')
	Arc_plot(lon_obs, lat_obs, zAWTmax_depthI, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap)
	
	plt.tight_layout()
	
	plt.savefig(zCONFCASE+zfile_ext+'y'+zclimyear+'.pdf')

	savefile=True
        if savefile:
		# FWC field 
		#######################
		cmd_ddate="date"  ;  get_output = subprocess.check_output(cmd_ddate)
                nc_f = './NETCDF/'+zCONFIG+'-'+zCASE+'_AWTmaxClim_'+'y'+zclimyear+'.nc'
                w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
                w_nc_fid.description = "Diagnostics has been calculated using the Arctic monitoring tool "
		w_nc_fid.date=get_output.decode("utf-8")
                w_nc_fid.createDimension('x', zAWTmax1.shape[1])
                w_nc_fid.createDimension('y', zAWTmax1.shape[0])
                w_nc_fid.createDimension('xobs', lat_obs.shape[1])
                w_nc_fid.createDimension('yobs', lat_obs.shape[0])
                #w_nc_fid.createDimension('time_counter', None)

                w_nc_var = w_nc_fid.createVariable('AWTmax_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Model AWT max calculated using a salinity criteria value Sref=33.5'
                w_nc_var.units="DegC"
                w_nc_fid.variables['AWTmax_mod'][:,:] = zAWTmax1

                w_nc_var = w_nc_fid.createVariable('AWTmax_init', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='Model Initial AWT max calculated using a salinity criteria value Sref=33.5'
                w_nc_var.units="DegC"
                w_nc_fid.variables['AWTmax_init'][:,:] = zAWTmaxI

                w_nc_var = w_nc_fid.createVariable('AWTmaxDepth_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Model AWT max depth calculated using a salinity criteria value Sref=33.5'
                w_nc_var.units="m"
                w_nc_fid.variables['AWTmaxDepth_mod'][:,:] = zAWTmax_depth1

                w_nc_var = w_nc_fid.createVariable('AWTmaxDepth_init', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='Model AWT max depth calculated using a salinity criteria value Sref=33.5'
                w_nc_var.units="m"
                w_nc_fid.variables['AWTmaxDepth_init'][:,:] = zAWTmax_depthI

                w_nc_var = w_nc_fid.createVariable('lat_obs', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='Degrees north'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lat_obs'][:,:] = lat_obs

                w_nc_var = w_nc_fid.createVariable('lon_obs', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='Degrees east'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lon_obs'][:,:] = lon_obs
                
                w_nc_var = w_nc_fid.createVariable('lat_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Degrees north'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lat_mod'][:,:] = zlat

                w_nc_var = w_nc_fid.createVariable('lon_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Degrees east'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lon_mod'][:,:] = zlon
                
                w_nc_fid.close()  # close the file

	return

def FWC_maps( zlon, zlat, zMy_var1S, zMy_varSinit, zMy_var1ssh, zCONFIG, zCASE, zclimyear, ze3, ztmask, ztime_dim) :

	# FWC calculation over the year
	###########################################
	Sref=34.80 
	 
	print '			FWC calculation & plot '

	# Use the initial state as reference
	ze33Dtime=ze3[:,:,:]*ztmask[:,:,:]
	ze33Dtime_msk = ze33Dtime.copy()
	ze33Dtime_mskI = ze33Dtime.copy()
	ze33Dtime_msk[npy.where(zMy_var1S > Sref )] = 0.e0
	ze33Dtime_mskI[npy.where(zMy_varSinit > Sref )] = 0.e0
	fwc2D = npy.zeros((zlon.shape[0],zlon.shape[1]))    ;    fwc2D_init = npy.zeros((zlon.shape[0],zlon.shape[1]))
	
	zMy_varSinit=CREG_INIT(zCONFIG,zCASE,zlon,ze3)

	#print a specific point to debug within the Beaufoirt Gyre
	dbg=False
	if zCONFIG == 'CREG025.L75' :
		idbg=185 ; jdbg=515 ; kdbg=50   # CREG025.L75
	elif zCONFIG == 'CREG12.L75' :
		idbg=562 ; jdbg=1558 ; kdbg=50   # CREG12.L75
	
	for jk in set(npy.arange(75)) :
	    fwc2D[:,:] = fwc2D[:,:] + ( Sref - zMy_var1S[jk,:,:] ) / Sref * ze33Dtime_msk[jk,:,:] 
	    fwc2D_init[:,:] = fwc2D_init[:,:] + ( Sref - zMy_varSinit[jk,:,:] ) / Sref * ze33Dtime_mskI[jk,:,:] 
	
	if dbg:
		print 'ze33Dtime_msk[ti,0:kdbg,jdbg,idbg]',ze33Dtime_msk[0,0:kdbg,jdbg,idbg], zMy_var1S[0,0:kdbg,jdbg,idbg]
		print 'fwc2D[jdbg,idbg]',fwc2D[jdbg,idbg]
	
	tmask2D=npy.squeeze(ztmask[0,:,:].copy())
	fwc2D=npy.ma.masked_where((tmask2D == 0),fwc2D)
	fwc2D_init=npy.ma.masked_where((tmask2D == 0),fwc2D_init)

	# Use obs. from Proshutinsky GRL2018
	locpath='./DATA/'
	locfile='BeaufortGyreFWC-Obs-Proshutinsky_GRL2018_y2003-2017.nc'
	if chkfile(locpath+locfile,zstop=True) :
		fid=Dataset(locpath+locfile)
		fwc2D_obs=npy.squeeze(fid.variables['BGFWC'])
		lat2D_obs=npy.squeeze(fid.variables['lat'])
		lon2D_obs=npy.squeeze(fid.variables['lon'])

		lat2D_obs=npy.ma.masked_where(lat2D_obs == 9999.,lat2D_obs)
		lon2D_obs=npy.ma.masked_where(lon2D_obs == 9999.,lon2D_obs)

	if npy.int(zclimyear) >= 2003 and npy.int(zclimyear) <= 2017 :
		s_ind=npy.int(zclimyear)-2003  ; e_ind=s_ind+1
		mean_FWCObs=npy.squeeze(npy.array(fwc2D_obs[s_ind:e_ind,:,:]).copy())
		mean_FWCObs=npy.ma.masked_where(mean_FWCObs == 9999.,mean_FWCObs)
		obsper=zclimyear
	else:
		mean_FWCObs=npy.squeeze(npy.nanmean(npy.array(fwc2D_obs),axis=0))
		mean_FWCObs=npy.ma.masked_where(mean_FWCObs == 9999.,mean_FWCObs)
		obsper='2003-2017'
	
	# Get observations SSH
	obs_ssh,lon_obs,lat_obs,obs_ssh_per=SSH_OBS(t_year=npy.int(zclimyear[0:4]))
	obs_ssh=npy.ma.masked_where(obs_ssh >= 9e20, obs_ssh)

	plt.clf()
	fig=plt.figure()

	# Plot the FWC map mean over the year
	#####################################
	my_cblab=r'(m)'
	zfile_ext='_FWCClim_'
	ztitle=zCASE+' mean FWC (m) over \n'+zclimyear
	vmin=0. ; vmax=25. ; vint=2.
	contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	limits=[vmin,vmax,vint]                        # limits for eke
	myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)
	my_cmap=plt.cm.get_cmap('coolwarm')
	my_cmap=plt.cm.get_cmap('Spectral_r')
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
	plt.savefig(zCONFIG+'-'+zCASE+zfile_ext+'y'+zclimyear+'.pdf')


	savefile=True
        if savefile:
		# FWC field 
		#######################
		cmd_ddate="date"  ;  get_output = subprocess.check_output(cmd_ddate)
                nc_f = './NETCDF/'+zCONFIG+'-'+zCASE+'_FWCClim_'+'y'+zclimyear+'.nc'
                w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
                w_nc_fid.description = "Diagnostics have been calculated using the Arctic monitoring tool "
		w_nc_fid.date=get_output.decode("utf-8")
                w_nc_fid.createDimension('x', fwc2D.shape[1])
                w_nc_fid.createDimension('y', fwc2D.shape[0])
                w_nc_fid.createDimension('xobs', lat2D_obs.shape[1])
                w_nc_fid.createDimension('yobs', lat2D_obs.shape[0])
                #w_nc_fid.createDimension('time_counter', None)

                w_nc_var = w_nc_fid.createVariable('fwc_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Model FWC calculated using a salinity reference value Sref=31.8'
                w_nc_var.units="m"
		fwc2D=npy.ma.masked_where(ztmask[0,:,:] == 0, fwc2D)
                w_nc_fid.variables['fwc_mod'][:,:] = fwc2D

                w_nc_var = w_nc_fid.createVariable('fwc_init', 'f4', ('y','x'))
                w_nc_var.long_name='Model Initial state FWC calculated using a salinity reference value Sref=31.8'
                w_nc_var.units="m"
		fwc2D_init=npy.ma.masked_where(ztmask[0,:,:] == 0, fwc2D_init)
                w_nc_fid.variables['fwc_init'][:,:] = fwc2D_init

                w_nc_var = w_nc_fid.createVariable('fwc_obs', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='Obs. FWC calculated using a salinity reference value Sref=31.8 '+ \
				   ' from Proshutinsky et al. (GRL2018). Considered period: '+obsper
                w_nc_var.units="m"
                w_nc_fid.variables['fwc_obs'][:,:] = mean_FWCObs

                w_nc_var = w_nc_fid.createVariable('lat_obs', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='Degrees north'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lat_obs'][:,:] = lat2D_obs

                w_nc_var = w_nc_fid.createVariable('lon_obs', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='Degrees east'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lon_obs'][:,:] = lon2D_obs
                
                w_nc_var = w_nc_fid.createVariable('lat_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Degrees north'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lat_mod'][:,:] = zlat

                w_nc_var = w_nc_fid.createVariable('lon_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Degrees east'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lon_mod'][:,:] = zlon
                
                w_nc_fid.close()  # close the file

		# SSH field 
		#######################
		cmd_ddate="date"  ;  get_output = subprocess.check_output(cmd_ddate)
                nc_f = './NETCDF/'+zCONFIG+'-'+zCASE+'_SSHClim_'+'y'+zclimyear+'.nc'
                w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
                w_nc_fid.description = "Diagnostics have been calculated using the Arctic monitoring tool "
                w_nc_fid.date = get_output.decode("utf-8")
                w_nc_fid.createDimension('x', zMy_var1ssh.shape[1])
                w_nc_fid.createDimension('y', zMy_var1ssh.shape[0])
                w_nc_fid.createDimension('xobs', lon_obs.shape[1])
                w_nc_fid.createDimension('yobs', lon_obs.shape[0])
                #w_nc_fid.createDimension('time_counter', None)

                w_nc_var = w_nc_fid.createVariable('ssh_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Model SSH '
                w_nc_var.units="m"
		zMy_var1ssh=npy.ma.masked_where(ztmask[0,:,:] == 0, zMy_var1ssh)
                w_nc_fid.variables['ssh_mod'][:,:] = zMy_var1ssh

                w_nc_var = w_nc_fid.createVariable('lat_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Degrees north'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lat_mod'][:,:] = zlat

                w_nc_var = w_nc_fid.createVariable('lon_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Degrees east'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lon_mod'][:,:] = zlon
                
                w_nc_var = w_nc_fid.createVariable('ssh_obs', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='Observed DOT from Armitage et al. 2017 for the period: '+obs_ssh_per    
                w_nc_var.units="m"
                w_nc_fid.variables['ssh_obs'][:,:] = obs_ssh

                w_nc_var = w_nc_fid.createVariable('lat_obs', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='Degrees north'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lat_obs'][:,:] = lat_obs

                w_nc_var = w_nc_fid.createVariable('lon_obs', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='Degrees east'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lon_obs'][:,:] = lon_obs

                w_nc_fid.close()  # close the file

        return


def Arc_plot(lon,lat,tab,contours,limits,myticks=None,name=None,zmy_cblab=None,zmy_cmap=None,filename='test.pdf',zvar=None):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
#	rcParams['text.usetex']=True
#	rcParams['text.latex.unicode']=True
	rcParams['font.family']='serif'

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
		pal = plt.cm.get_cmap('coolwarm')
		#pal = plt.cm.get_cmap('terrain')

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
		if zvar == 'votemper' or zvar == 'vosaline' or zvar == 'sivolu' :
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
		fieldbat=Dataset(locpath+locfile)
		lon  = npy.squeeze(fieldbat.variables['nav_lon'])
		lat  = npy.squeeze(fieldbat.variables['nav_lat'])
		My_var = npy.squeeze(fieldbat.variables['Bathymetry'])
	
	spval = 0.
	My_var= npy.ma.masked_where(My_var <= spval,My_var)
	
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
	rcParams['text.latex.unicode']=True
	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
        m = Basemap(projection='npstere',boundinglat=65,lon_0=-60, resolution='i')
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	pal = plt.cm.get_cmap('binary')
	X,Y = m(lon,lat)

	# contour (optional)
        CS2 = m.contour(X, Y, My_var, linewidths=0.5,levels=contours, colors='grey', alpha=0.4)
        plt.clabel(CS2, CS2.levels, inline=True, fmt='%.0f', fontsize=3)

	return m, X, Y

def SSH_OBS(t_year=1959):

	locpath='./DATA/'
	locfile='EKE_DOT_based_2003-2014.nc'
	if chkfile(locpath+locfile,zstop=True) :
		field=Dataset(locpath+locfile)
		lon = npy.array(field.variables['lon'])
		lat = npy.array(field.variables['lat'])
		ssh_init = npy.squeeze(field.variables['DOT'])

        SSH_lon2D=npy.tile(lon,(lat.size,1))
        SSH_lat2D=npy.tile(lat,(lon.size,1)).T

	if t_year >= 2003 and t_year <=2014 :
		# Get the specific year
		s_ind=(t_year-2003)*12   
		out_ssh_OBS = npy.mean(ssh_init[s_ind:s_ind+12,:,:],axis=0).squeeze()
		ssh_OBS_obsper=str(t_year)
	else:
		# Compute the mean over the obs. monthly period 2003-2014
		out_ssh_OBS = npy.mean(ssh_init,axis=0).squeeze()
		ssh_OBS_obsper='2003-2014'

	# Remove the domain mean to get an anomaly
	out_ssh_OBS = out_ssh_OBS - npy.nanmean(ssh_init)

	return out_ssh_OBS, SSH_lon2D, SSH_lat2D, ssh_OBS_obsper

def MLD_OBS():

	locpath='./DATA/'
	locfile='MLD_MIMOC_based_monthlyClim_rhocrit0.01.nc'
	if chkfile(locpath+locfile,zstop=True) :
		field=Dataset(locpath+locfile)
		lon = npy.array(field.variables['lon'])
		lat = npy.array(field.variables['lat'])
		mld_init = npy.squeeze(field.variables['MLD'])

        MLD_lon2D=npy.tile(lon,(lat.size,1))
        MLD_lat2D=npy.tile(lat,(lon.size,1)).T

	mld_init=npy.ma.masked_where(mld_init > 1e9,mld_init)

	# Compute the mean over the obs. monthly period 2003-2014
	mld_m03 = npy.squeeze(mld_init[2,:,:])
	mld_m09 = npy.squeeze(mld_init[8,:,:])

	return mld_m03, mld_m09, MLD_lon2D, MLD_lat2D

def MLTS_OBS():

	locpath='./DATA/'
	locfile='MIMOC_ML_v2.2_PT_S_MLP_Clim.nc'
	if chkfile(locpath+locfile,zstop=True) :
		field=Dataset(locpath+locfile)
		lon2D = npy.array(field.variables['longitude'])
		lat2D = npy.array(field.variables['latitude'])
		mlT_init = npy.squeeze(field.variables['POTENTIAL_TEMPERATURE_MIXED_LAYER'])
		mlS_init = npy.squeeze(field.variables['SALINITY_MIXED_LAYER'])
		mld_init = npy.squeeze(field.variables['DEPTH_MIXED_LAYER'])

	mlT_init=npy.ma.masked_where(mlT_init > 1e9,mlT_init)
	mlS_init=npy.ma.masked_where(mlS_init > 1e9,mlS_init)
	mld_init=npy.ma.masked_where(mld_init > 1e9,mld_init)

	return mlT_init, mlS_init, lon2D, lat2D


def EKE_OBS(t_year=1959):

	locpath='./DATA/'
	locfile='EKE_DOT_based_2003-2014.nc'
	if chkfile(locpath+locfile,zstop=True) :
		field=Dataset(locpath+locfile)
		lon = npy.array(field.variables['lon'])
		lat = npy.array(field.variables['lat'])
		EKE_init = npy.squeeze(field.variables['EKE_yearly'])

        EKE_lon2D=npy.tile(lon,(lat.size,1))
        EKE_lat2D=npy.tile(lat,(lon.size,1)).T

	if t_year >= 2003 and t_year <=2014 :
		# Get the specific year
		s_ind=(t_year-2003)
		out_EKE_OBS = npy.squeeze(EKE_init[s_ind,:,:].copy())
	else:
		# Compute the mean over the obs. monthly period 2003-2014
		out_EKE_OBS = npy.mean(EKE_init,axis=0).squeeze()

	return out_EKE_OBS,EKE_lon2D,EKE_lat2D


def ICE_THICK_OBS(zconfig='CREG025.L75',t_year=1959):

	locpath='./DATA/'
	if zconfig == 'CREG025.L75' : 
		locfile='PIOMAS_icethic_interpCREG025.L75_1-12_1979-2018.nc'
	elif zconfig == 'CREG12.L75' :
		locfile='PIOMAS_icethic_interpCREG12.L75_1-12_1979-2018.nc'
	if chkfile(locpath+locfile,zstop=True) :
		field=Dataset(locpath+locfile)
		ICE_thick_init = npy.squeeze(field.variables['icethic'])

	if t_year >= 1979 and t_year <= 2018 :
		s_ind=(t_year-1979)*12
		mean_ICE_thick=ICE_thick_init[s_ind:s_ind+12,:,:].copy()
	else:
		mean_ICE_thick=ICE_thick_init.copy()
	
	# Annual or climatological mean
	out_ICE_thick=npy.mean(mean_ICE_thick,axis=0).squeeze()

	return out_ICE_thick


def ICE_CONCE_OBS(t_year=1959):

	locpath='./DATA/'
	locfile='NSIDC-0051_92585_monthly.nc'
	if chkfile(locpath+locfile,zstop=True) :
		field=Dataset(locpath+locfile)
		lon = npy.array(field.variables['longitude'])
		lat = npy.array(field.variables['latitude'])
		CONC_init = npy.squeeze(field.variables['Average_Sea_Ice_Concentration_with_Final_Version'])

	# Initial data are based on add_offset
	# 251 > missing pole
	# 252 > not used
	# 253 > coastline
	# 254 > land
	# 255 > missing value
	# data will be recovereed in dividing it by 250
	CONC_init=npy.ma.masked_where(CONC_init == 255,CONC_init)
	CONC_init=npy.ma.masked_where(CONC_init == 254,CONC_init)
	CONC_init=npy.ma.masked_where(CONC_init == 253,CONC_init)
	CONC_init=npy.ma.masked_where(CONC_init == 251,CONC_init)
        COR_CONC_init = CONC_init/250.

	CONC_init_land=npy.squeeze(CONC_init[0,:,:].copy())

	if t_year >= 1979 and t_year <= 2015 : 
		CONC_m03=npy.empty((1,CONC_init.shape[1],CONC_init.shape[2]))  ;    CONC_m09=npy.empty((1,CONC_init.shape[1],CONC_init.shape[2]))
		CONC_m03[:,:,:]=npy.nan   ;    CONC_m09[:,:,:]=npy.nan
		s_indm03=(t_year-1979)*12+4     ;    s_indm09=(t_year-1979)*12+10
		CONC_m03[0,:,:] = COR_CONC_init[s_indm03,:,:].copy()
		CONC_m09[0,:,:] = COR_CONC_init[s_indm09,:,:].copy()
	else:
		CONC_m03=npy.empty((37,CONC_init.shape[1],CONC_init.shape[2]))  ;    CONC_m09=npy.empty((37,CONC_init.shape[1],CONC_init.shape[2]))
		CONC_m03[:,:,:]=npy.nan   ;    CONC_m09[:,:,:]=npy.nan

		# Get all March month
		c_month=4  ; it=0
		while c_month <= 446 :
			CONC_m03[it,:,:] = COR_CONC_init[c_month,:,:]
			c_month+=12  ;  it+=1

		# Get all September month
		c_month=10  ; it=0
		while c_month <= 446 :
			CONC_m09[it,:,:] = COR_CONC_init[c_month,:,:]
			c_month+=12  ;  it+=1

	# Annual march/septemper mean
	mean_CONC_m03 =npy.nanmean(CONC_m03,axis=0).squeeze()
	mean_CONC_m09 =npy.nanmean(CONC_m09,axis=0).squeeze()

	mean_CONC_m03=npy.ma.masked_where(npy.squeeze(CONC_init_land) == 254 ,mean_CONC_m03)
	mean_CONC_m09=npy.ma.masked_where(npy.squeeze(CONC_init_land) == 254 ,mean_CONC_m09)

	return mean_CONC_m03, mean_CONC_m09, lon, lat


def PHC3_OBS():
	Ar_size=(33,180,360)
	My_varTinit=npy.zeros(Ar_size)   ; My_varSinit=npy.zeros(Ar_size)
	print ' 			Read initial state  '
	locpath='./DATA/'
	locfile='phc3.0_annual.nc'
	if chkfile(locpath+locfile) : 
		field=Dataset(locpath+locfile)
		lon_obs = npy.array(field.variables['lon'])
		lat_obs = npy.array(field.variables['lat'])
                PHC_lon2D=npy.tile(lon_obs,(lat_obs.size,1))
                PHC_lat2D=npy.tile(lat_obs,(lon_obs.size,1)).T

		My_varTinit[:,:,:] = npy.array(field.variables['temp'])
		My_varSinit[:,:,:] = npy.array(field.variables['salt'])
		My_varTinit=npy.ma.masked_where(My_varTinit > 1.e9, My_varTinit)
		My_varSinit=npy.ma.masked_where(My_varSinit > 1.e9, My_varSinit)

		My_varTinit=npy.ma.masked_where(npy.isnan(My_varTinit),My_varTinit)
		My_varSinit=npy.ma.masked_where(npy.isnan(My_varSinit),My_varSinit)
	
	
	return My_varTinit, My_varSinit, PHC_lon2D, PHC_lat2D
	
def CREG_INIT(zCONFIG,zCASE,llon,lzd):
	# Read initial state to compare with
	Ar_size=(lzd.shape[0],llon.shape[0],llon.shape[1])
	My_varSinit=npy.zeros(Ar_size)
	print ' 			Read initial state  '
	locpath=zCONFIG+'/'+zCONFIG+'-'+zCASE+'-MEAN/'
	locfile=zCONFIG+'-'+zCASE+'_init_gridT.nc'
	if chkfile(locpath+locfile) : 
		field=Dataset(locpath+locfile)
		My_varSinit[:,:,:] = npy.squeeze(field.variables['vosaline'])
	
	return My_varSinit

def Atl_plot(lon,lat,tab,contours,limits,myticks=None,name=None,zmy_cblab=None,zmy_cmap=None,filename='test.pdf',zvar=None, zarea=None, data_ref=False):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
#	rcParams['text.usetex']=True
	rcParams['font.family']='serif'

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
		pal = plt.cm.get_cmap('coolwarm')
		#pal = plt.cm.get_cmap('terrain')

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
		fieldbat=Dataset(locpath+locfile)
		lon  = npy.squeeze(fieldbat.variables['nav_lon'])
		lat  = npy.squeeze(fieldbat.variables['nav_lat'])
		My_var = npy.squeeze(fieldbat.variables['Bathymetry'])
	
	spval = 0.
	My_var= npy.ma.masked_where(My_var <= spval,My_var)
	
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
	rcParams['text.latex.unicode']=True
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
	pal = plt.cm.get_cmap('binary')
	X,Y = m(lon,lat)

	# contour (optional)
        CS2 = m.contour(X, Y, My_var, linewidths=0.5,levels=contours, colors=zcolorbat, alpha=zalpha)
        plt.clabel(CS2, CS2.levels, inline=True, fmt='%.0f', fontsize=4)

	return m, X, Y
