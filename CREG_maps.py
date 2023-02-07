#!/usr/bin/env python

import sys 
import matplotlib
matplotlib.use('Agg')
import PyRaf
import numpy as npy
import CREG_maps_func
from checkfile import *
import matplotlib.pylab as plt
from netCDF4 import Dataset
import matplotlib as mpl
import subprocess

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

print
print '                              Configuration :' , CONFCASE
print '                              Period        :' , str(s_year),' - ',str(e_year)
print


########################################
# Read GRID 
########################################
#------------------------------------------------------------------------------------------------------------------------
# Read only once for the same section type
locpath=grid_dir
locfile=CONFCASE+'_mask.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	fieldmask=Dataset(locpath+locfile)
	lon = npy.array(fieldmask.variables['nav_lon'])
	lat = npy.array(fieldmask.variables['nav_lat'])
	tmask = npy.squeeze(fieldmask.variables['tmask'])
	fmask = npy.squeeze(fieldmask.variables['fmask'])

locpath=grid_dir
locfile=CONFCASE+'_mesh_hgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	fieldhmesh=Dataset(locpath+locfile)
	e1t= npy.squeeze(fieldhmesh.variables['e1t'])
	e2t= npy.squeeze(fieldhmesh.variables['e2t'])

locpath=grid_dir
locfile=CONFCASE+'_mesh_zgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	fieldzmesh=Dataset(locpath+locfile)
	ze3 = npy.squeeze(fieldzmesh.variables['e3t_0'])
	z   = npy.squeeze(fieldzmesh.variables['gdept_1d'])
	gdept_0 = npy.squeeze(fieldzmesh.variables['gdept_0'])

tmask2D=npy.squeeze(tmask[0,:,:].copy())
fmask2D=npy.squeeze(fmask[0,:,:].copy())

e1te2t = (e1t*e2t).squeeze()
e1te2t=npy.ma.masked_where((tmask2D == 0.),e1te2t)

########################################
# Read DATA 
########################################
#------------------------------------------------------------------------------------------------------------------------
str_month=''

Ar_size=(time_dim,lon.shape[0],lon.shape[1])
My_var1=npy.zeros(Ar_size)  ; My_var1SeasM=npy.zeros(Ar_size)   ; My_var1SeasS=npy.zeros(Ar_size) 
        
if AW_Tmax_maps or FWC_maps or TSD_maps or ATL_maps :
	Ar_size=(time_dim,z.shape[0],lon.shape[0],lon.shape[1])
	My_var1T=npy.zeros(Ar_size)   ;   My_var1S=npy.zeros(Ar_size)

	# For SSH maps
	Ar_size=(time_dim,lon.shape[0],lon.shape[1])
	My_var1ssh=npy.zeros(Ar_size)  

if MOC_maps :
	Ar_size=(time_dim,z.shape[0],lon.shape[0])
	My_MOC=npy.zeros(Ar_size)   

if ICE_maps :
	Ar_size=(time_dim,lon.shape[0],lon.shape[1])
	My_var1frld_SeasM=npy.zeros(Ar_size)   ; My_var1frld_SeasS=npy.zeros(Ar_size)

if MTS_maps :
	Ar_size=(2,z.shape[0],lon.shape[0],lon.shape[1])
	My_varT=npy.zeros(Ar_size)   ; My_varS=npy.zeros(Ar_size)

time_dim=(e_year-s_year+1)*12
Ar_size=(time_dim,lon.shape[0],lon.shape[1])
Mdata_read=npy.zeros(Ar_size)
tur_month = 0
t_months=(npy.arange(12)*30.+15.)/365.   ;   start = 1
time_grid=npy.arange(s_year,e_year+2,1.,dtype=int)
newlocsx  = npy.array(time_grid,'f')
newlabelsx = npy.array(time_grid,'i')

# Read the whole time series
c_year=s_year
while c_year <= e_year:
	print
	print ' 	The concerned year :', c_year
	print
	locpath=data_dir+'/'+str(c_year)+'/'+xiosfreq+'/'

	#########################################################################################################################################
        if FWC_maps or ATL_maps :
             print ' 			Read SSH variable'

	     locfile=CONFCASE+'_y'+str(c_year)+str_month+'.'+xiosfreq+'_gridT.nc'

	     # Read annual mean 
	     if chkfile(locpath+locfile) :
	     		My_var1ssh[c_year-s_year,:,:] = PyRaf.readfullNC4(locpath+locfile,'ssh')
			My_var1ssh[c_year-s_year,:,:]=npy.ma.where(My_var1ssh[c_year-s_year,:,:] >= 1e15, 0. , My_var1ssh[c_year-s_year,:,:])
             		My_sshmean1 = npy.sum(My_var1ssh[c_year-s_year,:,:]*e1te2t[:,:])/npy.sum(e1te2t)
	     		My_var1ssh[c_year-s_year,:,:] = My_var1ssh[c_year-s_year,:,:] - My_sshmean1
			if c_year == e_year : 
	     			ztmask2D=tmask[0,:,:].copy()
	     			tmask2Dtime=npy.tile(ztmask2D,(e_year-s_year+1,1,1))
	     			amask2Dtime = tmask2Dtime.copy()
				My_var1ssh=npy.ma.masked_where((amask2Dtime == 0),My_var1ssh).squeeze()

	#########################################################################################################################################
        if ICE_maps :

             print ' 			Read Ice volume variable '

	     zMyvar='sivolu'
	     locfile=CONFCASE+'_y'+str(c_year)+str_month+'.'+xiosfreq+'_icemod.nc'
	     # Read annual mean Ice thickness
	     if chkfile(locpath+locfile)   :      
			My_var1[c_year-s_year,:,:] = PyRaf.readfullNC4(locpath+locfile,zMyvar)
			if c_year == e_year : 
	     			ztmask2D=tmask[0,:,:].copy()
	     			tmask2Dtime=npy.tile(ztmask2D,(e_year-s_year+1,1,1))
	     			amask2Dtime = tmask2Dtime.copy()
				My_var1=npy.ma.masked_where((amask2Dtime == 0),My_var1).squeeze()
             else: 
			My_var1[c_year-s_year,:,:] = npy.nan

             print ' 			Read Sea-Ice fraction variable '

	     zMyvar='siconc'
	     # Read March/September mean 
	     locfilem3=CONFCASE+'_y'+str(c_year)+'m03.'+xiosfreq+'_icemod.nc'
	     locfilem9=CONFCASE+'_y'+str(c_year)+'m09.'+xiosfreq+'_icemod.nc'
	     if chkfile(locpath+locfilem3) : 
		My_var1frld_SeasM[c_year-s_year,:,:] = PyRaf.readnc_2d_tNC4(locpath+locfilem3,zMyvar,0)
	     else:
		My_var1frld_SeasM[c_year-s_year,:,:] = npy.nan

	     if chkfile(locpath+locfilem9) : 
		My_var1frld_SeasS[c_year-s_year,:,:] = PyRaf.readnc_2d_tNC4(locpath+locfilem9,zMyvar,0)
	     else:
		My_var1frld_SeasS[c_year-s_year,:,:] = npy.nan
	     if c_year == e_year : 
	     		ztmask2D=tmask[0,:,:].copy()
	     		tmask2Dtime=npy.tile(ztmask2D,(e_year-s_year+1,1,1))
	     		amask2Dtime = tmask2Dtime.copy()
			My_var1frld_SeasM=npy.ma.masked_where((amask2Dtime == 0),My_var1frld_SeasM).squeeze()
			My_var1frld_SeasS=npy.ma.masked_where((amask2Dtime == 0),My_var1frld_SeasS).squeeze()

	#########################################################################################################################################
	if MLD_maps or MTS_maps:

		zlocfilem3=CONFCASE+'_y'+str(c_year)+'m03.'+xiosfreq+'_gridT.nc'
		zlocfilem9=CONFCASE+'_y'+str(c_year)+'m09.'+xiosfreq+'_gridT.nc'
		zMyvar = 'mldr10_1'

             	print ' 			Read Mixed layer depth variable '
		if chkfile(locpath+zlocfilem3) : My_var1SeasM[c_year-s_year,:,:] = PyRaf.readnc_2d_tNC4(locpath+zlocfilem3,zMyvar,0)
		if chkfile(locpath+zlocfilem9) : My_var1SeasS[c_year-s_year,:,:] = PyRaf.readnc_2d_tNC4(locpath+zlocfilem9,zMyvar,0)
		if c_year == e_year :
	     		ztmask2D=tmask[0,:,:].copy()
	     		tmask2Dtime=npy.tile(ztmask2D,(e_year-s_year+1,1,1))
	     		amask2Dtime = tmask2Dtime.copy()

	        	My_var1SeasM=npy.ma.masked_where((amask2Dtime == 0),My_var1SeasM).squeeze()
	        	My_var1SeasS=npy.ma.masked_where((amask2Dtime == 0),My_var1SeasS).squeeze()

	#########################################################################################################################################
	# MLD Time-series in the K1 mooring 
	if MLD_maps or ATL_maps :

		zMyvar = 'mldr10_1'
		cur_month=0
		while cur_month <= 11 :
		        str_month='m0'+str(cur_month+1) if cur_month <= 8 else 'm'+str(cur_month+1)
		        locpath=data_dir+'/'+str(c_year)+'/'+xiosfreq+'/'
		        locfile=CONFIG+'-'+CASE+'_y'+str(c_year)+str_month+'.'+xiosfreq+'_gridT.nc'
		        if chkfile(locpath+locfile) :
		                field = Dataset(locpath+locfile)
		                Mdata_read[tur_month,:,:] = npy.squeeze(field.variables[zMyvar])
		        cur_month = cur_month + 1
		        tur_month = tur_month + 1

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
             print ' 			Read PSI variable '

	     # Read annual mean 
	     zlocfile=CONFCASE+'_y'+str(c_year)+'.'+xiosfreq+'_PSI.nc'
	     if chkfile(locpath+zlocfile) : My_var1[c_year-s_year,:,:] = PyRaf.readfullNC4(locpath+zlocfile,zMyvar)
	     if c_year == e_year :
	     		fmask2Dtime=npy.tile(fmask2D,(e_year-s_year+1,1,1))
	     		amask2Dtime = fmask2Dtime.copy()
			My_var1=npy.ma.masked_where((amask2Dtime == 0),My_var1).squeeze()

	     zMyvar = 'voeke' 
             print ' 			Read EKE variable @ depth ', z[0], ' m & ', z[23], ' m' 
	     
	     # Read annual mean
	     zlocfile=CONFCASE+'_y'+str(c_year)+str_month+'.'+xiosfreq+'_EKE.nc'
	     if chkfile(locpath+zlocfile) : My_var1SeasM[c_year-s_year,:,:] = PyRaf.readnc_3d_tlevNC4(locpath+zlocfile,zMyvar,0,0 )
	     if chkfile(locpath+zlocfile) : My_var1SeasS[c_year-s_year,:,:] = PyRaf.readnc_3d_tlevNC4(locpath+zlocfile,zMyvar,0,23)
	     if c_year == e_year :
			tmask2Dtime=npy.tile(tmask2D,(e_year-s_year+1,1,1))
			amask2Dtime = tmask2Dtime.copy()
			My_var1SeasM=npy.ma.masked_where((amask2Dtime == 0),My_var1SeasM).squeeze()
			
			ztmask2D=tmask[23,:,:].copy()
			tmask2Dtime=npy.tile(ztmask2D,(e_year-s_year+1,1,1))
			amask2Dtime = tmask2Dtime.copy()
			My_var1SeasS=npy.ma.masked_where((amask2Dtime == 0),My_var1SeasS).squeeze()


	#########################################################################################################################################
	if ( AW_Tmax_maps or FWC_maps or TSD_maps or ATL_maps ) :

		   print ' 			Read full T/S 3D spatial variables to compute the AWTmax or FWC '

		   # Read annual mean 
	           locfile=CONFCASE+'_y'+str(c_year)+str_month+'.'+xiosfreq+'_gridT.nc'
		   if chkfile(locpath+locfile) : My_var1T[c_year-s_year,:,:,:] = PyRaf.readfullNC4(locpath+locfile,'votemper')
		   if chkfile(locpath+locfile) : My_var1S[c_year-s_year,:,:,:] = PyRaf.readfullNC4(locpath+locfile,'vosaline')
		   if c_year == e_year :
	           		tmask3Dtime=npy.tile(tmask,(e_year-s_year+1,1,1,1))
	           		My_var1T=npy.ma.masked_where((tmask3Dtime == 0),My_var1T).squeeze() 
	           		My_var1S=npy.ma.masked_where((tmask3Dtime == 0),My_var1S).squeeze()

	#########################################################################################################################################
	if ( MTS_maps ) :

		   print ' 			Read full T/S 3D spatial variables to compute the AWTmax or FWC '

		   zlocfilem3=CONFCASE+'_y'+str(c_year)+'m03.'+xiosfreq+'_gridT.nc'
		   zlocfilem9=CONFCASE+'_y'+str(c_year)+'m09.'+xiosfreq+'_gridT.nc'
		   # Read annual mean 
		   if chkfile(locpath+zlocfilem3) : 
			My_varT[0,:,:,:] = PyRaf.readfullNC4(locpath+zlocfilem3,'votemper')
			My_varS[0,:,:,:] = PyRaf.readfullNC4(locpath+zlocfilem3,'vosaline')
		   if chkfile(locpath+zlocfilem9) :
			My_varT[1,:,:,:] = PyRaf.readfullNC4(locpath+zlocfilem9,'votemper')
			My_varS[1,:,:,:] = PyRaf.readfullNC4(locpath+zlocfilem9,'vosaline')
		   if c_year == e_year :
	           		tmask3Dtime=npy.tile(tmask,(2,1,1,1))
	           		My_varT=npy.ma.masked_where((tmask3Dtime == 0),My_varT).squeeze() 
	           		My_varS=npy.ma.masked_where((tmask3Dtime == 0),My_varS).squeeze()

	#########################################################################################################################################
	if MOC_maps :

		   print ' 			Read full MOC 2D spatial variables '

		   # Read annual mean 
	           locfile=CONFCASE+'_y'+str(c_year)+str_month+'.'+xiosfreq+'_MOC.nc'
		   if chkfile(locpath+locfile) : 
			fieldmoc=Dataset(locpath+locfile)
			My_MOC[c_year-s_year,:,:] = npy.array(npy.squeeze(fieldmoc.variables['zomsfglo']))
	
	c_year = c_year + 1 

# Read initial state to compare with
if AW_Tmax_maps or FWC_maps or TSD_maps or ATL_maps :
   Ar_size=(z.shape[0],lon.shape[0],lon.shape[1])
   My_varTinit=npy.zeros(Ar_size)   ; My_varSinit=npy.zeros(Ar_size)
   print ' 			Read initial state  '
   locpath=data_dir+'/'
   locfile=CONFCASE+'_init_gridT.nc'
   if chkfile(locpath+locfile) : 
	My_varTinit[:,:,:] = PyRaf.readfullNC4(locpath+locfile,'votemper')
	My_varSinit[:,:,:] = PyRaf.readfullNC4(locpath+locfile,'vosaline')

# If many years have been read, compute the mean over the corresponding period
if e_year-s_year+1 > 1 :
	if AW_Tmax_maps or FWC_maps  or ATL_maps:
		My_var1T=npy.mean(My_var1T,axis=0)   
		My_var1S=npy.mean(My_var1S,axis=0)
		My_var1ssh=npy.mean(My_var1ssh,axis=0)
	
	if ICE_maps :
		My_var1=npy.mean(My_var1,axis=0)
		My_var1frld_SeasM=npy.mean(My_var1frld_SeasM,axis=0)
		My_var1frld_SeasS=npy.mean(My_var1frld_SeasS,axis=0)
	
	if MLD_maps or ATL_maps :
		My_var1SeasM=npy.mean(My_var1SeasM,axis=0)
		My_var1SeasS=npy.mean(My_var1SeasS,axis=0)
	
	if DYN_maps :
		My_var1=npy.mean(My_var1,axis=0)  
		My_var1SeasM=npy.mean(My_var1SeasM,axis=0)   
		My_var1SeasS=npy.mean(My_var1SeasS,axis=0) 

	if MOC_maps :
		My_MOC=npy.mean(My_MOC,axis=0)  

########################################
# Call function to make plots
########################################
#------------------------------------------------------------------------------------------------------------------------
# In the following call the appropriate function depending the diagnostic

# To plot the mean T/S in the ML
if MTS_maps : CREG_maps_func.MTS_maps(lon, lat, CONFIG, CASE, My_var1SeasM, My_var1SeasS, My_varT, My_varS, gdept_0, ze3, climyear, data_dir, grid_dir )

# To plot the Atlantic Water maximum temperature as the associated depth
# Use the salinity criteria S < 33.5
if AW_Tmax_maps : 
	zMyvar='votemper'
	CREG_maps_func.AWTmax_maps( lon, lat, My_var1T, My_var1S, z, zMyvar, CONFIG, CASE, climyear)

# To plot SSH and FWC (based on a salinity ref of 34.8 PSU)
if FWC_maps : 
	CREG_maps_func.FWC_maps( lon, lat, My_var1S, My_varSinit, My_var1ssh, CONFIG, CASE, climyear, ze3, tmask, time_dim )

# To plot ICE variables
if ICE_maps : 
	num_fram=320
	# Annual mean Ice thickness
	zMyvar='sivolu'   ; fram=num_fram+1
	My_var1=npy.ma.masked_where(My_var1 == 0., My_var1).squeeze()
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
			fieldmask=Dataset(locpath+locfile)
			lon = npy.array(fieldmask.variables['nav_lon'])
			lat = npy.array(fieldmask.variables['nav_lat'])
			zmask= npy.squeeze(fieldmask.variables['tmask'])
	else: 
		zmask = tmask.copy()

        zMyvar='sivolu'   ; fram=num_fram+2
	obs_thick=CREG_maps_func.ICE_THICK_OBS(zconfig=CONFIG,t_year=s_year)
	obs_thick=npy.ma.masked_where(npy.squeeze(zmask[0,:,:]) == 0., obs_thick)
	obs_thick=npy.ma.masked_where(obs_thick == 0., obs_thick)
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
	plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.pdf')

        if NCDF_OUT:
		# ICE fields
		#######################
		cmd_ddate="date"  ;  get_output = subprocess.check_output(cmd_ddate)
                nc_f = './NETCDF/'+CONFIG+'-'+CASE+'_ICEClim_'+'y'+climyear+'.nc'
                w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
                w_nc_fid.description = "Diagnostics has been calculated using the Arctic monitoring tool "
		w_nc_fid.date=get_output.decode("utf-8")
                w_nc_fid.createDimension('x', My_var1.shape[1])
                w_nc_fid.createDimension('y', My_var1.shape[0])
                w_nc_fid.createDimension('xobs', obs_lat.shape[1])
                w_nc_fid.createDimension('yobs', obs_lat.shape[0])
                #w_nc_fid.createDimension('time_counter', None)

                w_nc_var = w_nc_fid.createVariable('IceThick_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Model annual mean Ice thickness'
                w_nc_var.units="m"
                w_nc_fid.variables['IceThick_mod'][:,:] = My_var1

                w_nc_var = w_nc_fid.createVariable('IceConceM03_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Model monthly mean Ice concentration in march'
                w_nc_var.units="-"
                w_nc_fid.variables['IceConceM03_mod'][:,:] = My_var1frld_SeasM

                w_nc_var = w_nc_fid.createVariable('IceConceM09_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Model monthly mean Ice concentration in september'
                w_nc_var.units="-"
                w_nc_fid.variables['IceConceM09_mod'][:,:] = My_var1frld_SeasS

                w_nc_var = w_nc_fid.createVariable('IceThick_obs', 'f4', ('y','x'))
		if s_year >= 1979 and s_year <= 2013 :
			w_nc_var.long_name='PIOMAS annual mean Ice thickness over '+str(s_year)
		else :
			w_nc_var.long_name='PIOMAS climatological mean Ice thickness over 1979-2013'
                w_nc_var.units="m"
                w_nc_fid.variables['IceThick_obs'][:,:] = obs_thick

                w_nc_var = w_nc_fid.createVariable('IceConceM03_obs', 'f4', ('yobs','xobs'))
		if s_year >= 1979 and s_year <= 2015 :
			w_nc_var.long_name='NSDIC monthly mean Ice concentration in march '+str(s_year)
		else :
			w_nc_var.long_name='NSDIC climatological mean Ice concentration in march over 1979-2015'
                w_nc_var.units="-"
                w_nc_fid.variables['IceConceM03_obs'][:,:] = obs_conc_m03

                w_nc_var = w_nc_fid.createVariable('IceConceM09_obs', 'f4', ('yobs','xobs'))
		if s_year >= 1979 and s_year <= 2015 :
			w_nc_var.long_name='NSDIC monthly mean Ice concentration in september '+str(s_year)
		else :
			w_nc_var.long_name='NSDIC climatological mean Ice concentration in september over 1979-2015'
                w_nc_var.units="-"
                w_nc_fid.variables['IceConceM09_obs'][:,:] = obs_conc_m09


                w_nc_var = w_nc_fid.createVariable('lat_obs', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='Degrees north'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lat_obs'][:,:] = obs_lat

                w_nc_var = w_nc_fid.createVariable('lon_obs', 'f4', ('yobs','xobs'))
                w_nc_var.long_name='Degrees east'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lon_obs'][:,:] = obs_lon
                
                w_nc_var = w_nc_fid.createVariable('lat_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Degrees north'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lat_mod'][:,:] = lat

                w_nc_var = w_nc_fid.createVariable('lon_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Degrees east'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lon_mod'][:,:] = lon
                
                w_nc_fid.close()  # close the file

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
	plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.pdf')

        if NCDF_OUT:
		# MLD fields
		#######################
		cmd_ddate="date"  ;  get_output = subprocess.check_output(cmd_ddate)
                nc_f = './NETCDF/'+CONFIG+'-'+CASE+'_MLDClim_'+'y'+climyear+'.nc'
                w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
                w_nc_fid.description = "Diagnostics has been calculated using the Arctic monitoring tool "
		w_nc_fid.date=get_output.decode("utf-8")
                w_nc_fid.createDimension('x', My_var1SeasM.shape[1])
                w_nc_fid.createDimension('y', My_var1SeasM.shape[0])
                w_nc_fid.createDimension('xobs', lat_obs.shape[1])
                w_nc_fid.createDimension('yobs', lat_obs.shape[0])
                #w_nc_fid.createDimension('time_counter', None)

                w_nc_var = w_nc_fid.createVariable('MLDd001M03_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Model monthly mean MLD in march based on a density criterai of 0.01 kg/m^3'
                w_nc_var.units="m"
                w_nc_fid.variables['MLDd001M03_mod'][:,:] = My_var1SeasM

                w_nc_var = w_nc_fid.createVariable('MLDd001M09_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Model monthly mean MLD in september based on a density criterai of 0.01 kg/m^3'
                w_nc_var.units="m"
                w_nc_fid.variables['MLDd001M09_mod'][:,:] = My_var1SeasS

                w_nc_var = w_nc_fid.createVariable('MLDd001M03_obs', 'f4', ('yobs','xobs'))
		w_nc_var.long_name='MIMOC climatological mean in march over 2003-2014 based on a density criterai of 0.01 kg/m^3'
                w_nc_var.units="m"
                w_nc_fid.variables['MLDd001M03_obs'][:,:] = mld_obs_m03

                w_nc_var = w_nc_fid.createVariable('MLDd001M09_obs', 'f4', ('yobs','xobs'))
		w_nc_var.long_name='MIMOC climatological mean in september over 2003-2014 based on a density criterai of 0.01 kg/m^3'
                w_nc_var.units="m"
                w_nc_fid.variables['MLDd001M09_obs'][:,:] = mld_obs_m09

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
                w_nc_fid.variables['lat_mod'][:,:] = lat

                w_nc_var = w_nc_fid.createVariable('lon_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Degrees east'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lon_mod'][:,:] = lon
                
                w_nc_fid.close()  # close the file


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
	plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.pdf')

        if NCDF_OUT:
		# DYN fields
		#######################
		cmd_ddate="date"  ;  get_output = subprocess.check_output(cmd_ddate)
                nc_f = './NETCDF/'+CONFIG+'-'+CASE+'_DYNClim_'+'y'+climyear+'.nc'
                w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
                w_nc_fid.description = "Diagnostics has been calculated using the Arctic monitoring tool "
		w_nc_fid.date=get_output.decode("utf-8")
                w_nc_fid.createDimension('x', My_var1.shape[1])
                w_nc_fid.createDimension('y', My_var1.shape[0])
                w_nc_fid.createDimension('xobs', lat_obs.shape[1])
                w_nc_fid.createDimension('yobs', lat_obs.shape[0])
                #w_nc_fid.createDimension('time_counter', None)

                w_nc_var = w_nc_fid.createVariable('PSI_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Model annual mean barotropic streamfunction '
                w_nc_var.units="Sv"
                w_nc_fid.variables['PSI_mod'][:,:] = My_var1

                w_nc_var = w_nc_fid.createVariable('EKESurf_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Model annual mean EKE at the surface'
                w_nc_var.units="-"
                w_nc_fid.variables['EKESurf_mod'][:,:] = My_var1SeasM

                w_nc_var = w_nc_fid.createVariable('EKEz100_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Model annual mean EKE @ ~100m depth'
                w_nc_var.units="cm2/s2"
                w_nc_fid.variables['EKEz100_mod'][:,:] = My_var1SeasS

                w_nc_var = w_nc_fid.createVariable('EKESurf_obs', 'f4', ('yobs','xobs'))
		if s_year >= 2003 and s_year <= 2014 :
			w_nc_var.long_name='EKE annual mean derived from DOT field (Armitage et al. 2017) in '+str(s_year)
		else :
			w_nc_var.long_name='EKE climatological mean derived from DOT field (Armitage et al. 2017) over 2003-2014'
                w_nc_var.units="cm2/s2"
                w_nc_fid.variables['EKESurf_obs'][:,:] = obs_eke

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
                w_nc_fid.variables['lat_mod'][:,:] = lat

                w_nc_var = w_nc_fid.createVariable('lon_mod', 'f4', ('y','x'))
                w_nc_var.long_name='Degrees east'
                w_nc_var.units="Deg"
                w_nc_fid.variables['lon_mod'][:,:] = lon
                
                w_nc_fid.close()  # close the file


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
	plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.pdf')


# To plot ATL variables
if ATL_maps : 
	# MLD IN THE LABRADOR SEA IN MARCH
	###################################
	##	plt.figure()
	##	num_fram=110
	##	# March mean MLD
	##	zMyvar='mldr10_1'   ; fram=num_fram+1
	##	my_cblab=r'MLD (m)'   ;   my_cmap=plt.cm.get_cmap('Blues')
	##	ztitle=CASE +' mean MLD01 over \n'+climyear+'  m03'
	##	vmin=0. ; vmax=2400. ; vint=100.   ;   contours=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]
	##	limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]

	##	plt.subplot(fram)
	##	zoutmap=CREG_maps_func.Atl_Bat(ztype='isol1000',zarea='labsea')
	##	CREG_maps_func.Atl_plot(lon, lat, My_var1SeasM, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='labsea')

	##	zfile_ext='_LAB_MLDClimm03_'
        ##	plt.tight_layout()
	##	plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.pdf')

	# Plot the Time-serie for MLD at a specific location K1 mooring in the Labrador Sea and in Irminger Sea
	# After Schott et al. DSRI2009 56.33N, -52.40W
	plt.clf()
	plt.figure()
	i_K1=518   ;   j_K1=502   # CREG12.L75 C-type indices
	ax=plt.subplot(211)
	# In Lab. Sea
	plt.plot(time_axis,-1.*npy.squeeze(Mdata_read[:,j_K1:j_K1+1,i_K1:i_K1+1]),linewidth=0.8, color='k', label='Lab Sea K1')
	# Plot obs. MLD in March
	year_obs=npy.arange(1995,2006,1)+0.20547945   ; mld_obs=[-2300.,-1300.,-1400.,-1000.,-1000.,-1100.,-1100.,-1200.,-1400.,-700.,-1300.]
	plt.scatter(year_obs,mld_obs)

	# In Irm. Sea
	i_K1=697   ;   j_K1=577   # CREG12.L75 C-type indices
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
	plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.pdf')

	# Add an artificial mooring within the deepest convection area 
	# -54W 58N
	#  dl_dis=    1.634 km
	#      507       507       541       541
	# -54.0272  -54.0272   57.9970   57.9970
	plt.clf()
	plt.figure()
	i_K1=506   ;   j_K1=540   # CREG12.L75 C-type indices
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
	plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.pdf')


	###	plt.clf()
	###	plt.figure()
	###	# MLD IN THE IRMINGER SEA IN MARCH
	###	###################################
	###	num_fram=110
	###	# March mean MLD
	###	zMyvar='mldr10_1'   ; fram=num_fram+1
	###	my_cblab=r'MLD (m)'   ;   my_cmap=plt.cm.get_cmap('Blues')
	###	ztitle=CASE +' mean MLD01 over \n'+climyear+'  m03'
	###	vmin=0. ; vmax=1600. ; vint=100.   ;   contours=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
	###	limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
	###	plt.subplot(fram)
	###	zoutmap=CREG_maps_func.Atl_Bat(ztype='isol1000',zarea='irmsea')
	###	CREG_maps_func.Atl_plot(lon, lat, My_var1SeasM, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='irmsea')

	###	zfile_ext='_IRM_MLDClimm03_'
        ###	plt.tight_layout()
	###	plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.pdf')

	###	plt.clf()
	###	plt.figure()
	###	# MLD IN THE GIN SEAS IN MARCH
	###	###################################
	###	num_fram=110
	###	# March mean MLD
	###	zMyvar='mldr10_1'   ; fram=num_fram+1
	###	my_cblab=r'MLD (m)'   ;   my_cmap=plt.cm.get_cmap('Blues')
	###	ztitle=CASE +' mean MLD01 over \n'+climyear+'  m03'
	###	vmin=0. ; vmax=1600. ; vint=100.   ;   contours=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
	###	limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,600.,800.,1000.,1200.,1600.]
	###	plt.subplot(fram)
	###	zoutmap=CREG_maps_func.Atl_Bat(ztype='isol1000',zarea='ginsea')
	###	CREG_maps_func.Atl_plot(lon, lat, My_var1SeasM, contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='ginsea')

	###	zfile_ext='_GIN_MLDClimm03_'
        ###	plt.tight_layout()
	###	plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.pdf')

	# PLOT ISOTHERME 17 Deg off CAPE HATTERAS
	#########################################
	##	num_fram=110
	##	zMyvar='votemper'   ; fram=num_fram+1
	##	my_cblab=r'ISO 17 (DegC)'   ;   my_cmap=plt.cm.get_cmap('jet')
	##	ztitle=CASE +' mean Iso 17 DegC over \n'+climyear
	##	vmin=0. ; vmax=2400. ; vint=100.   ;   contours=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]
	##	limits=[vmin,vmax,vint]  ;              myticks=[0.,100.,200.,400.,800.,1200.,1600.,2000.,2400.]

	##	plt.subplot(fram)
	##	klev=29
	##	zoutmap=CREG_maps_func.Atl_Bat(ztype='isol1000', zarea='GulfS')
	##	CREG_maps_func.Atl_plot(lon, lat, npy.squeeze(My_var1T[klev,:,:])   , contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='GulfS')
	##	CREG_maps_func.Atl_plot(lon, lat, npy.squeeze(My_varTinit[klev,:,:]), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='GulfS', data_ref=True)

	##	zfile_ext='_ISO17Clim_'
        ##	plt.tight_layout()
	##	plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.pdf')

	# PLOT SSH OVER THE ATLANTIC AREA
	#########################################
	##	num_fram=110
	##	zMyvar='sossheig'   ; fram=num_fram+1
	##	my_cblab=r'SSH (cm)'   ;   my_cmap=plt.cm.get_cmap('coolwarm')
	##	ztitle=CASE +' mean SSH over \n'+climyear
	##	vmin=-100. ; vmax=100. ; vint=5.  ;   contours=npy.arange(vmin,vmax+vint,vint)
	##	limits=[vmin,vmax,vint]           ;   myticks=npy.arange(vmin,vmax+vint,vint)

	##	plt.figure()
	##	#plt.subplot(fram)
	##	zoutmap=CREG_maps_func.Atl_Bat(ztype='isol1000')
	##	CREG_maps_func.Atl_plot(lon, lat, My_var1ssh*100.  , contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar)
	##	#CREG_maps_func.Atl_plot(lon, lat, npy.squeeze(My_varTinit[klev,:,:]), contours, limits, myticks, name=ztitle, zmy_cblab=my_cblab, zmy_cmap=my_cmap, zvar=zMyvar, zarea='GulfS', data_ref=True)

	##	zfile_ext='_SSHClim_'
        ##	#plt.tight_layout()
	##	plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.pdf')


# To plot AMOC 
if MOC_maps : 
	plt_AMOCTS=True

	if not plt_AMOCTS: 
		# AAMOC 
		#######
		plt.figure()
		num_fram=210
		fram=num_fram+1
		my_cblab=r'AMOC (Sv)'   ;   my_cmap=plt.cm.get_cmap('jet')
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
			fieldmask=Dataset(locpath+locfile)
			full_lat = npy.array(fieldmask.variables['nav_lat'])
    		select_ylat=full_lat[:,862]
		select_ylat_reshape=npy.reshape(select_ylat,(select_ylat.size,1))
		ypltz=npy.repeat(select_ylat_reshape,z.shape[0],axis=1).T

		z2dt=npy.reshape(z,(z.size,1))
		zplt = npy.repeat(z2dt,lat.shape[0],axis=1)

		plt.title(ztitle,fontsize=6)
		plt.contourf(ypltz,zplt*(-1.e-3),npy.squeeze(My_MOC),contours,cmap=my_cmap,norm=norm,extend='both')
		contours=npy.arange(vmin,vmax+vint,2.*vint)
		C=plt.contour(ypltz,zplt*(-1.e-3),npy.squeeze(My_MOC),linewidths=0.5,levels=contours, colors='k')
		plt.clabel(C, C.levels, inline=True, fmt='%3.0f', fontsize=6)

		zfile_ext='_AAMOCClim_'
        	#plt.tight_layout()
		plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.pdf')


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
		plt.plot( time_grid, npy.max(My_MOC[:,:,jloc],axis=1)   , 'g',linewidth=1., label=CONFIG)
        	plt.title(CASE+' Max AAMOCz @ '+zlat+'\n '+str(s_year)+str(e_year),size=9)
        	plt.ylabel('Max AAMOCz  \n'+r'(Sv)', size=7)
        	plt.ylim([10.,15.])
        	plt.xticks(newlocsx,newlabelsx,size=5)
        	plt.setp(ax.get_xticklabels(),rotation=90)
        	plt.yticks(size=6)
        	plt.grid(True, linewidth=0.7,linestyle='--',alpha=0.7,color='grey')

		ax=plt.subplot(212)
		plt.plot( time_grid, npy.max(My_MOC[:,:,:],axis=(1,2))   , 'k',linewidth=1., label=CONFIG)
        	#plt.title(CASE+' Max AAMOCz '+'\n '+str(s_year)+str(e_year),size=9)
        	plt.ylabel('Global Max AAMOCz  \n'+r'(Sv)', size=7)
        	plt.ylim([10.,20.])
        	plt.xticks(newlocsx,newlabelsx,size=5)
        	plt.setp(ax.get_xticklabels(),rotation=90)
        	plt.yticks(size=6)
        	plt.grid(True, linewidth=0.7,linestyle='--',alpha=0.7,color='grey')


		zfile_ext='_MaxAAMOCz_'+zlat+'_TiSe_'
        	#plt.tight_layout()
		plt.savefig(CONFIG+'-'+CASE+zfile_ext+'y'+climyear+'.pdf')
