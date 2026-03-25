#!/usr/bin/env python

import sys 
import matplotlib
matplotlib.use('Agg')
import numpy as npy
from checkfile import *
from CREG_sections_func import *
import subprocess
import xarray as xr 
from fsspec.implementations.local import LocalFileSystem
fs = LocalFileSystem()

s_year=XXSYEAXX
e_year=XXEYEAXX
lgTS_ys=XXLGTSSXX
lgTS_ye=XXLGTSEXX
xiosfreq=XXXIOSFREQXX
NCDF_OUT=XXNCDFOUTXX
main_dir='./'
CONFIG='XXCONFXX'   ; CASE='XXCASEXX'     ;  CASE2='XXCASE2XX'
CONFCASE=CONFIG+'-'+CASE   ; CONFCASE2=CONFIG+'-'+CASE2
data_dir=main_dir+CONFIG+'/'+CONFCASE+'-MEAN/'	  ; data_dir2=main_dir+CONFIG+'/'+CONFCASE2+'-MEAN/'
grid_dir=main_dir+CONFIG+'/GRID/'

DIR_FIG_OUT='./'

# Infos concernant les climatologies sur la periode de la simulation
if e_year-s_year == 0 :
	climyear=str(s_year)
else:
	climyear=str(s_year)+str(e_year)

print()
print()
print('				       Configuration :' + CONFCASE)
print('				       Period	     :' + str(s_year),' - ',str(e_year))
print()

#####################################################################
var_temp={'name':"votemper",'units':u"degC",'_FillValue': 0.,'fext':"T",'igrd':1 ,'ze3':"e3t_0",'gdep':"gdept",'long_name':"Temperature"}
var_sali={'name':"vosaline",'units':u"psu" ,'_FillValue': 0.,'fext':"S",'igrd':1 ,'ze3':"e3t_0",'gdep':"gdept",'long_name':"Salinity"}
var_crtx={'name':"vozocrtx",'units':u"m/s" ,'_FillValue': 0.,'fext':"U",'igrd':3 ,'ze3':"e3u_0",'gdep':"gdepu",'long_name':"Zonal velocity"}
var_crty={'name':"vomecrty",'units':u"m/s" ,'_FillValue': 0.,'fext':"V",'igrd':3 ,'ze3':"e3v_0",'gdep':"gdepv",'long_name':"Meridional velocity"}
var_avtz={'name':"votkeavt",'units':u"m2/s" ,'_FillValue': 0.,'fext':"W",'igrd':3 ,'ze3':"e3w_0",'gdep':"gdepw",'long_name':"Vertical diffusivity"}
var_hthi={'name':"sivolu"  ,'units':u"m"   ,'_FillValue': 0.,'fext':"I",'igrd':1 ,'ze3':"e3t_0",'gdep':"gdept",'long_name':"Ice thickness "}
var_vehi={'name':"sivelv"  ,'units':u"m/s" ,'_FillValue': 0.,'fext':"I",'igrd':3 ,'ze3':"e3v_0",'gdep':"gdepv",'long_name':"Ice meridional velocity "}
var_dens={'name':"rhop_sig0",'units':u"kg/m3"	,'_FillValue': 0.,'fext':"D",'igrd':1 ,'ze3':"e3t_0",'gdep':"gdept",'long_name':"Potential density "}
########################################

# Initialize dictionnaries
Bering={'name':"Bering"}   ;   SouthG={'name':"SouthG"}   ;   FramS={'name':"FramS"}	 ;   FramObs={'name':"FramObs"}       ;   Davis={'name':"Davis"}
Beauf={'name':"Beaufort"}  ;   ArcAn={'name':"ArcAnna"}   ;   Baren={'name':"Barents"}	 ;   Kara={'name':"Kara"}   

########################################
# Read GRID 
########################################
#------------------------------------------------------------------------------------------------------------------------

# Read only once for the same section type
locpath=grid_dir
locfile=CONFCASE+'_mask.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	ds_grd = xr.open_dataset(locpath+locfile)[['glamt','gphit','tmask','umask','vmask']]
	lon = ds_grd['glamt'].squeeze()
	lat = ds_grd['gphit'].squeeze()
	tmask = ds_grd['tmask'].isel(time_counter=0)
	tmask = tmask.rename({'nav_lev':'z'}) 
	umask = ds_grd['umask'].isel(time_counter=0)
	umask = umask.rename({'nav_lev':'z'}) 
	vmask = ds_grd['vmask'].isel(time_counter=0)
	vmask = vmask.rename({'nav_lev':'z'}) 

# Read appropriate vertical scale factor 
locpath=grid_dir
locfile=CONFCASE+'_mesh_zgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	infieldV = var_crty   ;   infieldT = var_temp
	ds_zgr = xr.open_dataset(locpath+locfile)[[infieldT['ze3'],infieldV['ze3'],'gdept_0']]
	e3t = ds_zgr[infieldT['ze3']].squeeze()
	e3t = e3t.rename({'nav_lev':'z'})
	e3v = ds_zgr[infieldV['ze3']].squeeze()
	e3v = e3v.rename({'nav_lev':'z'})
	gdept_0 = ds_zgr['gdept_0'].squeeze()
	gdept_0 = gdept_0.rename({'nav_lev':'z'})

locpath=grid_dir
locfile=CONFCASE+'_mesh_hgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	ds_mes = xr.open_dataset(locpath+locfile)[['e1t','e2t','e1v']]
	e1t = ds_mes['e1t'].squeeze()
	e1v = ds_mes['e1v'].squeeze()

#time_dim = (e_year-s_year+1)*12

#vmask3Dtime=npy.tile(vmask[:,:,:],(time_dim,1,1,1))

#------------------------------------------------------------------------------------------------------------------------
########################################
# Read DATA
########################################
#------------------------------------------------------------------------------------------------------------------------
print('				READ INPUT DATA FROM '+CONFIG+'-'+CASE)
print() 

for selsec in XXSECPLOTXX:
	if selsec['name'] == 'Beaufort': 
		infield = var_crtx
		print()
		print('						 >>>>>>>>>  The curent variable processed is :	      '+ infield['fext'])
		print()
		# List files to be read
		locpath = data_dir+'/'+str(s_year)+'/'+xiosfreq+'/'
		locfile = CONFCASE+'_y'+str(s_year)+'m*.'+xiosfreq+'_gridU.nc'
		UOCE_files = [f for f in fs.glob(locpath+locfile)]
		
		if len(UOCE_files) == 12 :
		   ds_Udata = xr.open_mfdataset(locpath+locfile, engine="netcdf4", concat_dim=["time_counter"], combine='nested', parallel=True)[[infield['name']]]
		   ds_Udata = ds_Udata.rename({'depthu':'z'})
		   Udata_read = ds_Udata[infield['name']].squeeze()

		   # Annual mean 
		   Udata_read_mean = Udata_read.mean('time_counter').squeeze()
		   
		   # Mask DataArray 
		   Udata_read_mean = xr.where( umask < 1., npy.nan, Udata_read_mean )

#########################################################################################################################################
infield = var_crty
print()
print('						 >>>>>>>>>  The curent variable processed is :	      '+ infield['fext'])
print()
# List files to be read
locpath = data_dir+'/'+str(s_year)+'/'+xiosfreq+'/'
locfile = CONFCASE+'_y'+str(s_year)+'m*.'+xiosfreq+'_gridV.nc'
VOCE_files = [f for f in fs.glob(locpath+locfile)]

if len(VOCE_files) == 12 :
   ds_Vdata = xr.open_mfdataset(locpath+locfile, engine="netcdf4", concat_dim=["time_counter"], combine='nested', parallel=True)[[infield['name']]]
   ds_Vdata = ds_Vdata.rename({'depthv':'z'})
   Vdata_read = ds_Vdata[infield['name']].squeeze()

   # Annual mean 
   Vdata_read_mean = Vdata_read.mean('time_counter').squeeze()
   
   # Mask DataArray 
   Vdata_read_mean = xr.where( vmask < 1., npy.nan, Vdata_read_mean )

#########################################################################################################################################
infieldT = var_temp   ;  infieldS = var_sali
print()
print('						 >>>>>>>>>  The curent variable processed is :	      '+ infieldT['fext']+' & '+infieldS['fext'])
print()
# List files to be read
locpath = data_dir+'/'+str(s_year)+'/'+xiosfreq+'/'
locfile = CONFCASE+'_y'+str(s_year)+'m*.'+xiosfreq+'_gridT.nc'
TS_files = [f for f in fs.glob(locpath+locfile)]

if len(TS_files) == 12 :
   ds_TSdata = xr.open_mfdataset(locpath+locfile, engine="netcdf4", concat_dim=["time_counter"], combine='nested', parallel=True)[[infieldT['name'],infieldS['name']]]
   ds_TSdata = ds_TSdata.rename({'deptht':'z'})
   Tdata_read = ds_TSdata[infieldT['name']].squeeze()
   Sdata_read = ds_TSdata[infieldS['name']].squeeze()

   # Annual mean 
   Tdata_read_mean = Tdata_read.mean('time_counter').squeeze()
   Sdata_read_mean = Sdata_read.mean('time_counter').squeeze()

   # Mask DataArray 
   Tdata_read_mean = xr.where( tmask < 1., npy.nan, Tdata_read_mean )
   Sdata_read_mean = xr.where( tmask < 1., npy.nan, Sdata_read_mean )

#########################################################################################################################################
infield = var_avtz
print()
print('						 >>>>>>>>>  The curent variable processed is :	      '+ infield['fext'])
print()
# List files to be read
locpath = data_dir+'/'+str(s_year)+'/'+xiosfreq+'/'
locfile = CONFCASE+'_y'+str(s_year)+'m*.'+xiosfreq+'_gridW.nc'
WOCE_files = [f for f in fs.glob(locpath+locfile)]

if len(VOCE_files) == 12 :
   ds_Wdata = xr.open_mfdataset(locpath+locfile, engine="netcdf4", concat_dim=["time_counter"], combine='nested', parallel=True)[[infield['name']]]
   ds_Wdata = ds_Wdata.rename({'depthw':'z'})
   Kdata_read = ds_Wdata[infield['name']].squeeze()

   # Annual mean
   Kdata_read_mean = Kdata_read.mean('time_counter').squeeze()

   # Mask DataArray 
   #Kdata_read_mean = xr.where( tmask < 1., npy.nan, Kdata_read_mean )

#########################################################################################################################################
infieldIH = var_hthi   ;   infieldIV = var_vehi
print()
print('						 >>>>>>>>>  The curent variable processed is :	      '+ infieldIH['fext']+ ' & '+infieldIV['fext'])
print()
# List files to be read
locpath = data_dir+'/'+str(s_year)+'/'+xiosfreq+'/'
locfile = CONFCASE+'_y'+str(s_year)+'m*.'+xiosfreq+'_icemod.nc'
ICE_files = [f for f in fs.glob(locpath+locfile)]

if len(ICE_files) == 12 :
   ds_Idata = xr.open_mfdataset(locpath+locfile, engine="netcdf4", concat_dim=["time_counter"], combine='nested', parallel=True)[[infieldIH['name'],infieldIV['name']]]
   IHdata_read = ds_Idata[infieldIH['name']].squeeze()
   IVdata_read = ds_Idata[infieldIV['name']].squeeze()

#########################################################################################################################################
print()
print('						 >>>>>>>>>  Read initial state 	     ')
print()
# Read initial state to compare with
locpath = data_dir+'/'
locfile = CONFCASE+'_init_gridT.nc'
if chkfile(locpath+locfile) : 
	ds_TSinit = xr.open_dataset(locpath+locfile, engine="netcdf4")[['votemper','vosaline']]
	ds_TSinit = ds_TSinit.rename({'nav_lev':'z'})
	My_varTinit = ds_TSinit['votemper'].squeeze()
	My_varSinit = ds_TSinit['vosaline'].squeeze()
	
	# Mask DataArray 
	My_varTinit_mean = xr.where( tmask < 1., npy.nan, My_varTinit )
	My_varSinit_mean = xr.where( tmask < 1., npy.nan, My_varSinit )

#------------------------------------------------------------------------------------------------------------------------
########################################
# Compute DIAGNOSTICS
########################################
#------------------------------------------------------------------------------------------------------------------------

print()
print('				#############################################################################  ')
print('				#############################################################################  ')
print('				###### COMPUTE VOLUME, TEMPERATURE & SALINITY TRANSPORT THROUGH SECTIONS ####  ')
print('				#############################################################################  ')
print('				#############################################################################  ')
print()

strait='' ; icount_str=0
for selsec in XXDIAGSSECXX :
	strait = DEF_LOC_SEC( CONFIG, selsec )

	print('			        >>>>   The concerned section is: '+ strait['name'])

	# Compute mean temp/sal & Ice thickness at V-point 
	##################################################
	# The calculation is done using data on the boundary (external part) and the first row (internal part)
	Tdata_read_V = 0.5 * ( Tdata_read.roll(shifts={'y': 1}) + Tdata_read ) * vmask
	Sdata_read_V = 0.5 * ( Sdata_read.roll(shifts={'y': 1}) + Sdata_read ) * vmask
	IHdata_read_V = 0.5 * ( IHdata_read.roll(shifts={'y': 1}) + IHdata_read ) * vmask.isel(z=0).squeeze()
	
	# Now select just the section to consider ( Back to Pyhton arrays indices starting at zero )
	jjloc = strait['jloc']-1
	iis = strait['is']-1
	iie = strait['ie']-1
	print('							jloc   :'+ str(strait['jloc']))
	print('							istart :'+ str(strait['is']))
	print('							iend   :'+ str(strait['ie']))
	
	if strait['name'] == 'Bering' :
	    vmask[:,jjloc,iie::] = 0   # The vmask is changed a little bit to remove a small open ocean area eastward to Bering
	
	Tsec = Tdata_read_V.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
	Ssec = Sdata_read_V.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
	HIsec = IHdata_read_V.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
	Vsec = Vdata_read.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
	VIsec = IVdata_read.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
	fasec = ( e3v * e1v ).isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
	e1vsec = e1v.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
	vmsksec = vmask.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()

	print()
	print('							Data selection done' )
	print()

	# Call the function dedicated to the calculation
	CAL_VOLHEATSALTICE( data_dir, CONFIG, CASE, s_year, strait, Tsec, Ssec, Vsec, HIsec, VIsec, fasec, e1vsec, vmsksec ) 

	#------------------------------------------------------------------------------------------------------------------------
	########################################
	# Compute mean temperature and velocity along a short section close to FRAM strait
	########################################
	#------------------------------------------------------------------------------------------------------------------------
	FramObs = {'name':"FramObs"}
	
	if strait['name'] == 'FramS' : 
		strait = DEF_LOC_SEC( CONFIG, FramObs )
		
		# Back to Pyhton arrays indices starting at zero
		jjloc = strait['jloc']-1
		iis = strait['is']-1
		iie = strait['ie']-1
		print(' 				>>>>   Special treatment for section : '+ strait['name'])
		print() 
		
		# Compute vertical surfaces along the section
		e1te3t = e1t * e3t * tmask
		e1ve3v = e1v * e3v * vmask
		
		# Set to zero values below 700m & out of the range [ 5degE- 8.4degE] along Fram Strait @ 79DegN
		e1te3t = xr.where( gdept_0 > 700. , 0., e1te3t )
		e1ve3v = xr.where( gdept_0 > 700. , 0., e1ve3v )
		
		# Now select just the section to consider ( Back to Pyhton arrays indices starting at zero )
		Tsec = Tdata_read.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
		Vsec = Vdata_read.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
		e1te3tsec = e1te3t.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze() 
		e1ve3vsec = e1ve3v.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
		
		mean_T_FraOb = ( e1te3tsec * Tsec ).sum(dim=['z','x']) / e1te3tsec.sum(dim=['z','x']) 
		mean_V_FraOb = ( e1ve3vsec * Vsec ).sum(dim=['z','x']) / e1ve3vsec.sum(dim=['z','x'])
		
		# Save fields to be able to reload them later
		npy.savez(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_FramObs_STRAITSTrans_y'+str(s_year), mean_T_FraOb=mean_T_FraOb.values, mean_V_FraOb=mean_V_FraOb.values)

		# Recover the Fram strait for following steps 
		strait = DEF_LOC_SEC( CONFIG, {'name':"FramS"} )

	#------------------------------------------------------------------------------------------------------------------------
	########################################
	# Plot LONG TIME-SERIES
	########################################
	#------------------------------------------------------------------------------------------------------------------------

	if lgTS_ye-lgTS_ys+1 > 1 :
		print()
		print('				##################################################################  ')
		print('				##################################################################  ')
		print('				###### PLOT CLASSICAL DIAGS LONG TIME-SERIES THROUGH SECTIONS ####  ')
		print('				##################################################################  ')
		print('				##################################################################  ')
		print()
	
		PLOT_TRANS_TISE( data_dir, CONFIG, CASE, lgTS_ys, lgTS_ye, strait, NCDF_OUT ) 

#------------------------------------------------------------------------------------------------------------------------
########################################
# Plot SECTIONS
########################################
#------------------------------------------------------------------------------------------------------------------------

for selsec in XXSECPLOTXX :
	strait = DEF_LOC_SEC( CONFIG, selsec )

	print()
	print('				##################################################################  ')
	print('				##################################################################  ')
	print('				##################  SECTIONS PLOTS  ##############################  ')
	print('				##################################################################  ')
	print('				##################################################################  ')
	print()

	# Back to Pyhton arrays indices starting at zero
	jjloc = strait['jloc']-1
	iis = strait['is']-1
	iie = strait['ie']-1
	print('				>>>> The strait plotted is :'+ strait['name'])
	print('				      		jloc   :'+ str(strait['jloc']))
	print('				      		istart :'+ str(strait['is']))
	print('				      		iend   :'+ str(strait['ie']))

	# Need to select either column or lines depending the section
	if strait['name'] == 'Beaufort' : 
		Tsec = Tdata_read_mean.isel(y=slice(iis,iie+1),x=jjloc).compute().squeeze()
		Ssec = Sdata_read_mean.isel(y=slice(iis,iie+1),x=jjloc).compute().squeeze()
		Usec = Udata_read_mean.isel(y=slice(iis,iie+1),x=jjloc).compute().squeeze()
		Ksec = Kdata_read_mean.isel(y=slice(iis,iie+1),x=jjloc).compute().squeeze()
		Tsec_init = My_varTinit_mean.isel(y=slice(iis,iie+1),x=jjloc).compute().squeeze()
		Ssec_init = My_varSinit_mean.isel(y=slice(iis,iie+1),x=jjloc).compute().squeeze()

		z2D = gdept_0.isel(y=slice(iis,iie+1),x=jjloc).compute().values.squeeze()
		lon2D = npy.tile(lon.isel(y=slice(iis,iie+1),x=jjloc).values,(z2D.shape[0],1))
	else :
		Tsec = Tdata_read_mean.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
		Ssec = Sdata_read_mean.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
		Usec = Vdata_read_mean.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
		Ksec = Kdata_read_mean.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
		Tsec_init = My_varTinit_mean.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()
		Ssec_init = My_varSinit_mean.isel(x=slice(iis,iie+1),y=jjloc).compute().squeeze()

		if strait['name'] == 'Beaufort': 
			z2D = gdept_0.isel(y=slice(iis,iie+1),x=jjloc).compute().values.squeeze()
			lon2D = npy.tile(lon.isel(y=slice(iis,iie+1),x=jjloc).values,(z2D.shape[0],1))
		elif strait['name'] == 'Kara' : 
			z2D = gdept_0.isel(x=slice(iis,iie+1),y=jjloc).compute().values.squeeze()
			lon2D = npy.tile(lat.isel(x=slice(iis,iie+1),y=jjloc).values,(z2D.shape[0],1))
		else :
			z2D = gdept_0.isel(x=slice(iis,iie+1),y=jjloc).compute().values.squeeze()
			lon2D = npy.tile(lon.isel(x=slice(iis,iie+1),y=jjloc).values,(z2D.shape[0],1))

	PLOT_SECTION( CONFIG, CASE, strait, Tsec, Ssec, Usec, Tsec_init, Ssec_init, Ksec, z2D, lon2D, climyear, NCDF_OUT )
