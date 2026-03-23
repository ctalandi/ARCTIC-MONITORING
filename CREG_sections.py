#!/usr/bin/env python

import sys 
import matplotlib
matplotlib.use('Agg')
import numpy as npy
import matplotlib.pylab as plt
from netCDF4 import Dataset
import matplotlib as mpl
import scipy.io as sio
from checkfile import *
from CREG_sections_func import *
import subprocess
from datetime import datetime
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
	infield = var_crty
	ds_zgr = xr.open_dataset(locpath+locfile)[[infield['ze3'],'gdept_0']]
	ze3 = ds_zgr[infield['ze3']].squeeze()
	ze3 = ze3.rename({'nav_lev':'z'})
	gdept_0 = ds_zgr['gdept_0'].squeeze()
	gdept_0 = gdept_0.rename({'nav_lev':'z'})

locpath=grid_dir
locfile=CONFCASE+'_mesh_hgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	ds_mes = xr.open_dataset(locpath+locfile)[['e1t','e2t','e1v']]
	ze1t = ds_mes['e1t'].squeeze()
	ze2t = ds_mes['e2t'].squeeze()
	ze1v = ds_mes['e1v'].squeeze()

e1te2t = ze1t * ze2t * tmask[0,:,:]

time_dim = (e_year-s_year+1)*12

vmask3Dtime=npy.tile(vmask[:,:,:],(time_dim,1,1,1))

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
   ITdata_read = ds_Idata[infieldIH['name']].squeeze()
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

zstrait='' ; icount_str=0
for selsec in XXDIAGSSECXX :
	strait = DEF_LOC_SEC( CONFIG, selsec )

	print('			        >>>>   The concerned section is: '+ strait['name'])
	print()
	# Call the function dedicated to the calculation
#	CAL_VOLHEATSALTICE( data_dir, CONFIG, CASE, s_year, strait, Tdata_read, Sdata_read, Vdata_read, ITdata_read, IVdata_read, ze3, ze1t, ze1v, vmask3Dtime, time_dim ) 

	###########################################################################################################################
	###########################################################################################################################
	############################################ PLOT THE LONG TIME-SERIES ####################################################
	###########################################################################################################################
	###########################################################################################################################

	if lgTS_ye-lgTS_ys+1 > 1 :
		print()
		print('				##################################################################  ')
		print('				##################################################################  ')
		print('				###### PLOT CLASSICAL DIAGS LONG TIME-SERIES THROUGH SECTIONS ####  ')
		print('				##################################################################  ')
		print('				##################################################################  ')
		print()
	

#		PLOT_TRANS_TISE( data_dir, CONFIG, CASE, lgTS_ys, lgTS_ye, strait, NCDF_OUT ) 

#		# Plot the long Time-series if requested
#		########################################
#		lgtstime_dim=(lgTS_ye-lgTS_ys+1)*12
#
#		LongTS_volu_trans= []	  ; LongTS_heat_trans=[]     ; LongTS_salt_trans=[]  ; LongTS_icet_trans=[]
#		LongTS_voluIN_trans= []   ; LongTS_heatIN_trans=[]   ; LongTS_saltIN_trans=[]
#		LongTS_FramObs_T=[]   ;  LongTS_FramObs_V=[]
#		# Start to read all yearly files
#		################################
#		lgts_year=lgTS_ys    ;	  t_months=(npy.arange(12)*30.+15.)/365.   ;   start = 1
#		while  lgts_year <= lgTS_ye  :
#			print('   >>>>	 Read year:'+ str(lgts_year))
#
#			# Read FWC & Sea-ice extent
#			locpath=data_dir+'/DATA/'
#			locfile=CONFIG+'-'+CASE+'_'+strait['name']+'_STRAITSTrans_y'+str(lgts_year)+'.npz'
#			if chkfile(locpath+locfile) :
#				open_npzfile = npy.load(locpath+locfile,mmap_mode='r')
#				LongTS_volu_trans=npy.append(LongTS_volu_trans,open_npzfile['net_volu_trans'])
#				LongTS_heat_trans=npy.append(LongTS_heat_trans,open_npzfile['net_heat_trans'])
#				LongTS_salt_trans=npy.append(LongTS_salt_trans,open_npzfile['net_salt_trans'])
#				LongTS_icet_trans=npy.append(LongTS_icet_trans,open_npzfile['net_icet_trans'])
#				LongTS_voluIN_trans=npy.append(LongTS_voluIN_trans,open_npzfile['net_voluIN_trans'])
#				LongTS_heatIN_trans=npy.append(LongTS_heatIN_trans,open_npzfile['net_heatIN_trans'])
#				LongTS_saltIN_trans=npy.append(LongTS_saltIN_trans,open_npzfile['net_saltIN_trans'])
#			else:
#				LongTS_volu_trans=npy.append(LongTS_volu_trans,npy.arange(12)+npy.nan)
#				LongTS_heat_trans=npy.append(LongTS_heat_trans,npy.arange(12)+npy.nan)
#				LongTS_salt_trans=npy.append(LongTS_salt_trans,npy.arange(12)+npy.nan)
#				LongTS_icet_trans=npy.append(LongTS_icet_trans,npy.arange(12)+npy.nan)
#				LongTS_voluIN_trans=npy.append(LongTS_voluIN_trans,npy.arange(12)+npy.nan)
#				LongTS_heatIN_trans=npy.append(LongTS_heatIN_trans,npy.arange(12)+npy.nan)
#				LongTS_saltIN_trans=npy.append(LongTS_saltIN_trans,npy.arange(12)+npy.nan)
#
#			if strait['name'] == 'FramS' :
#				# Read mean T and V to compare to obs.
#				locpath=data_dir+'/DATA/'
#				locfile=CONFIG+'-'+CASE+'_FramObs_STRAITSTrans_y'+str(lgts_year)+'.npz'
#				if chkfile(locpath+locfile) :
#					open_npzfile = npy.load(locpath+locfile,mmap_mode='r')
#					LongTS_FramObs_T=npy.append(LongTS_FramObs_T,open_npzfile['mean_T_FraOb'])
#					LongTS_FramObs_V=npy.append(LongTS_FramObs_V,open_npzfile['mean_V_FraOb'])
#				else:
#					LongTS_FramObs_T=npy.append(LongTS_FramObs_T,npy.arange(12)+npy.nan)
#					LongTS_FramObs_V=npy.append(LongTS_FramObs_V,npy.arange(12)+npy.nan)
#
#			# Set the time axis
#			y_years=npy.tile(lgts_year,12)+t_months
#			if start == 1:
#				time_axis=y_years
#				start=0
#			else:
#				time_axis=npy.append(time_axis,y_years)
#
#			lgts_year+=1
#
#		# Plot the time-series over SEVERAL YEARS
#		########################################################
#
#		time_grid=npy.arange(lgTS_ys,lgTS_ye+2,1.,dtype=int)
#		newlocsx  = npy.array(time_grid,'f')
#		newlabelsx = npy.array(time_grid,'i')
#		
#		lgtsclimyear=str(lgTS_ys)+str(lgTS_ye)
#
#		if strait['name'] == 'SouthG': xwind = 310  
#		if strait['name'] == 'FramS' or strait['name'] == 'Davis' or strait['name'] == 'Bering' : xwind = 410	 
#		
#		plt.clf()
#
#		# Plot the net volume transport
#		###############################
#		ax=plt.subplot(xwind+1)
#		plt.plot(time_axis, LongTS_volu_trans*1e-6   , 'r', linewidth=0.7)
#		if strait['name'] == 'FramS' :
#			plt.plot(time_axis, LongTS_voluIN_trans*1e-6   , 'g', linewidth=0.7)
#			plt.plot(time_axis, (LongTS_volu_trans-LongTS_voluIN_trans)*1e-6   , 'k', linewidth=0.7)
#			plt.text(lgTS_ye+1.,-2. ,str(npy.round(npy.nanmean(LongTS_volu_trans*1e-6),decimals=2)),color='r',size=8)
#			plt.text(lgTS_ye+1.,3.	,str(npy.round(npy.nanmean(LongTS_voluIN_trans*1e-6),decimals=2)),color='g',size=8)
#			plt.text(lgTS_ye+1.,-6. ,str(npy.round(npy.nanmean((LongTS_volu_trans-LongTS_voluIN_trans)*1e-6),decimals=2)),color='k',size=8)
#			#plt.text(lgTS_ys,7. ,'Obs: -1.6$\pm$3.9',color='m',size=8)
#			plt.ylim([-10.,10.])
#		elif strait['name'] == 'Davis' :
#			plt.plot(time_axis, LongTS_voluIN_trans*1e-6   , 'g', linewidth=0.7)
#			plt.plot(time_axis, (LongTS_volu_trans-LongTS_voluIN_trans)*1e-6   , 'k', linewidth=0.7)
#			plt.text(lgTS_ye+1.,-0.5 ,str(npy.round(npy.nanmean(LongTS_volu_trans*1e-6),decimals=2)),color='r',size=8)
#			plt.text(lgTS_ye+1.,2.	 ,str(npy.round(npy.nanmean(LongTS_voluIN_trans*1e-6),decimals=2)),color='g',size=8)
#			plt.text(lgTS_ye+1.,-2. ,str(npy.round(npy.nanmean((LongTS_volu_trans-LongTS_voluIN_trans)*1e-6),decimals=2)),color='k',size=8)
#			#plt.text(lgTS_ys,2. ,'Obs: -3.1$\pm$0.7',color='m',size=8)
#			plt.ylim([-4.,4.])
#		elif strait['name'] == 'Bering' :
#			plt.text(lgTS_ye+1.,0.75 ,str(npy.round(npy.nanmean(LongTS_volu_trans*1e-6),decimals=2)),color='r',size=8)
#			#plt.text(lgTS_ys,0.5 ,r'Obs: 1.$\pm$0.2',color='m',size=8)
#			plt.ylim([0.,2.])
#		elif strait['name'] == 'SouthG' :
#			plt.text(lgTS_ye+1.,-2 ,str(npy.round(npy.nanmean(LongTS_volu_trans*1e-6),decimals=2)),color='r',size=8)
#			plt.ylim([-4.,2.])
#		
#		plt.title(CASE+' integrated quantities @ '+strait['name']+' \n '+str(lgtsclimyear),size=9)
#		plt.xticks(newlocsx,newlabelsx,size=6)
#		plt.ylabel('Volume flux \n'+'(Sv)',size=7)
#		plt.yticks(size=6)
#		plt.setp(ax.get_xticklabels(),visible=False)
#		plt.grid(True,linestyle='--',color='grey',alpha=0.7)
#		
#		# Plot the heat transport
#		#########################
#		if strait['name'] == 'SouthG' :
#		   z_alpha=1e-15
#		else:
#		   z_alpha=1e-12
#		ax=plt.subplot(xwind+2)
#		curve_tot = plt.plot(time_axis, LongTS_heat_trans*z_alpha   , 'r', linewidth=0.7)
#		if strait['name'] == 'FramS' :
#			curve_in = plt.plot(time_axis, LongTS_heatIN_trans*z_alpha   , 'g', linewidth=0.7)
#			curve_ou = plt.plot(time_axis, (LongTS_heat_trans-LongTS_heatIN_trans)*z_alpha	 , 'k', linewidth=0.7)
#			plt.text(lgTS_ye+1.,30	,str(npy.round(npy.nanmean(LongTS_heat_trans*z_alpha),decimals=2)),color='r',size=8)
#			plt.text(lgTS_ye+1.,40	,str(npy.round(npy.nanmean(LongTS_heatIN_trans*z_alpha),decimals=2)),color='g',size=8)
#			plt.text(lgTS_ye+1.,10	,str(npy.round(npy.nanmean((LongTS_heat_trans-LongTS_heatIN_trans)*z_alpha),decimals=2)),color='k',size=8)
#			#plt.text(lgTS_ys,40 ,r'Obs: 62.$\pm$17',color='m',size=8)
#			plt.ylim([-10,50])
#		elif strait['name'] == 'Davis' :
#			curve_in = plt.plot(time_axis, LongTS_heatIN_trans*z_alpha   , 'g', linewidth=0.7)
#			curve_ou = plt.plot(time_axis, (LongTS_heat_trans-LongTS_heatIN_trans)*z_alpha	 , 'k', linewidth=0.7)
#			plt.text(lgTS_ye+1.,0 ,str(npy.round(npy.nanmean(LongTS_heat_trans*z_alpha),decimals=2)),color='r',size=8)
#			plt.text(lgTS_ye+1.,20	,str(npy.round(npy.nanmean(LongTS_heatIN_trans*z_alpha),decimals=2)),color='g',size=8)
#			plt.text(lgTS_ye+1.,10	,str(npy.round(npy.nanmean((LongTS_heat_trans-LongTS_heatIN_trans)*z_alpha),decimals=2)),color='k',size=8)
#			#plt.text(lgTS_ys,23 ,r'Obs: 28.$\pm$3',color='m',size=8)
#			plt.ylim([-20,30])
#		elif strait['name'] == 'Bering' :
#			plt.text(lgTS_ye+1.,30 ,str(npy.round(npy.nanmean(LongTS_heat_trans*z_alpha),decimals=2)),color='r',size=8)
#			#plt.text(lgTS_ys,40 ,r'Obs: 13.$\pm$2',color='m',size=8)
#			plt.ylim([-10,50.])
#		elif strait['name'] == 'SouthG' :
#			plt.text(lgTS_ye+1.,0.6 ,str(npy.round(npy.nanmean(LongTS_heat_trans*z_alpha),decimals=2)),color='r',size=8)
#			plt.ylim([0,1.])
#		
#		plt.xticks(newlocsx,newlabelsx,size=6)
#		plt.ylabel('Heat flux \n'+'(TW)',size=7)
#		plt.setp(ax.get_xticklabels(),rotation=45)
#		plt.yticks(size=6)
#		if strait['name'] == 'Bering' or strait['name'] == 'FramS' or strait['name'] == 'Davis' or strait['name'] == 'SouthG' : plt.setp(ax.get_xticklabels(),visible=False)
#		plt.grid(True,linestyle='--',color='grey',alpha=0.7)
#		
#		if strait['name'] == 'Bering' or strait['name'] == 'FramS' or strait['name'] == 'Davis' or strait['name'] == 'SouthG':
#			# Plot the salt transport
#			##########################
#			ax=plt.subplot(xwind+3)
#			plt.plot(time_axis, LongTS_salt_trans*1e-3   , 'r', label='Net', linewidth=0.7)
#			if strait['name'] == 'FramS' :
#				plt.plot(time_axis, LongTS_saltIN_trans*1e-3   , 'g', label='Northward', linewidth=0.7)
#				plt.plot(time_axis, (LongTS_salt_trans-LongTS_saltIN_trans)*1e-3   , 'k', label='Southward', linewidth=0.7)
#				plt.ylim([-160.,50.])
#				plt.text(lgTS_ye+1.,-80 ,str(npy.round(npy.nanmean(LongTS_salt_trans*1e-3),decimals=2)),color='r',size=8)
#				plt.text(lgTS_ye+1.,-10 ,str(npy.round(npy.nanmean(LongTS_saltIN_trans*1e-3),decimals=2)),color='g',size=8)
#				plt.text(lgTS_ye+1.,-120,str(npy.round(npy.nanmean((LongTS_salt_trans-LongTS_saltIN_trans)*1e-3),decimals=2)),color='k',size=8)
#				#plt.text(lgTS_ys,10 ,r'Obs: 70.$\pm$37',color='m',size=8)
#			elif strait['name'] == 'Davis' :
#				plt.plot(time_axis, LongTS_saltIN_trans*1e-3   , 'g', label='Northward', linewidth=0.7)
#				plt.plot(time_axis, (LongTS_salt_trans-LongTS_saltIN_trans)*1e-3   , 'k', label='Southward', linewidth=0.7)
#				plt.ylim([-200.,80.])
#				plt.text(lgTS_ye+1.,-60 ,str(npy.round(npy.nanmean(LongTS_salt_trans*1e-3),decimals=2)),color='r',size=8)
#				plt.text(lgTS_ye+1., 0 ,str(npy.round(npy.nanmean(LongTS_saltIN_trans*1e-3),decimals=2)),color='g',size=8)
#				plt.text(lgTS_ye+1.,-100,str(npy.round(npy.nanmean((LongTS_salt_trans-LongTS_saltIN_trans)*1e-3),decimals=2)),color='k',size=8)
#				#plt.text(lgTS_ys,25 ,r'Obs: 119.$\pm$14',color='m',size=8)
#			elif strait['name'] == 'Bering' :
#				plt.ylim([0.,140.])
#				plt.text(lgTS_ye+1.,50 ,str(npy.round(npy.nanmean(LongTS_salt_trans*1e-3),decimals=2)),color='r',size=8)
#				#plt.text(lgTS_ys,120 ,r'Obs: -72.$\pm$14',color='m',size=8)
#			elif strait['name'] == 'SouthG' :
#				plt.ylim([-200.,50.])
#				plt.text(lgTS_ye+1.,-100 ,str(npy.round(npy.nanmean(LongTS_salt_trans*1e-3),decimals=2)),color='r',size=8)
#			
#			plt.ylabel('Freshwater \n'+'(mSv)',size=7)
#			plt.xticks(newlocsx,newlabelsx,size=5)
#			#plt.setp(ax.get_xticklabels(),visible=False)
#			plt.setp(ax.get_xticklabels(),rotation=90)
#			plt.yticks(size=6)
#			if strait['name'] == 'FramS' or strait['name'] == 'Davis' : plt.setp(ax.get_xticklabels(),visible=False)
#			plt.grid(True,linestyle='--',color='grey',alpha=0.7)
#
#			plt.legend(loc='lower left',ncol=3)
#			leg = plt.gca().get_legend()
#			ltext = leg.get_texts()
#			plt.setp(ltext, fontsize=5.)
#
#
#		if strait['name'] == 'FramS' or strait['name'] == 'Davis' or strait['name'] == 'Bering' :
#			# Plot the Ice export
#			#####################
#			ax=plt.subplot(xwind+4)
#			z_ialpha=1.e-4
#			plt.plot(time_axis, LongTS_icet_trans*z_ialpha , 'r', linewidth=0.7)
#			if strait['name'] == 'FramS' :
#			   plt.text(lgTS_ye+1.,-3.  ,str(npy.round(npy.nanmean(LongTS_icet_trans*z_ialpha),decimals=2)),color='r',size=8)
#			   #plt.text(lgTS_ys,-10. ,r'Obs: -5.$\pm$2',color='m',size=8)
#			   plt.ylim([-25.,1])
#			elif strait['name'] == 'Davis':
#			   plt.text(lgTS_ye+1.,-1.5  ,str(npy.round(npy.nanmean(LongTS_icet_trans*z_ialpha),decimals=2)),color='r',size=8)
#			   plt.ylim([-7.,1])
#			elif strait['name'] == 'Bering':
#			   plt.text(lgTS_ye+1.,1.5  ,str(npy.round(npy.nanmean(LongTS_icet_trans*z_ialpha),decimals=2)),color='r',size=8)
#			   plt.ylim([-1.,4.])
#			plt.ylabel('Ice export \n '+r'(x10$^{-2}$Sv)',size=7)
#			plt.xticks(newlocsx,newlabelsx,size=5)
#			plt.setp(ax.get_xticklabels(),rotation=90)
#			plt.yticks(size=6)
#			plt.grid(True,linestyle='--',color='grey',alpha=0.7)
#
#		plt.savefig(CONFIG+'-'+CASE+'_STRAITS_'+strait['name']+'_NetVoluHeatSalt_LGTS_y'+str(lgTS_ys)+'LASTy.pdf')
#
#		if NCDF_OUT:
#			# FWC field 
#			#######################
#			cmd_ddate="date"  ;  get_output = subprocess.check_output(cmd_ddate)
#			nc_f = './NETCDF/'+CONFIG+'-'+CASE+'_STRAITS_'+strait['name']+'_NetVoluHeatSalt_LGTS_y'+str(lgTS_ys)+'LASTy.nc'
#			w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
#			w_nc_fid.description = "Diagnostics have been calculated using the Arctic monitoring tool "
#			w_nc_fid.date=get_output.decode("utf-8")
#			w_nc_fid.createDimension('time', time_axis.shape[0])
#			#w_nc_fid.createDimension('time_counter', None)
#
#			w_nc_var = w_nc_fid.createVariable('LongTS_volu_trans', 'f4', ('time'))
#			w_nc_var.long_name='Model net volume flux '
#			w_nc_var.units="m3/s"
#			w_nc_fid.variables['LongTS_volu_trans'][:] = LongTS_volu_trans
#
#			w_nc_var = w_nc_fid.createVariable('LongTS_voluIN_trans', 'f4', ('time'))
#			w_nc_var.long_name='Model northward volume flux '
#			w_nc_var.units="m3/s"
#			w_nc_fid.variables['LongTS_voluIN_trans'][:] = LongTS_voluIN_trans
#
#			w_nc_var = w_nc_fid.createVariable('LongTS_voluOU_trans', 'f4', ('time'))
#			w_nc_var.long_name='Model southward volume flux '
#			w_nc_var.units="m3/s"
#			w_nc_fid.variables['LongTS_voluOU_trans'][:] = LongTS_volu_trans-LongTS_voluIN_trans
#
#			w_nc_var = w_nc_fid.createVariable('LongTS_heat_trans', 'f4', ('time'))
#			w_nc_var.long_name='Model net heat flux '
#			w_nc_var.units="W"
#			w_nc_fid.variables['LongTS_heat_trans'][:] = LongTS_heat_trans
#
#			w_nc_var = w_nc_fid.createVariable('LongTS_heatIN_trans', 'f4', ('time'))
#			w_nc_var.long_name='Model northward heat flux '
#			w_nc_var.units="W"
#			w_nc_fid.variables['LongTS_heatIN_trans'][:] = LongTS_heatIN_trans
#
#			w_nc_var = w_nc_fid.createVariable('LongTS_heatOU_trans', 'f4', ('time'))
#			w_nc_var.long_name='Model southward flux '
#			w_nc_var.units="W"
#			w_nc_fid.variables['LongTS_heatOU_trans'][:] = LongTS_heat_trans-LongTS_heatIN_trans
#
#			w_nc_var = w_nc_fid.createVariable('LongTS_salt_trans', 'f4', ('time'))
#			w_nc_var.long_name='Model net salt flux using a reference salinity of 34.8 PSU'
#			w_nc_var.units="Sv"
#			w_nc_fid.variables['LongTS_salt_trans'][:] = LongTS_salt_trans
#
#			w_nc_var = w_nc_fid.createVariable('LongTS_saltIN_trans', 'f4', ('time'))
#			w_nc_var.long_name='Model northward salt flux using a reference salinity of 34.8 PSU'
#			w_nc_var.units="Sv"
#			w_nc_fid.variables['LongTS_saltIN_trans'][:] = LongTS_saltIN_trans
#
#			w_nc_var = w_nc_fid.createVariable('LongTS_saltOU_trans', 'f4', ('time'))
#			w_nc_var.long_name='Model southward salt flux using a reference salinity of 34.8 PSU'
#			w_nc_var.units="Sv"
#			w_nc_fid.variables['LongTS_saltOU_trans'][:] = LongTS_salt_trans-LongTS_saltIN_trans
#
#			w_nc_var = w_nc_fid.createVariable('LongTS_icet_trans', 'f4', ('time'))
#			w_nc_var.long_name='Model southward ice export flux '
#			w_nc_var.units="Sv"
#			w_nc_fid.variables['LongTS_icet_trans'][:] = LongTS_icet_trans
#
#
#			w_nc_fid.close()  # close the file
#
#
#
#		if strait['name'] == 'FramS':
#			# Plot the mean Temperature through FRAM section
#			################################################
#			# Read observations 
#			locpath='./'
#			locfile='FRAM_inflow.mat'
#			if chkfile(locpath+locfile) :
#				FramObservations_read = sio.loadmat(locpath+locfile,squeeze_me=True)
#				select_FramObs_T = npy.array(FramObservations_read['T_5degE_700m_obs'])
#				select_FramObs_V = npy.array(FramObservations_read['V_5degE_700m_obs'])
#			else:
#				select_FramObs_T = npy.arange((2004-1997+1)*12)+npy.nan
#				select_FramObs_V = npy.arange((2004-1997+1)*12)+npy.nan
#
#			# Set the time axis for observations
#			cyear=1997  ;  start=1
#			while cyear <= 2004 :
#				y_years=npy.tile(cyear,12)+t_months
#				if start == 1:
#					time_axis_obs=y_years
#					start=0
#				else:
#					time_axis_obs=npy.append(time_axis_obs,y_years)
#				cyear+=1
#
#
#			xwind = 210
#			
#			plt.clf()
#			ax=plt.subplot(xwind+1)
#			plt.plot(time_axis    , LongTS_FramObs_T, 'k', linewidth=0.7)
#			plt.plot(time_axis_obs, select_FramObs_T, 'g', label='Obs', linewidth=0.7)
#			plt.text(lgTS_ye+1.,2.	,str(npy.round(npy.nanmean(LongTS_FramObs_T),decimals=2)),color='k',size=8)
#			plt.text(2005.,3.5	,str(npy.round(npy.nanmean(select_FramObs_T),decimals=2)),color='g',size=8)
#			plt.title(CASE+' Mean Temp & Velocity above 700m & 5degE @ FramObs \n '+str(lgtsclimyear),size=9)
#			plt.ylabel('Mean temperature \n'+r'($^{\circ}$C)', size=7)
#			plt.ylim([0.,4.])
#			plt.xticks(newlocsx,newlabelsx,size=5)
#			plt.setp(ax.get_xticklabels(),rotation=90)
#			plt.yticks(size=6)
#			plt.grid(True)
#			plt.legend(loc='upper left')
#			leg = plt.gca().get_legend()
#			ltext = leg.get_texts()
#			plt.setp(ltext, fontsize=5.)
#			
#			ax=plt.subplot(xwind+2)
#			plt.plot(time_axis    , LongTS_FramObs_V*100 , 'k', linewidth=0.7)
#			plt.plot(time_axis_obs, select_FramObs_V, 'g', label='Obs', linewidth=0.7)
#			plt.text(lgTS_ye+1,10.	 ,str(npy.round(npy.nanmean(LongTS_FramObs_V)*100.,decimals=1)),color='k',size=8)
#			plt.text(2005.,10.   ,str(npy.round(npy.nanmean(select_FramObs_V)     ,decimals=1)),color='g',size=8)
#			plt.ylabel('Mean velocity \n '+r'(cm $s^{-1}$)', size=7)
#			plt.ylim([-2.,12.])
#			plt.xticks(newlocsx,newlabelsx,size=5)
#			plt.setp(ax.get_xticklabels(),rotation=90)
#			plt.yticks(size=6)
#			plt.grid(True)
#
#			plt.savefig(CONFIG+'-'+CASE+'_STRAITS_FramObs_meanVT_LGTS_y'+str(lgTS_ys)+'LASTy.pdf')
#
#			if NCDF_OUT:
#				# FWC field 
#				#######################
#				cmd_ddate="date"  ;  get_output = subprocess.check_output(cmd_ddate)
#				nc_f = './NETCDF/'+CONFIG+'-'+CASE+'_STRAITS_FramObs_meanVT_LGTS_y'+str(lgTS_ys)+'LASTy.nc'
#				w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
#				w_nc_fid.description = "Diagnostics have been calculated using the Arctic monitoring tool "
#				w_nc_fid.data = "Observations time series have been retrieved from the following URL: http://www.whoi.edu/page.do?pid=30914"
#				w_nc_fid.date=get_output.decode("utf-8")
#				w_nc_fid.createDimension('time_mod', time_axis.shape[0])
#				w_nc_fid.createDimension('time_obs', time_axis_obs.shape[0])
#				#w_nc_fid.createDimension('time_counter', None)
#
#				w_nc_var = w_nc_fid.createVariable('LongTS_FramObs_T', 'f4', ('time_mod'))
#				w_nc_var.long_name='Model mean temperature above 700m & in the range [ 5degE- 8.4degE] along Fram Strait @ 79DegN'
#				w_nc_var.units="DegC"
#				w_nc_fid.variables['LongTS_FramObs_T'][:] = LongTS_FramObs_T
#
#				w_nc_var = w_nc_fid.createVariable('select_FramObs_T', 'f4', ('time_obs'))
#				w_nc_var.long_name='Observed monthly means temperature above 700m & in the range [ 5degE- 8.4degE] along Fram Strait @ 79DegN'
#				w_nc_var.units="DegC"
#				w_nc_fid.variables['select_FramObs_T'][:] = select_FramObs_T
#
#				w_nc_var = w_nc_fid.createVariable('LongTS_FramObs_V', 'f4', ('time_mod'))
#				w_nc_var.long_name='Model mean velocity above 700m & in the range [ 5degE- 8.4degE] along Fram Strait @ 79DegN'
#				w_nc_var.units="m/s"
#				w_nc_fid.variables['LongTS_FramObs_V'][:] = LongTS_FramObs_V
#
#				w_nc_var = w_nc_fid.createVariable('select_FramObs_V', 'f4', ('time_obs'))
#				w_nc_var.long_name='Observed monthly means velocities above 700m & in the range [ 5degE- 8.4degE] along Fram Strait @ 79DegN'
#				w_nc_var.units="cm/s"
#				w_nc_fid.variables['select_FramObs_V'][:] = select_FramObs_V
#
#				w_nc_fid.close()  # close the file

################################################################################################################
################################################################################################################
#########################################  SECTIONS PLOTS  #####################################################
################################################################################################################
################################################################################################################

for selsec in XXSECPLOTXX:
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

