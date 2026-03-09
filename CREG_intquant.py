#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import numpy as npy
import matplotlib.pylab as plt
from matplotlib import rcParams
import matplotlib as mpl
import sys 
from checkfile import *
from CREG_intquant_func import *
import xarray as xr
from datetime import datetime
from fsspec.implementations.local import LocalFileSystem
fs = LocalFileSystem()

#rcParams['text.usetex']=True
rcParams['font.family']='serif'

s_year=XXSYEAXX
e_year=XXEYEAXX
lgTS_ys=XXLGTSSXX
lgTS_ye=XXLGTSEXX
xiosfreq=XXXIOSFREQXX
NCDF_OUT=XXNCDFOUTXX
main_dir='./'
CONFIG='XXCONFXX'   ; CASE='XXCASEXX'	  ;  CASE2='XXCASE2XX'
CONFCASE=CONFIG+'-'+CASE   ; CONFCASE2=CONFIG+'-'+CASE2
data_dir=main_dir+CONFIG+'/'+CONFCASE+'-MEAN/'	  ; data_dir2=main_dir+CONFIG+'/'+CONFCASE2+'-MEAN/'
grid_dir=main_dir+CONFIG+'/GRID/'

DIR_FIG_OUT='./'

# Infos concernant les climatologies sur la periode de la simulation
climyear=str(s_year)+str(e_year)

print() 
print("				     Configuration :" + CONFCASE)
print("				     Period	   :" + str(s_year)+" - "+str(e_year))
print()


#####################################################################
var_temp={'name':"votemper",'units':u"degC",'_FillValue': 0.,'fext':"T",'igrd':1 ,'ze3':"e3t_0",'gdep':"gdept",'long_name':"Temperature"}
var_sali={'name':"vosaline",'units':u"psu" ,'_FillValue': 0.,'fext':"S",'igrd':1 ,'ze3':"e3t_0",'gdep':"gdept",'long_name':"Salinity"}
var_crty={'name':"vomecrty",'units':u"m/s" ,'_FillValue': 0.,'fext':"V",'igrd':3 ,'ze3':"e3v_0",'gdep':"gdepv",'long_name':"Meridional velocity"}
var_hthi={'name':"sivolu"  ,'units':u"m"   ,'_FillValue': 0.,'fext':"I",'igrd':1 ,'ze3':"e3t_0",'gdep':"gdept",'long_name':"Ice thickness "}
var_Wcur={'name':"socurl"  ,'units':u"m/s" ,'_FillValue': 0.,'fext':"Xi",'igrd':1 ,'ze3':"e3f_0",'gdep':"gdepw",'long_name':"Surface ocean stress curl "}
var_ISfx={'name':"sfxice"  ,'units':u"kg/m2/s" ,'_FillValue': 0.,'fext':"sfxice",'igrd':1 ,'ze3':"e3t_0",'gdep':"gdept",'long_name':"Surface ice salt flux "}
var_IVfx={'name':"vfxice"  ,'units':u"kg/m2/s" ,'_FillValue': 0.,'fext':"vfxice",'igrd':1 ,'ze3':"e3t_0",'gdep':"gdept",'long_name':"Surface ice mass flux"}
MassSaltFLX=True

#------------------------------------------------------------------------------------------------------------------------
########################################
# Read GRID 
########################################

locpath=grid_dir
locfile=CONFCASE+'_mask.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	ds_grd=xr.open_dataset(locpath+locfile)
	tmask = ds_grd['tmask'].isel(time_counter=0)
	fmask = ds_grd['fmask'].isel(time_counter=0)
	tmask = tmask.rename({'nav_lev':'z'})

locfile=CONFCASE+'_coordinates.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	ds_coo=xr.open_dataset(locpath+locfile)
	lon = ds_coo['glamt'].squeeze()
	lat = ds_coo['gphit'].squeeze()
	ffCor = ds_coo['ff_t'].squeeze()

# Read appropriate scale factors 
infield=var_sali
locpath=grid_dir
locfile=CONFCASE+'_mesh_zgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	infield=var_temp
	ds_zgr=xr.open_dataset(locpath+locfile)
	ze3 = ds_zgr[infield['ze3']].squeeze()
	ze3 = ze3.rename({'nav_lev':'z'})
	ze33D= ze3 * tmask 

locpath=grid_dir
locfile=CONFCASE+'_mesh_hgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	ds_mes=xr.open_dataset(locpath+locfile)
	ze1t = ds_mes['e1t'].squeeze()
	ze2t = ds_mes['e2t'].squeeze()
	ze1f = ds_mes['e1f'].squeeze()
	ze2f = ds_mes['e2f'].squeeze()

areaT = ze1t*ze2t*tmask[0,:,:]
areaF = ze1f*ze2f*fmask[0,:,:]

#------------------------------------------------------------------------------------------------------------------------
########################################
# Read input  data 
########################################
infield=var_sali
print("						   ")
print("					       >>>>>>>>>  The curent variable processed is :	    "+ infield['fext'])
print("						   ")
# List files to be read
locpath=data_dir+'/'+str(s_year)+'/'+xiosfreq+'/'
locfile=CONFIG+'-'+CASE+'_y'+str(s_year)+'m??.'+xiosfreq+'_gridT.nc'
SAL_files = [f for f in fs.glob(locpath+locfile)]

drp_var=["time_centered", "deptht_bounds","time_centered_bounds","time_counter_bounds"]
if len(SAL_files) == 12 :
   ds_Sdata = xr.open_mfdataset(locpath+locfile, engine="netcdf4", concat_dim=["time_counter"], combine='nested', parallel=True, drop_variables=drp_var)
   Sdata_read = ds_Sdata[infield['name']]
   # Fill NaN with zero
   Sdata_read = Sdata_read.fillna(0)
   Sdata_read = Sdata_read.rename({'deptht':'z'})
   # Mask land grid points
   Sdata_read = Sdata_read * tmask

print("						   ")
print("					       >>>>>>>>>  The curent variable processed is :	    at_i(:,:)"	)
print("						   ")
# List files to be read
locpath=data_dir+'/'+str(s_year)+'/'+xiosfreq+'/'
locfile=CONFIG+'-'+CASE+'_y'+str(s_year)+'m??.'+xiosfreq+'_icemod.nc'
ICE_files = [f for f in fs.glob(locpath+locfile)]

drp_var=["time_centered_bounds","time_counter_bounds","siages"]
if len(ICE_files) == 12 :
   ds_Idata = xr.open_mfdataset(locpath+locfile, engine="netcdf4", concat_dim=["time_counter"], combine='nested', parallel=True, drop_variables=drp_var)
   Ifrdata_read = ds_Idata['siconc']
   Ivedata_read = ds_Idata['sivelo']
   Ithdata_read = ds_Idata['sivolu']
   if MassSaltFLX: 
      infield1=var_ISfx ; infield2=var_IVfx
      print("						   ")
      print("					       >>>>>>>>>  The curent variables processed are :	    "+ infield1['fext']+' & '+infield2['fext'])
      print("						   ")
      SFXdata_read = ds_Idata[infield1['name']]
      VFXdata_read = ds_Idata[infield2['name']]

   # Fill NaN with zero
   Ifrdata_read= Ifrdata_read.fillna(0)
   Ivedata_read= Ivedata_read.fillna(0)
   Ithdata_read= Ithdata_read.fillna(0)
   if MassSaltFLX: 
      SFXdata_read = SFXdata_read.fillna(0)
      VFXdata_read = VFXdata_read.fillna(0)
      
      SFXdata_read = xr.where(npy.abs(SFXdata_read) > 2e10, 0., SFXdata_read)
      VFXdata_read = xr.where(npy.abs(VFXdata_read) > 1e10, 0., VFXdata_read)
	
   # Mask land grid points
   Ifrdata_read= Ifrdata_read * tmask[0,:,:]
   Ivedata_read= Ivedata_read * tmask[0,:,:]
   Ithdata_read= Ithdata_read * tmask[0,:,:]
   if MassSaltFLX: 
      SFXdata_read = SFXdata_read * tmask[0,:,:]
      VFXdata_read = VFXdata_read * tmask[0,:,:]

infield=var_Wcur
print("						   ")
print("					       >>>>>>>>>  The curent variable processed is :	    "+ infield['fext'])
print("						   ")
# List files to be read
locpath=data_dir+'/'+str(s_year)+'/'+xiosfreq+'/'
locfile=CONFIG+'-'+CASE+'_y'+str(s_year)+'.'+xiosfreq+'_SurOceCurl.nc'
CUR_files = [f for f in fs.glob(locpath+locfile)]

if chkfile(locpath+locfile) :
   ds_Cdata = xr.concat([xr.open_dataset(f) for f in CUR_files], dim='time_counter', data_vars=None )
   Scurldata_read = ds_Cdata[infield['name']]
   Scurldata_read = xr.where(npy.abs(Scurldata_read) > 1e10 , 0., Scurldata_read)
   # Fill NaN with zero
   Scurldata_read = Scurldata_read.fillna(0)
   # Mask land grid points
   Scurldata_read = Scurldata_read * fmask[0,:,:]


#------------------------------------------------------------------------------------------------------------------------
################################################################################################################
################################################################################################################
##################################  FWC, ICE VOLUME & CONCENTRATION etc ...Diags ###############################
################################################################################################################
################################################################################################################
print()
print("			FWC, ICE VOLUME, AREA, EXTENSION, DRIFT , WIND STRESS MODULE Diags")
print()

plt.clf()

# Absolute Salinity 
Sref=34.80*1.004715
# Practical Salinity UNIT 
#Sref=34.80 

########################################
# Diagnostics in the Beaufort Gyre area 
########################################

# Define a mask over this area on T-points
tmskBFG=npy.ones((tmask.shape[1],tmask.shape[2]))
tmskBFG[npy.where(lon[:,:] > -130.)]=0.
tmskBFG[npy.where(lon[:,:] < -170.)]=0.
tmskBFG[npy.where(lat[:,:] >  80.5)]=0.
tmskBFG[npy.where(lat[:,:] <  70.5)]=0.
areaTBG = areaT * tmskBFG

# Still the same area but now for the vorticity field or F-point
fmskBFG=npy.ones((fmask.shape[1],fmask.shape[2]))
fmskBFG[npy.where(lon[:,:] > -130.)]=0.
fmskBFG[npy.where(lon[:,:] < -170.)]=0.
fmskBFG[npy.where(lat[:,:] >  80.5)]=0.
fmskBFG[npy.where(lat[:,:] <  70.5)]=0.
areaFBG = areaF * fmskBFG

# Freshwater content in the Beaufort Gyre box 
######################################################################################
fwcmask = xr.where(Sdata_read > Sref, 0., ze33D)
FW4D = (Sref - Sdata_read) / Sref * fwcmask
# Sum over depth 
fwc2D = FW4D.sum(dim="z")
#fwc2D = fwc2D * tmask[0,:,:]
# Sum over the BG area 
FWC = npy.sum(fwc2D*areaTBG,axis=(1,2))

#print a specific point to debug within the Beaufort Gyre
dbg=False
if CONFIG == 'CREG025.L75' :
	idbg=185 ; jdbg=515 ; kdbg=50	# CREG025.L75
elif CONFIG == 'CREG12.L75' :
	idbg=562 ; jdbg=1558 ; kdbg=50	 # CREG12.L75

if dbg:
   print("ze33D[6,0:kdbg,jdbg,idbg]:"+str(ze33D[0:kdbg,jdbg,idbg])+"  "+str( Sdata_read[6,0:kdbg,jdbg,idbg]))
   print("fwc2D[6,jdbg,idbg]:"+str(fwc2D[6,jdbg,idbg]))


# Ice drift in the Beaufort Gyre box 
######################################################################################
Vice_mean = npy.sum(Ivedata_read*areaTBG,axis=(1,2))/npy.sum(areaTBG,axis=(0,1))

# Sea Ice extent over Arctic 
######################################################################################
sice_ext3D = xr.where(Ifrdata_read < 0.15, 0., areaT)
sice_ext = npy.sum(sice_ext3D, axis=(1,2))


rho=1027.5  # Same value as in Gianluca et al. JPO2017, Observations of Seasonal Upwelling and Downwelling in the Beaufort Sea Mediated by Sea Ice
WEkm_mean = npy.sum(Scurldata_read*areaFBG/(rho*ffCor),axis=(1,2))/npy.sum(areaFBG,axis=(0,1))

# Sea Ice volume and sal fluxes in the BG area
######################################################################################
if MassSaltFLX :
   SFX_mean = npy.sum(SFXdata_read*areaTBG,axis=(1,2))/npy.sum(areaTBG,axis=(0,1))
   VFX_mean = npy.sum(VFXdata_read*areaTBG,axis=(1,2))/npy.sum(areaTBG,axis=(0,1))

# Save diagnostics of the current year 
######################################################################################
npy.save(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_FWCTS_y'+str(s_year),npy.array(FWC))
npy.save(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_ICEExtTS_y'+str(s_year),npy.array(sice_ext))
npy.save(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_VICETS_y'+str(s_year),npy.array(Vice_mean))
npy.save(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_WEKMATS_y'+str(s_year),npy.array(WEkm_mean))
if MassSaltFLX :
   npy.save(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_ICESFXTS_y'+str(s_year),npy.array(SFX_mean))
   npy.save(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_ICEVFXTS_y'+str(s_year),npy.array(VFX_mean))

##########################################
# Diagnostics over the entire Arctic area
##########################################
# Define a proper Arctic box following Haine et al. GPC2015
# including all Arctic area + CAA + north to Davis, FRAM and Bering straits as east to BSO

tmskARCTIC=npy.ones((tmask.shape[1],tmask.shape[2]))

tmskARCTIC[npy.where(lat[:,:] <  62.)] = 0.
# Remove the Baffin bay area
tmskARCTIC[npy.where((lon[:,:] >= -90.) & (lon[:,:] < -45.) & (lat[:,:] <= 80.) )]=0.
# Adjust the Baffin bay area north to Davis Strait 66.25N
tmskARCTIC[npy.where((lon[:,:] >= -70.) & (lon[:,:] < -45.) & (lat[:,:] >= 62.) & (lat[:,:] < 66.25) )] = 0.
tmskARCTIC[npy.where((lon[:,:] >= -70.) & (lon[:,:] < -65.) & (lat[:,:] >= 66.25) & (lat[:,:] < 67.) )] = 0.
# Remove all GIN seas area
tmskARCTIC[npy.where((lon[:,:] < 17.) & (lon[:,:] > -45.) & (lat[:,:] < 79.))] = 0.
# Remove Hudson area
tmskARCTIC[npy.where((lon[:,:] < -70.) & (lon[:,:] > -82. ) & (lat[:,:] < 70.5))] = 0.
tmskARCTIC[npy.where((lon[:,:] <= -82.) & (lon[:,:] > -100.) & (lat[:,:] < 67.))] = 0.
# Adriatic area 
tmskARCTIC[npy.where((lon[:,:] < 30.) & (lon[:,:] > 17.) & (lat[:,:] < 68.))] = 0.

areaARCTIC = areaT * tmskARCTIC

BigFWC = npy.sum(fwc2D*areaARCTIC,axis=(1,2))

npy.save(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_BIGFWCTS_y'+str(s_year),npy.array(BigFWC))


#------------------------------------------------------------------------------------------------------------------------
################################################################################################################
################################################################################################################
#################################   Beaufort Gyre FWC, Ice thick & Drift Time Series  ##########################
################################################################################################################
################################################################################################################

if lgTS_ye-lgTS_ys+1 > 1 :
	print()
	print("				##################################################################  " )
	print("				##################################################################  " ) 
	print("				######### PLOT DIAGS LONG TIME-SERIES WITHIN BEAUFORT GYRE #######  " ) 
	print("				##################################################################  " ) 
	print("				##################################################################  " ) 
	print()

	lgtstime_dim=(lgTS_ye-lgTS_ys+1)*12

	LongTS_FWC= []	   ; LongTS_sice_ext=[]   ; LongTS_ibgvoltot=[]  ; LongTS_ibgarea=[]
	LongTS_VEICE= []   ; LongTS_HIICE=[]	  ; LongTS_BIGFWC= []  ; LongTS_WEkm=[] ;  LongTS_SFX= []   ; LongTS_VFX=[]	 

	# Start to read all yearly files
	#######################################################################################################
	lgts_year=lgTS_ys    ;	  t_months=(npy.arange(12)*30.+15.)/365.   ;   start = 1
	while  lgts_year <= lgTS_ye  :
		print ()
		print("				>>>>   Read year:"+str(lgts_year))

		# Read FWC & Sea-ice extent
		locpath=data_dir+'/DATA/'
		locfile=CONFIG+'-'+CASE+'_FWCTS_y'+str(lgts_year)+'.npy'
		LongTS_FWC = READ_MOD_LGTS_DATA(locpath,locfile,LongTS_FWC)

		locfile=CONFIG+'-'+CASE+'_BIGFWCTS_y'+str(lgts_year)+'.npy'
		LongTS_BIGFWC = READ_MOD_LGTS_DATA(locpath,locfile,LongTS_BIGFWC)

		locfile=CONFIG+'-'+CASE+'_ICEExtTS_y'+str(lgts_year)+'.npy'
		#locfile=CONFIG+'-'+CASE+'_ICEExtTS_y'+str(s_year)+'.npy'
		LongTS_sice_ext = READ_MOD_LGTS_DATA(locpath,locfile,LongTS_sice_ext)

		locfile=CONFIG+'-'+CASE+'_WEKMATS_y'+str(lgts_year)+'.npy'
		LongTS_WEkm = READ_MOD_LGTS_DATA(locpath,locfile,LongTS_WEkm)

		if MassSaltFLX :
			locfile=CONFIG+'-'+CASE+'_ICESFXTS_y'+str(lgts_year)+'.npy'
			LongTS_SFX = READ_MOD_LGTS_DATA(locpath,locfile,LongTS_SFX)

			locfile=CONFIG+'-'+CASE+'_ICEVFXTS_y'+str(lgts_year)+'.npy'
			LongTS_VFX = READ_MOD_LGTS_DATA(locpath,locfile,LongTS_VFX)

		# Read Sea-ice volume & area
		locpath=data_dir+'/'+str(lgts_year)+'/1m/'
		locfile=CONFIG+'-'+CASE+'_y'+str(lgts_year)+'.1m_SBC_scalar.nc'
		if chkfile(locpath+locfile) :
			ds_sbc=xr.open_dataset(locpath+locfile)
			year_ibgvoltot = ds_sbc['ibgvoltot']
			year_ibgarea = ds_sbc['ibgarea']
			LongTS_ibgvoltot=npy.append(LongTS_ibgvoltot,year_ibgvoltot)
			LongTS_ibgarea=npy.append(LongTS_ibgarea,year_ibgarea)
		else:
			ibgvoltot, ibgarea = CAL_ICE_VOL_AREA(CONFIG,CASE,lgts_year,data_dir,xiosfreq, dom_area=areaT,tmask2D=tmask[0,:,:].squeeze())
			LongTS_ibgvoltot=npy.append(LongTS_ibgvoltot,ibgvoltot*1e-9)
			LongTS_ibgarea=npy.append(LongTS_ibgarea,ibgarea*1.e-6)

		# Read Mean drift & ice thickness and heat flux
		locpath=data_dir+'/DATA/'
		locfile=CONFIG+'-'+CASE+'_VICETS_y'+str(lgts_year)+'.npy'
		LongTS_VEICE = READ_MOD_LGTS_DATA(locpath,locfile,LongTS_VEICE)


		# Set the time axis
		y_years=npy.tile(lgts_year,12)+t_months
		if start == 1:
			time_axis=y_years
			start=0
		else:
			time_axis=npy.append(time_axis,y_years)

		lgts_year+=1

	# Read observations DATA from PIOMAS, NSIDC and IABP
	LongTS_OBS_icevol, LongTS_OBS_iceext, LongTS_OBS_iceare, IABPObservations, time_axis_obs, time_axis_PIO, time_axis_NSIDC = READ_OBS_LGTS_DATA(CONFIG,lgTS_ys,lgTS_ye)
	# Read observations DATA from Proshutinsky et al. GRL2018
	LongTS_OBS_CRFFWC, time_axis_CRFFWC = READ_OBS_LGTS_CRFFWC(lgTS_ys,lgTS_ye)
	# Read observations DATA from Gianluca et al. JPO2017
	LongTS_OBS_CRFEkm, time_axis_CRFEkm = READ_OBS_LGTS_CRFEkm(lgTS_ys,lgTS_ye)

	# Set time axis properly 
	#######################################################################################################
	time_grid=npy.arange(lgTS_ys,2018.,1.,dtype=int)
	newlocsx  = npy.array(time_grid,'f')
	newlabelsx = npy.array(time_grid,'i')

	# IABPObservations time axis
	newlocsx_obs  = npy.array(time_axis_obs,'f')
	newlabelsx_obs = npy.array(time_axis_obs,'i')
	full_time_axis = npy.arange(lgTS_ys,2017.,1.,dtype=int)
	full_newlocsx=npy.array(full_time_axis,'f')
	full_newlabelsx=npy.array(full_time_axis,'i')

	# PIOMAS time axis
	newlocsx_piomas  = npy.array(time_axis_PIO,'f')
	newlabelsx_piomas = npy.array(time_axis_PIO,'i')
	piomas_newlocsx=npy.array(npy.append(newlocsx,npy.unique(newlabelsx_piomas)),'f')
	piomas_newlabelsx=npy.array(npy.append(newlabelsx,npy.unique(newlabelsx_piomas)),'i')
	
	# NSIDC time axis
	newlocsx_nsidc	= npy.array(time_axis_NSIDC,'f')
	newlabelsx_nsidc = npy.array(time_axis_NSIDC,'i')
	nsidc_newlocsx=npy.array(npy.append(newlocsx,npy.unique(newlabelsx_nsidc)),'f')
	nsidc_newlabelsx=npy.array(npy.append(newlabelsx,npy.unique(newlabelsx_nsidc)),'i')

	# Ekman pumping time axis
	time_axis_CRFEkm=[]
	for pyear in npy.arange(2003,2015):
		y_years=npy.tile(pyear,12)+t_months
		time_axis_CRFEkm=npy.append(time_axis_CRFEkm,y_years)

		pyear+=1
	newlocsx_ekman	= npy.array(time_axis_CRFEkm,'f')
	newlabelsx_ekman = npy.array(time_axis_CRFEkm,'i')
	ekman_newlocsx=npy.array(npy.append(newlocsx,npy.unique(newlabelsx_ekman)),'f')
	ekman_newlabelsx=npy.array(npy.append(newlabelsx,npy.unique(newlabelsx_ekman)),'i')
	
	if dbg: print("new obs time"+str(newlocsx_ekman) )
	
	lgtsclimyear=str(lgTS_ys)+str(lgTS_ye)

	# Make plots 
	#######################################################################################################
	plt.clf()
	xwind=410
	
	# Ice volume 
	################
	ax=plt.subplot(xwind+1)
	plt.title(CASE+' Arctic Ice quantities over \n '+str(lgtsclimyear),size=9)
	plt.plot(time_axis, LongTS_ibgvoltot*1e-3 , 'k', label=CASE ,linewidth=0.7)
	plt.text(lgTS_ye+1.,32.,str(npy.round(npy.nanmean(LongTS_ibgvoltot)*1e-3 ,decimals=1)),color='k',size=8)
	plt.plot(time_axis_PIO, LongTS_OBS_icevol*1e-12 , 'g',label='Obs.' , linewidth=0.7 )
	plt.text(2013+1.,5.,str(npy.round(npy.nanmean(LongTS_OBS_icevol)*1e-12 ,decimals=1)),color='g',size=8)
	plt.xlim([lgTS_ys-1.,2018.])
	plt.ylim([0,40])
	plt.xticks(piomas_newlocsx,piomas_newlabelsx,size=5)
	#plt.setp(ax.get_xticklabels(),rotation=90)
	plt.setp(ax.get_xticklabels(),visible=False)
	plt.yticks(size=6)
	plt.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
	plt.ylabel('Ice volume \n'+r'(x$10^3$ $km^{3}$)',size=6)
	
	# Ice area 
	################
	ax=plt.subplot(xwind+2)
	plt.plot(time_axis, LongTS_ibgarea*1e-6 , 'k', linewidth=0.7  )
	plt.text(lgTS_ye+1.,17.,str(npy.round(npy.nanmean(LongTS_ibgarea)*1e-6 ,decimals=1)),color='k',size=8)
	plt.plot(time_axis_NSIDC, LongTS_OBS_iceare*1e-12 , 'g', label='Obs.', linewidth=0.7  )
	plt.text(2015+1.,1.,str(npy.round(npy.nanmean(LongTS_OBS_iceare)*1e-12 ,decimals=1)),color='g',size=8)
	plt.xlim([lgTS_ys-1.,2018.])
	plt.ylim([0,20])
	plt.xticks(nsidc_newlocsx,nsidc_newlabelsx,size=5)
	#plt.setp(ax.get_xticklabels(),rotation=90)
	plt.setp(ax.get_xticklabels(),visible=False)
	plt.yticks(size=6)
	plt.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
	plt.ylabel('Ice area \n'+r'(x$10^6$ $km^{2}$)',size=6)
	if lgTS_ys >= 1979 :
		plt.legend(loc='upper left')
		leg = plt.gca().get_legend()
		ltext = leg.get_texts()
		plt.setp(ltext, fontsize=5.)
	
	# Ice Extent 
	################
	ax=plt.subplot(xwind+3)
	curve_case = plt.plot(time_axis, LongTS_sice_ext*1e-12 , 'k', linewidth=0.7  )
	plt.text(lgTS_ye+1.,17.,str(npy.round(npy.nanmean(LongTS_sice_ext)*1e-12 ,decimals=1)),color='k',size=8)
	curve_obs = plt.plot(time_axis_NSIDC, LongTS_OBS_iceext*1e-12 , 'g', linewidth=0.7  )
	plt.text(2015+1.,1.,str(npy.round(npy.nanmean(LongTS_OBS_iceext)*1e-12 ,decimals=1)),color='g',size=8)
	plt.xlim([lgTS_ys-1.,2018.])
	plt.ylim([0,20])
	plt.xticks(nsidc_newlocsx,nsidc_newlabelsx,size=5)
	#plt.setp(ax.get_xticklabels(),rotation=90)
	plt.setp(ax.get_xticklabels(),visible=False)
	plt.yticks(size=6)
	plt.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
	plt.ylabel('Ice extent \n'+r'(x$10^6$ $km^{2}$)',size=6)
	
	# Mean Ice drift 
	################
	ax=plt.subplot(xwind+4)
	plt.title(' Sea-Ice drift',size=9)
	plt.plot(time_axis, LongTS_VEICE*100 , 'k', label='CRF box' , linewidth=0.7 )
	plt.plot(time_axis_obs, IABPObservations*100 , 'g', label='Obs.', linewidth=0.7  )
	plt.text(lgTS_ye+1.,15.   ,str(npy.round(npy.nanmean(LongTS_VEICE)*100. ,decimals=1)),color='k',size=8)
	plt.text(2011+1,1.  ,str(npy.round(npy.nanmean(IABPObservations)*100. ,decimals=1)),color='g',size=8)
	plt.xlim([lgTS_ys-1.,2018.])
	plt.ylim([0,20])
	plt.xticks(full_newlocsx,full_newlabelsx,size=5)
	plt.setp(ax.get_xticklabels(),rotation=90)
	plt.yticks(size=6)
	plt.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
	plt.ylabel('Ice drift '+r'(cm $s^{-1}$)',size=7)
	
	plt.legend(loc='upper left',ncol=3)
	leg = plt.gca().get_legend()
	ltext = leg.get_texts()
	plt.setp(ltext, fontsize=5.)
	plt.tight_layout()
	
	plt.savefig(CONFIG+'-'+CASE+'_ICEVolExtDrift-LGTS_y'+str(lgTS_ys)+'LASTy.png',dpi=400)
	
	# Plot the mean FWC & Ice drift 
	###############################

	plt.clf()
	xwind=210
	# Fresh water content (FWC)
	###########################
	ax=plt.subplot(xwind+1)
	ax2=ax.twinx()
	plt.title(CASE+' FWC \n '+str(lgtsclimyear),size=9)
	ax.text(lgTS_ye+0.5,17.   ,str(npy.round(npy.nanmean(LongTS_FWC)*1e-12 ,decimals=1)),color='k',size=8)
	ax.plot(time_axis, LongTS_FWC*1e-12 , 'k', label='CRF box model', linewidth=0.5  )
	ax.plot(time_axis_CRFFWC, LongTS_OBS_CRFFWC , 'g', linestyle='None', marker='8', markersize=3, label='CRF box Obs.', linewidth=0.5  )
	plt.xlim([lgTS_ys-1.,2018.])
	ax.set_ylim(12,24)
	ax.set_ylabel('FWC '+r'(x$10^3$ $km^3$)',size=6)
	ax.tick_params('y', labelsize=6)
	ax.set_xticks(newlocsx)
	#ax.tick_params('x', labelsize=5, labelrotation=90, )
	ax.set_xticklabels(list(newlabelsx),size=5,rotation=90)
	ax.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)

	ax2.plot(time_axis, LongTS_BIGFWC*1e-12 , 'r', label='Arctic box', linewidth=0.5  )
	ax2.set_ylim(70,130)
	ax2.tick_params('y', colors='r',labelsize=6)
	ax2.text(lgTS_ye-0.5,85. ,str(npy.round(npy.nanmean(LongTS_BIGFWC)*1e-12 ,decimals=1)),color='r',size=8)
	ax2.set_xlim([lgTS_ys-1.,2018.])
	
	ax.legend(loc='upper left',ncol=2)
	ax2.legend(loc='lower right',ncol=2)


	# Ekman pumping in the CRF box 
	###############################
	xwind=211
	ax=plt.subplot(xwind+1)
	plt.title(CASE+' Ekman pumping \n '+str(lgtsclimyear),size=9)
	ax.plot(time_axis, LongTS_WEkm*86400.*365. , 'k', label='Model CRF box', linewidth=0.5	)
	ax.plot(time_axis_CRFEkm, LongTS_OBS_CRFEkm , 'g', label='Obs. CRF box', linewidth=0.5	)
	plt.xlim([lgTS_ys-1.,2018.])
	ax.set_ylim(-40.,40.)
	ax.set_ylabel('WEk '+r'( m $year^{-1}$)',size=6)
	ax.tick_params('y', labelsize=6)
	ax.set_xticks(newlocsx)
	ax.set_xticklabels(list(newlabelsx),size=5,rotation=90)
	ax.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
	ax.legend(loc='upper left',ncol=2)


	plt.tight_layout()
	plt.savefig(CONFIG+'-'+CASE+'_FWC-WEK-LGTS_y'+str(lgTS_ys)+'LASTy.png',dpi=400)

	if MassSaltFLX :
		# Ice salt/mass flux in the CRF box 
		####################################
		plt.clf()
		xwind=210
		ax=plt.subplot(xwind+1)
		plt.title(CASE+' Salt flux \n '+str(lgtsclimyear),size=9)
		ax.plot(time_axis, LongTS_SFX*86400. , 'k', label='Model CRF box', linewidth=0.5  )
		plt.xlim([lgTS_ys-1.,2018.])
		#ax.set_ylim(-40.,40.)
		ax.set_ylabel('Sfx '+r'( kg $m^{-2}$ $day^{-1}$)',size=6)
		ax.tick_params('y', labelsize=6)
		ax.set_xticks(newlocsx)
		ax.set_xticklabels(list(newlabelsx),size=5,rotation=90)
		ax.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
		ax.legend(loc='upper left',ncol=2)

		ax=plt.subplot(xwind+2)
		plt.title(CASE+' Volume flux \n '+str(lgtsclimyear),size=9)
		ax.plot(time_axis, LongTS_VFX*86400. , 'k', label='Model CRF box', linewidth=0.5  )
		plt.xlim([lgTS_ys-1.,2018.])
		#ax.set_ylim(-40.,40.)
		ax.set_ylabel('Sfx '+r'( kg $m^{-2}$ $day^{-1}$)',size=6)
		ax.tick_params('y', labelsize=6)
		ax.set_xticks(newlocsx)
		ax.set_xticklabels(list(newlabelsx),size=5,rotation=90)
		ax.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
		ax.legend(loc='upper left',ncol=2)

		plt.tight_layout()
		plt.savefig(CONFIG+'-'+CASE+'_SVFX-LGTS_y'+str(lgTS_ys)+'LASTy.png',dpi=400)

	# Output diagnostics into a Netcdf file
	#######################################################################################################
	if NCDF_OUT:
		ds_outTS = xr.Dataset()
		# Define time axis
		ds_outTS.coords['time_axis']           = (('time_axis')          , time_axis.astype('float32'))
		ds_outTS.coords['time_axis_NSIDC']     = (('time_axis_NSIDC')    , time_axis_NSIDC.astype('float32'))
		ds_outTS.coords['time_axis_PIO']       = (('time_axis_PIO')      , time_axis_PIO.astype('float32'))
		ds_outTS.coords['time_axis_IABPobs']   = (('time_axis_IABPobs')  , time_axis_obs.astype('float32'))
		ds_outTS.coords['time_axis_CRFFWCobs'] = (('time_axis_CRFFWCobs'), time_axis_CRFFWC.astype('float32'))

		# Save diagnotics 
		ds_outTS['LongTS_ibgvoltot']= (('time_axis'), LongTS_ibgvoltot.astype('float32'))
		ds_outTS['LongTS_ibgvoltot'].attrs['long_name']='Sea-ice volume in the Arctic basin'
		ds_outTS['LongTS_ibgvoltot'].attrs['units']="km3"
		ds_outTS['LongTS_ibgarea']= (('time_axis'), LongTS_ibgarea.astype('float32'))
		ds_outTS['LongTS_ibgarea'].attrs['long_name']='Sea-ice area in the Arctic basin'
		ds_outTS['LongTS_ibgarea'].attrs['units']="km2"
		ds_outTS['LongTS_sice_ext']= (('time_axis'), LongTS_sice_ext.astype('float32'))
		ds_outTS['LongTS_sice_ext'].attrs['long_name']='Sea-ice extent in the Arctic basin'
		ds_outTS['LongTS_sice_ext'].attrs['units']="km2"
		ds_outTS['LongTS_VEICE']= (('time_axis'), LongTS_VEICE.astype('float32'))
		ds_outTS['LongTS_VEICE'].attrs['long_name']='Sea-ice drift calculated in the CRF box centered over the Beaufort Gyre'
		ds_outTS['LongTS_VEICE'].attrs['units']="m/s"
		ds_outTS['LongTS_FWC']= (('time_axis'), LongTS_FWC.astype('float32'))
		ds_outTS['LongTS_FWC'].attrs['long_name']='Freshwater content over the CRF box'
		ds_outTS['LongTS_FWC'].attrs['units']="km3"
		ds_outTS['LongTS_BIGFWC']= (('time_axis'), LongTS_BIGFWC.astype('float32'))
		ds_outTS['LongTS_BIGFWC'].attrs['long_name']='Freshwater content over the whole Arctic basin '
		ds_outTS['LongTS_BIGFWC'].attrs['units']="km3"
		ds_outTS['LongTS_WEkm']= (('time_axis'), (LongTS_WEkm*86400.*365.).astype('float32'))
		ds_outTS['LongTS_WEkm'].attrs['long_name']='Ekman pumping over the CRF box'
		ds_outTS['LongTS_WEkm'].attrs['units']="m/day"
		ds_outTS['LongTS_OBS_icevol']= (('time_axis_PIO'   ), LongTS_OBS_icevol.astype('float32')) 
		ds_outTS['LongTS_OBS_icevol'].attrs['long_name']='Sea-ice volume from PIOMAS reanalysis'
		ds_outTS['LongTS_OBS_icevol'].attrs['units']="km3"
		ds_outTS['LongTS_OBS_iceare']= (('time_axis_NSIDC' ), LongTS_OBS_iceare.astype('float32'))
		ds_outTS['LongTS_OBS_iceare'].attrs['long_name']='Sea-ice area from NSIDC '
		ds_outTS['LongTS_OBS_iceare'].attrs['units']="km2"
		ds_outTS['LongTS_OBS_iceext']= (('time_axis_NSIDC'), LongTS_OBS_iceext.astype('float32'))
		ds_outTS['LongTS_OBS_iceext'].attrs['long_name']='Sea-ice extent from NSIDC '
		ds_outTS['LongTS_OBS_iceext'].attrs['units']="km2"
		ds_outTS['LongTS_OBS_CRFFWC']= (('time_axis_CRFFWCobs'), LongTS_OBS_CRFFWC.astype('float32'))
		ds_outTS['LongTS_OBS_CRFFWC'].attrs['long_name']='Freshwater content over the CRF box from Proshutinsky et al. GRL2018 observations'
		ds_outTS['LongTS_OBS_CRFFWC'].attrs['units']="m/s"
		ds_outTS['IABPObservations']= (('time_axis_obs'), IABPObservations.astype('float32')) 
		ds_outTS['IABPObservations'].attrs['long_name']='Sea-ice drift from IABP observation system'
		ds_outTS['IABPObservations'].attrs['units']="m/s"

		if MassSaltFLX :
			ds_outTS['LongTS_SFX']= (('time_axis'), LongTS_SFX.astype('float32'))
			ds_outTS['LongTS_SFX'].attrs['long_name']='Ice salt flux'
			ds_outTS['LongTS_SFX'].attrs['units']="kg/m2/s"
			ds_outTS['LongTS_VFX']= (('time_axis'), LongTS_VFX.astype('float32'))
			ds_outTS['LongTS_VFX'].attrs['long_name']='Ice mass flux'
			ds_outTS['LongTS_VFX'].attrs['units']="kg/m2/s"

		# Write the NetCDF file 
		ds_outTS.attrs['History'] = "Diagnostics have been calculated using the Arctic monitoring tool "
		ds_outTS.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
		nc_f = './NETCDF/'+CONFIG+'-'+CASE+'_ICEFWCWEK-LGTS_'+'y'+str(lgTS_ys)+'LASTy.nc'
		ds_outTS.to_netcdf(nc_f,engine='netcdf4')
