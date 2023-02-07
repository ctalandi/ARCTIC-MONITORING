#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import numpy as npy
import matplotlib.pylab as plt
from matplotlib import rcParams
from netCDF4 import Dataset
import matplotlib as mpl
import sys 
from checkfile import *
from CREG_intquant_func import *
import subprocess

#rcParams['text.usetex']=True
rcParams['font.family']='serif'

s_year=XXSYEAXX
e_year=XXEYEAXX
lgTS_ys=XXLGTSSXX
lgTS_ye=XXLGTSEXX
xiosfreq=XXXIOSFREQXX
NCDF_OUT=XXNCDFOUTXX
main_dir='./'
CONFIG='XXCONFXX'   ; CASE='XXCASEXX'     ;  CASE2='XXCASE2XX'
CONFCASE=CONFIG+'-'+CASE   ; CONFCASE2=CONFIG+'-'+CASE2
data_dir=main_dir+CONFIG+'/'+CONFCASE+'-MEAN/'    ; data_dir2=main_dir+CONFIG+'/'+CONFCASE2+'-MEAN/'
grid_dir=main_dir+CONFIG+'/GRID/'

DIR_FIG_OUT='./'

# Infos concernant les climatologies sur la periode de la simulation
climyear=str(s_year)+str(e_year)

print 
print '                              Configuration :' , CONFCASE
print '                              Period        :' , str(s_year),' - ',str(e_year)
print 


#####################################################################
var_temp={'name':"votemper",'units':u"degC",'_FillValue': 0.,'fext':"T",'igrd':1 ,'ze3':"e3t_0",'gdep':"gdept",'long_name':"Temperature"}
var_sali={'name':"vosaline",'units':u"psu" ,'_FillValue': 0.,'fext':"S",'igrd':1 ,'ze3':"e3t_0",'gdep':"gdept",'long_name':"Salinity"}
var_crty={'name':"vomecrty",'units':u"m/s" ,'_FillValue': 0.,'fext':"V",'igrd':3 ,'ze3':"e3v_0",'gdep':"gdepv",'long_name':"Meridional velocity"}
var_hthi={'name':"sivolu"  ,'units':u"m"   ,'_FillValue': 0.,'fext':"I",'igrd':1 ,'ze3':"e3t_0",'gdep':"gdept",'long_name':"Ice thickness "}
var_Wcur={'name':"socurl"  ,'units':u"m/s" ,'_FillValue': 0.,'fext':"Xi",'igrd':1 ,'ze3':"e3f_0",'gdep':"gdepw",'long_name':"Surface ocean stress curl "}
var_ISfx={'name':"sfx"     ,'units':u"kg/m2/s" ,'_FillValue': 0.,'fext':"sfx",'igrd':1 ,'ze3':"e3t_0",'gdep':"gdept",'long_name':"Surface ice salt flux "}
var_IVfx={'name':"fmmflx"  ,'units':u"kg/m2/s" ,'_FillValue': 0.,'fext':"fmmflx",'igrd':1 ,'ze3':"e3t_0",'gdep':"gdept",'long_name':"Surface ice mass flux"}
MassSaltFLX=False

########################################
# Read GRID 
########################################
#------------------------------------------------------------------------------------------------------------------------

# Read only once for the same section type
locpath=grid_dir
locfile=CONFCASE+'_mask.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	fieldmask=Dataset(locpath+locfile)
	tmask = npy.squeeze(fieldmask.variables['tmask'])
	fmask = npy.squeeze(fieldmask.variables['fmask'])

locfile=CONFCASE+'_coordinates.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	fieldcoor=Dataset(locpath+locfile)
	lon = npy.squeeze(fieldcoor.variables['glamt'])
	lat = npy.squeeze(fieldcoor.variables['gphit'])
	ffCor = npy.squeeze(fieldcoor.variables['ff_f'])
	#ffCor = npy.squeeze(fieldcoor.variables['ff'])

# Read appropriate vertical scale factor 
infield=var_sali
locpath=grid_dir
locfile=CONFCASE+'_mesh_zgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	infield=var_temp
	fieldzmesh=Dataset(locpath+locfile)
	ze3 = npy.squeeze(fieldzmesh.variables[infield['ze3']])

locpath=grid_dir
locfile=CONFCASE+'_mesh_hgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	fieldhmesh=Dataset(locpath+locfile)
	ze1t= npy.squeeze(fieldhmesh.variables['e1t'])
	ze2t= npy.squeeze(fieldhmesh.variables['e2t'])
	ze1f= npy.squeeze(fieldhmesh.variables['e1f'])
	ze2f= npy.squeeze(fieldhmesh.variables['e2f'])

e1te2t= ze1t*ze2t*tmask[0,:,:]
e1fe2f= ze1f*ze2f*fmask[0,:,:]

time_dim=(e_year-s_year+1)*12

#------------------------------------------------------------------------------------------------------------------------
########################################
# Read input  data 
########################################
infield=var_sali
print "                                            "
print "                                        >>>>>>>>>  The curent variable processed is :        ", infield['fext']
print "                                            "
Ar_size=(time_dim,tmask.shape[0],lon.shape[0],lon.shape[1])
Sdata_read=npy.zeros(Ar_size)
year=s_year
cur_month=0
while cur_month <= 11 :
	str_month='m0'+str(cur_month+1) if cur_month <= 8 else 'm'+str(cur_month+1)
        locpath=data_dir+'/'+str(s_year)+'/'+xiosfreq+'/'
        locfile=CONFIG+'-'+CASE+'_y'+str(year)+str_month+'.'+xiosfreq+'_gridT.nc'
	if chkfile(locpath+locfile) :
        	field = Dataset(locpath+locfile)
        	Sdata_read[cur_month,:,:,:] = npy.squeeze(field.variables[infield['name']])
	cur_month = cur_month + 1 

print "                                            "
print "                                        >>>>>>>>>  The curent variable processed is :        at_i(:,:)", 
print "                                            "
Ar_size=(time_dim,lon.shape[0],lon.shape[1])
Ifrdata_read=npy.zeros(Ar_size)
Ivedata_read=npy.zeros(Ar_size)
Ithdata_read=npy.zeros(Ar_size)
year=s_year
cur_month=0
while cur_month <= 11 :
	str_month='m0'+str(cur_month+1) if cur_month <= 8 else 'm'+str(cur_month+1)
        locpath=data_dir+'/'+str(s_year)+'/'+xiosfreq+'/'
        locfile=CONFIG+'-'+CASE+'_y'+str(year)+str_month+'.'+xiosfreq+'_icemod.nc'
	if chkfile(locpath+locfile) :
        	field = Dataset(locpath+locfile)
        	Ifrdata_read[cur_month,:,:] = npy.squeeze(field.variables['siconc'])
        	Ivedata_read[cur_month,:,:] = npy.squeeze(field.variables['sivelo'])
        	Ithdata_read[cur_month,:,:] = npy.squeeze(field.variables['sivolu'])
	else:
        	Ifrdata_read[cur_month,:,:] = npy.nan
        	Ivedata_read[cur_month,:,:] = npy.nan
        	Ithdata_read[cur_month,:,:] = npy.nan
	cur_month = cur_month + 1 

infield=var_Wcur
print "                                            "
print "                                        >>>>>>>>>  The curent variable processed is :        ", infield['fext']
print "                                            "
Ar_size=(time_dim,lon.shape[0],lon.shape[1])
Scurldata_read=npy.zeros(Ar_size)
year=s_year
cur_month=0
locpath=data_dir+'/'+str(s_year)+'/'+xiosfreq+'/'
locfile=CONFIG+'-'+CASE+'_y'+str(year)+'.'+xiosfreq+'_SurOceCurl.nc'
if chkfile(locpath+locfile) :
	field = Dataset(locpath+locfile)
	Scurldata_read[:,:,:] = npy.squeeze(field.variables[infield['name']])
	Scurldata_read[:,:,:] = npy.where( npy.abs(Scurldata_read) > 1e10 , 0., Scurldata_read )
	Scurldata_read[:,:,:] = npy.where( npy.isnan(Scurldata_read) , 0., Scurldata_read )
	Scurldata_read[:,:,:] = npy.where( fmask[0,:,:] < 1. , 0., Scurldata_read )

if MassSaltFLX: 
	infield=var_ISfx
	print "                                            "
	print "                                        >>>>>>>>>  The curent variable processed is :        ", infield['fext']
	print "                                            "
	Ar_size=(time_dim,lon.shape[0],lon.shape[1])
	SFXdata_read=npy.zeros(Ar_size)
	year=s_year
	cur_month=0
	while cur_month <= 11 :
		str_month='m0'+str(cur_month+1) if cur_month <= 8 else 'm'+str(cur_month+1)
	        locpath=data_dir+'/'+str(s_year)+'/'+xiosfreq+'/'
	        locfile=CONFIG+'-'+CASE+'_y'+str(year)+str_month+'.'+xiosfreq+'_flxT.nc'
		if chkfile(locpath+locfile) :
			field = Dataset(locpath+locfile)
			SFXdata_read[cur_month,:,:] = npy.squeeze(field.variables[infield['name']])
		else :
			SFXdata_read[cur_month,:,:] = npy.nan
	
		cur_month = cur_month + 1
	
	SFXdata_read = npy.where( npy.abs(SFXdata_read) > 1e10 , 0., SFXdata_read )
	SFXdata_read = npy.where( tmask[0,:,:] < 1. , 0., SFXdata_read )
	
	infield=var_IVfx
	print "                                            "
	print "                                        >>>>>>>>>  The curent variable processed is :        ", infield['fext']
	print "                                            "
	Ar_size=(time_dim,lon.shape[0],lon.shape[1])
	VFXdata_read=npy.zeros(Ar_size)
	year=s_year
	cur_month=0
	while cur_month <= 11 :
		str_month='m0'+str(cur_month+1) if cur_month <= 8 else 'm'+str(cur_month+1)
	        locpath=data_dir+'/'+str(s_year)+'/'+xiosfreq+'/'
	        locfile=CONFIG+'-'+CASE+'_y'+str(year)+str_month+'.'+xiosfreq+'_flxT.nc'
		if chkfile(locpath+locfile) :
			field = Dataset(locpath+locfile)
			VFXdata_read[cur_month,:,:] = npy.squeeze(field.variables[infield['name']])
		else :
			VFXdata_read[cur_month,:,:] = npy.nan
	
		cur_month = cur_month + 1
	
	#VFXdata_read = npy.where( npy.abs(VFXdata_read) > 1e10 , 0., VFXdata_read )
	VFXdata_read = npy.where( tmask[0,:,:] < 1. , 0., VFXdata_read )


################################################################################################################
################################################################################################################
##################################  FWC, ICE VOLUME & CONCENTRATION etc ...Diags ###############################
################################################################################################################
################################################################################################################
print 
print '			FWC, ICE VOLUME, AREA, EXTENSION, DRIFT , WIND STRESS MODULE Diags'
print 

plt.clf()

Sref=34.80 
ze33D=ze3[:,:,:]*tmask[:,:,:]
ze33D=npy.ma.masked_where( tmask[:,:,:] == 0., ze33D )
fwc2D = npy.zeros((time_dim,lon.shape[0],lon.shape[1]))

#print a specific point to debug within the Beaufort Gyre
dbg=False
if CONFIG == 'CREG025.L75' :
	idbg=185 ; jdbg=515 ; kdbg=50   # CREG025.L75
elif CONFIG == 'CREG12.L75' :
	idbg=562 ; jdbg=1558 ; kdbg=50   # CREG12.L75

for ti in set(npy.arange(12)) :
    ze33D_msk = ze33D.copy()
    ze33D_msk=npy.ma.masked_where(npy.squeeze(Sdata_read[ti,:,:,:]) > Sref, ze33D_msk )
    for jk in set(npy.arange(75)) :
        fwc2D[ti,:,:] = fwc2D[ti,:,:] + ( Sref - Sdata_read[ti,jk,:,:] ) / Sref * ze33D_msk[jk,:,:] 

    if dbg:
    	print 'ze33Dtime_msk[ti,0:kdbg,jdbg,idbg]',ze33Dtime_msk[ti,0:kdbg,jdbg,idbg], Sdata_read[ti,0:kdbg,jdbg,idbg]
    	print 'fwc2D[ti,jdbg,idbg]',fwc2D[ti,jdbg,idbg]

tmask2Dtime=npy.tile(tmask[0,:,:],((e_year-s_year+1)*12,1,1))
fwc2D=npy.ma.masked_where((tmask2Dtime == 0),fwc2D)

# Beaufort Gyre as defined in J. Scott et al. 2015
e1te2tBG=e1te2t.copy()
area=npy.tile(e1te2t,(time_dim,1,1))
tmskBFG=npy.ones((tmask.shape[1],tmask.shape[2]))
tmskBFG[npy.where(lon[:,:] > -130.)]=0.
tmskBFG[npy.where(lon[:,:] < -170.)]=0.
tmskBFG[npy.where(lat[:,:] >  80.5)]=0.
tmskBFG[npy.where(lat[:,:] <  70.5)]=0.
e1te2tBG[:,:] = e1te2tBG * tmskBFG
e1te2tBGtime=npy.tile(e1te2tBG[:,:],(time_dim,1,1))

FWC = npy.sum(fwc2D*e1te2tBGtime,axis=(1,2))
Vice_mean = npy.sum(Ivedata_read*e1te2tBGtime,axis=(1,2))/npy.sum(e1te2tBGtime,axis=(1,2))

area[npy.where(Ifrdata_read < 0.15 )] = 0.e0
sice_ext = npy.sum(area[:,:,:], axis=(1,2))

# Still the same area but now for the vorticity field or F-point
e1fe2fBG=e1fe2f.copy()
fmskBFG=npy.ones((fmask.shape[1],fmask.shape[2]))
fmskBFG[npy.where(lon[:,:] > -130.)]=0.
fmskBFG[npy.where(lon[:,:] < -170.)]=0.
fmskBFG[npy.where(lat[:,:] >  80.5)]=0.
fmskBFG[npy.where(lat[:,:] <  70.5)]=0.
e1fe2fBG[:,:] = e1fe2fBG * fmskBFG
e1fe2fBGtime=npy.tile(e1fe2fBG[:,:],(time_dim,1,1))
ffCortime=npy.tile(ffCor[:,:],(time_dim,1,1))

rho=1027.5  # Same value as in Gianluca et al. JPO2017, Observations of Seasonal Upwelling and Downwelling in the Beaufort Sea Mediated by Sea Ice
WEkm_mean = npy.sum(Scurldata_read*e1fe2fBGtime/(rho*ffCortime),axis=(1,2))/npy.sum(e1fe2fBGtime,axis=(1,2))

if MassSaltFLX :
	SFX_mean = npy.sum(SFXdata_read*e1te2tBGtime,axis=(1,2))/npy.sum(e1te2tBGtime,axis=(1,2))
	VFX_mean = npy.sum(VFXdata_read*e1te2tBGtime,axis=(1,2))/npy.sum(e1te2tBGtime,axis=(1,2))

npy.save(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_FWCTS_y'+str(s_year),npy.array(FWC))
npy.save(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_ICEExtTS_y'+str(s_year),npy.array(sice_ext))
npy.save(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_VICETS_y'+str(s_year),npy.array(Vice_mean))
npy.save(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_WEKMATS_y'+str(s_year),npy.array(WEkm_mean))
if MassSaltFLX :
	npy.save(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_ICESFXTS_y'+str(s_year),npy.array(SFX_mean))
	npy.save(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_ICEVFXTS_y'+str(s_year),npy.array(VFX_mean))

# Beaufort Gyre, define a bigger box than the previous one
e1te2tBGBIG=e1te2t.copy()
tmskBFGBIG=npy.ones((tmask.shape[1],tmask.shape[2]))
# Define a proper Arctic box following Haine et al. GPC2015
# including all Arctic area + CAA + north to Davis, FRAM and Bering straits as east to BSO
tmskBFGBIG[npy.where(lat[:,:] <  62.)] = 0.
# Remove the Baffin bay area
tmskBFGBIG[npy.where((lon[:,:] >= -90.) & (lon[:,:] < -45.) & (lat[:,:] <= 80.) )]=0.
# Adjust the Baffin bay area north to Davis Strait 66.25N
tmskBFGBIG[npy.where((lon[:,:] >= -70.) & (lon[:,:] < -45.) & (lat[:,:] >= 62.) & (lat[:,:] < 66.25) )] = 0.
tmskBFGBIG[npy.where((lon[:,:] >= -70.) & (lon[:,:] < -65.) & (lat[:,:] >= 66.25) & (lat[:,:] < 67.) )] = 0.
# Remove all GIN seas area
tmskBFGBIG[npy.where((lon[:,:] < 17.) & (lon[:,:] > -45.) & (lat[:,:] < 79.))] = 0.
# Remove Hudson area
tmskBFGBIG[npy.where((lon[:,:] < -70.) & (lon[:,:] > -82. ) & (lat[:,:] < 70.5))] = 0.
tmskBFGBIG[npy.where((lon[:,:] <= -82.) & (lon[:,:] > -100.) & (lat[:,:] < 67.))] = 0.
# Adriatic area 
tmskBFGBIG[npy.where((lon[:,:] < 30.) & (lon[:,:] > 17.) & (lat[:,:] < 68.))] = 0.

e1te2tBGBIG[:,:] = e1te2tBGBIG * tmskBFGBIG
e1te2tBGBIGtime=npy.tile(e1te2tBGBIG[:,:],(time_dim,1,1))

BigFWC = npy.sum(fwc2D*e1te2tBGBIGtime,axis=(1,2))

npy.save(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_BIGFWCTS_y'+str(s_year),npy.array(BigFWC))


######################################################################################################
######################################################################################################

################################################################################################################
################################################################################################################
#################################   Beaufort Gyre FWC, Ice thick & Drift Time Series  ##########################
################################################################################################################
################################################################################################################

if lgTS_ye-lgTS_ys+1 > 1 :
	print
	print " 			##################################################################  " 
	print " 			##################################################################  " 
	print " 			######### PLOT DIAGS LONG TIME-SERIES WITHIN BEAUFORT GYRE #######  " 
	print " 			##################################################################  " 
	print " 			##################################################################  " 
	print

	lgtstime_dim=(lgTS_ye-lgTS_ys+1)*12

	LongTS_FWC= []     ; LongTS_sice_ext=[]   ; LongTS_ibgvoltot=[]  ; LongTS_ibgarea=[]
	LongTS_VEICE= []   ; LongTS_HIICE=[]      ; LongTS_BIGFWC= []  ; LongTS_WEkm=[] ;  LongTS_SFX= []   ; LongTS_VFX=[]      

	# Start to read all yearly files
	################################
	lgts_year=lgTS_ys    ;    t_months=(npy.arange(12)*30.+15.)/365.   ;   start = 1
	while  lgts_year <= lgTS_ye  :
		print 
		print " 			>>>>   Read year:", lgts_year

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
			field = Dataset(locpath+locfile)
			year_ibgvoltot = npy.squeeze(field.variables['ibgvoltot'])
			year_ibgarea = npy.squeeze(field.variables['ibgarea'])
                	LongTS_ibgvoltot=npy.append(LongTS_ibgvoltot,year_ibgvoltot)
                	LongTS_ibgarea=npy.append(LongTS_ibgarea,year_ibgarea)
		else:
			ibgvoltot, ibgarea = CAL_ICE_VOL_AREA(CONFIG,CASE,lgts_year,data_dir,xiosfreq, dom_area=e1te2t.copy(),tmask2D=tmask[0,:,:].squeeze())
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

	# Plot the FWC & Sea-ice time-series over SEVERAL YEARS
	########################################################
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
        newlocsx_nsidc  = npy.array(time_axis_NSIDC,'f')
        newlabelsx_nsidc = npy.array(time_axis_NSIDC,'i')
	nsidc_newlocsx=npy.array(npy.append(newlocsx,npy.unique(newlabelsx_nsidc)),'f')
	nsidc_newlabelsx=npy.array(npy.append(newlabelsx,npy.unique(newlabelsx_nsidc)),'i')

	# Ekman pumping time axis
	time_axis_CRFEkm=[]
	for pyear in npy.arange(2003,2015):
        	y_years=npy.tile(pyear,12)+t_months
        	time_axis_CRFEkm=npy.append(time_axis_CRFEkm,y_years)

		pyear+=1
        newlocsx_ekman  = npy.array(time_axis_CRFEkm,'f')
        newlabelsx_ekman = npy.array(time_axis_CRFEkm,'i')
	ekman_newlocsx=npy.array(npy.append(newlocsx,npy.unique(newlabelsx_ekman)),'f')
	ekman_newlabelsx=npy.array(npy.append(newlabelsx,npy.unique(newlabelsx_ekman)),'i')
	
	print('new obs time',newlocsx_ekman )
	
	lgtsclimyear=str(lgTS_ys)+str(lgTS_ye)

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
	plt.ylabel('Ice volume \n'+r'(x$10^3$ $km^{3}}$)',size=6)
	
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
	plt.ylabel('Ice area \n'+r'(x$10^6$ $km^{2}}$)',size=6)
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
	plt.ylabel('Ice extent \n'+r'(x$10^6$ $km^{2}}$)',size=6)
	
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
	
	plt.savefig(CONFIG+'-'+CASE+'_ICEVolExtDrift-LGTS_y'+str(lgTS_ys)+'LASTy.pdf')
	
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
	ax.plot(time_axis, LongTS_WEkm*86400.*365. , 'k', label='Model CRF box', linewidth=0.5  )
	ax.plot(time_axis_CRFEkm, LongTS_OBS_CRFEkm , 'g', label='Obs. CRF box', linewidth=0.5  )
	plt.xlim([lgTS_ys-1.,2018.])
	ax.set_ylim(-40.,40.)
	ax.set_ylabel('WEk '+r'( m $year^{-1}$)',size=6)
	ax.tick_params('y', labelsize=6)
	ax.set_xticks(newlocsx)
	ax.set_xticklabels(list(newlabelsx),size=5,rotation=90)
	ax.grid(True, linestyle='--', which='both', color='grey', alpha=0.50)
	ax.legend(loc='upper left',ncol=2)


	plt.tight_layout()
	plt.savefig(CONFIG+'-'+CASE+'_FWC-WEK-LGTS_y'+str(lgTS_ys)+'LASTy.pdf',dpi=400)

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
		plt.savefig(CONFIG+'-'+CASE+'_SVFX-LGTS_y'+str(lgTS_ys)+'LASTy.pdf',dpi=400)

        if NCDF_OUT:
		# FWC field 
		#######################
		cmd_ddate="date"  ;  get_output = subprocess.check_output(cmd_ddate)
                nc_f = './NETCDF/'+CONFIG+'-'+CASE+'_ICEFWCWEK-LGTS_'+'y'+str(lgTS_ys)+'LASTy.nc'
                w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
                w_nc_fid.description = "Diagnostics have been calculated using the Arctic monitoring tool "
		w_nc_fid.date=get_output.decode("utf-8")
                w_nc_fid.createDimension('time_axis_mod', time_axis.shape[0])
                w_nc_fid.createDimension('time_axis_NSIDC', time_axis_NSIDC.shape[0])
                w_nc_fid.createDimension('time_axis_PIO', time_axis_PIO.shape[0])
                w_nc_fid.createDimension('time_axis_IABPobs', time_axis_obs.shape[0])
                w_nc_fid.createDimension('time_axis_CRFFWCobs', time_axis_CRFFWC.shape[0])
                #w_nc_fid.createDimension('time_counter', None)

                w_nc_var = w_nc_fid.createVariable('LongTS_ibgvoltot', 'f4', ('time_axis_mod'))
                w_nc_var.long_name='Sea-ice volume in the Arctic basin'
                w_nc_var.units="km3"
                w_nc_fid.variables['LongTS_ibgvoltot'][:] = LongTS_ibgvoltot

                w_nc_var = w_nc_fid.createVariable('LongTS_ibgarea', 'f4', ('time_axis_mod'))
                w_nc_var.long_name='Sea-ice area in the Arctic basin'
                w_nc_var.units="km2"
                w_nc_fid.variables['LongTS_ibgarea'][:] = LongTS_ibgarea

                w_nc_var = w_nc_fid.createVariable('LongTS_sice_ext', 'f4', ('time_axis_mod'))
                w_nc_var.long_name='Sea-ice extent in the Arctic basin'
                w_nc_var.units="km2"
                w_nc_fid.variables['LongTS_sice_ext'][:] = LongTS_sice_ext

                w_nc_var = w_nc_fid.createVariable('LongTS_VEICE', 'f4', ('time_axis_mod'))
                w_nc_var.long_name='Sea-ice drift calculated in the CRF box centered over the Beaufort Gyre'
                w_nc_var.units="m/s"
                w_nc_fid.variables['LongTS_VEICE'][:] = LongTS_VEICE

                w_nc_var = w_nc_fid.createVariable('LongTS_FWC', 'f4', ('time_axis_mod'))
                w_nc_var.long_name='Freshwater content over the CRF box'
                w_nc_var.units="km3"
                w_nc_fid.variables['LongTS_FWC'][:] = LongTS_FWC

                w_nc_var = w_nc_fid.createVariable('LongTS_BIGFWC', 'f4', ('time_axis_mod'))
                w_nc_var.long_name='Freshwater content over the whole Arctic basin '
                w_nc_var.units="km3"
                w_nc_fid.variables['LongTS_BIGFWC'][:] = LongTS_BIGFWC

                w_nc_var = w_nc_fid.createVariable('LongTS_WEkm', 'f4', ('time_axis_mod'))
                w_nc_var.long_name='Ekman pumping over the CRF box'
                w_nc_var.units="m/day"
                w_nc_fid.variables['LongTS_WEkm'][:] = LongTS_WEkm*86400.*365.

		if MassSaltFLX :
                	w_nc_var = w_nc_fid.createVariable('LongTS_SFX', 'f4', ('time_axis_mod'))
                	w_nc_var.long_name='Ice salt flux'
                	w_nc_var.units="kg/m2/s"
                	w_nc_fid.variables['LongTS_SFX'][:] = LongTS_SFX

                	w_nc_var = w_nc_fid.createVariable('LongTS_VFX', 'f4', ('time_axis_mod'))
                	w_nc_var.long_name='Ice mass flux'
                	w_nc_var.units="kg/m2/s"
                	w_nc_fid.variables['LongTS_VFX'][:] = LongTS_VFX

                w_nc_var = w_nc_fid.createVariable('LongTS_OBS_icevol', 'f4', ('time_axis_PIO'))
                w_nc_var.long_name='Sea-ice volume from PIOMAS model'
                w_nc_var.units="km3"
                w_nc_fid.variables['LongTS_OBS_icevol'][:] = LongTS_OBS_icevol

                w_nc_var = w_nc_fid.createVariable('LongTS_OBS_iceare', 'f4', ('time_axis_NSIDC'))
                w_nc_var.long_name='Sea-ice area from NSIDC '
                w_nc_var.units="km2"
                w_nc_fid.variables['LongTS_OBS_iceare'][:] = LongTS_OBS_iceare

                w_nc_var = w_nc_fid.createVariable('LongTS_OBS_iceext', 'f4', ('time_axis_NSIDC'))
                w_nc_var.long_name='Sea-ice extent from NSIDC '
                w_nc_var.units="km2"
                w_nc_fid.variables['LongTS_OBS_iceext'][:] = LongTS_OBS_iceext

                w_nc_var = w_nc_fid.createVariable('IABPObservations', 'f4', ('time_axis_IABPobs'))
                w_nc_var.long_name='Sea-ice drift from IABP observation system'
                w_nc_var.units="m/s"
                w_nc_fid.variables['IABPObservations'][:] = IABPObservations

                w_nc_var = w_nc_fid.createVariable('LongTS_OBS_CRFFWC', 'f4', ('time_axis_CRFFWCobs'))
                w_nc_var.long_name='Freshwater content over the CRF box from Proshutinsky et al. GRL2018 observations'
                w_nc_var.units="m/s"
                w_nc_fid.variables['LongTS_OBS_CRFFWC'][:] = LongTS_OBS_CRFFWC
                
                w_nc_fid.close()  # close the file





