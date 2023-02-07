import numpy as npy
import scipy.io as sio
from netCDF4 import Dataset
from checkfile import *

def READ_MOD_LGTS_DATA(locpath,locfile,zLongTS,npz=None,zvarname=None) :
	# Function dedicated to read model data in numpy npy format (Default)
	if chkfile(locpath+locfile) :
		year_fld=npy.load(locpath+locfile,mmap_mode='r')
		zLongTS=npy.append(zLongTS,year_fld) if npz == None else npy.append(zLongTS,year_fld[zvarname])
	else:
		zLongTS=npy.append(zLongTS,npy.arange(12)+npy.nan)

	return  zLongTS


def READ_OBS_LGTS_DATA(CONFIG,lgTS_ys,lgTS_ye) :
	# Function dedicated to read observations data that don't change 
	t_months=(npy.arange(12)*30.+15.)/365.

	# Ice volume obs
        locpath='./'
	if CONFIG == 'CREG025.L75' :
		locfile = 'PIOMAS_icevol_maskedBeringSea_interpCREG025.L75_1-12_1979-2018.nc'
	elif CONFIG == 'CREG12.L75' :
		locfile = 'PIOMAS_icevol_maskedBeringSea_interpCREG12.L75_1-12_1979-2018.nc'
	if chkfile(locpath+locfile) :
		field = Dataset(locpath+locfile)
		LongTS_OBS_icevol = npy.squeeze(field.variables['icevol-BS'])
	else:
		LongTS_OBS_icevol = npy.arange((2018-1979+1)*12)+npy.nan
	# Set the time axis for PIOMAS observations
	lgts_year=1979    ;      start = 1
	while  lgts_year <= 2018  :
        	y_years=npy.tile(lgts_year,12)+t_months
        	if start == 1:
        	        time_axis_PIO=y_years
        	        start=0
        	else:
        	        time_axis_PIO=npy.append(time_axis_PIO,y_years)
		lgts_year+=1

	
	# Ice extent obs
	# Data start in November 1978, so there is a 2 months shift
        locpath='./'
	locfile='NSIDC_ice_area_and_extent_maskBeringSea_fullPoleGap.nc'
	if chkfile(locpath+locfile) :
		field = Dataset(locpath+locfile)
		LongTS_OBS_iceext = npy.squeeze(field.variables['ice_extent'][2::])
		LongTS_OBS_iceare = npy.squeeze(field.variables['ice_area'][2::])
	else:
		LongTS_OBS_iceext = npy.arange((2015-1979+1)*12)+npy.nan
		LongTS_OBS_iceare = npy.arange((2015-1979+1)*12)+npy.nan
	# Set the time axis for NSIDC observations
	lgts_year=1979    ;      start = 1
	while  lgts_year <= 2015  :
        	y_years=npy.tile(lgts_year,12)+t_months
        	if start == 1:
        	        time_axis_NSIDC=y_years
        	        start=0
        	else:
        	        time_axis_NSIDC=npy.append(time_axis_NSIDC,y_years)
		lgts_year+=1


	# Ice drift from IABP
	# Data start in 18/01/1979, with 784 Buoys and 2 smapling / day : 0 & 12 
        locpath='./'
	locfile='ice_drift_BG_1979-2011.mat'
	if chkfile(locpath+locfile) :
		IABPObservations_read = sio.loadmat(locpath+locfile,squeeze_me=True)
		IABPObservations = npy.array(IABPObservations_read['time_series'])
	else:
		IABPObservations = npy.arange(396)+npy.nan

	# Set the time axis for observations
	lgts_year=1979    ;      start = 1
	while  lgts_year <= 2011  :
        	y_years=npy.tile(lgts_year,12)+t_months
        	if start == 1:
        	        time_axis_obs=y_years
        	        start=0
        	else:
        	        time_axis_obs=npy.append(time_axis_obs,y_years)
		lgts_year+=1

	return LongTS_OBS_icevol, LongTS_OBS_iceext, LongTS_OBS_iceare, IABPObservations, time_axis_obs, time_axis_PIO, time_axis_NSIDC

def CAL_ICE_VOL_AREA(CONFIG,CASE,lgts_year,data_dir,xiosfreq,dom_area,tmask2D) :

	Ar_size=(12,dom_area.shape[0],dom_area.shape[1]) ;  dom_areatime=npy.tile(dom_area,(12,1,1))   ;  tmask2Dtime=npy.tile(tmask2D,(12,1,1))
	Ithdata_read=npy.zeros(Ar_size)  ;     Ifrdata_read=npy.zeros(Ar_size)
	cur_month=0
	while cur_month <= 11 :
		str_month='m0'+str(cur_month+1) if cur_month <= 8 else 'm'+str(cur_month+1)
	        locpath=data_dir+'/'+str(lgts_year)+'/'+xiosfreq+'/'
	        locfile=CONFIG+'-'+CASE+'_y'+str(lgts_year)+str_month+'.'+xiosfreq+'_icemod.nc'
		if chkfile(locpath+locfile) :
	        	field = Dataset(locpath+locfile)
	        	Ithdata_read[cur_month,:,:] = npy.where(npy.squeeze(field.variables['sivolu']) > 1e18 , 0, npy.squeeze(field.variables['sivolu']))
	        	Ifrdata_read[cur_month,:,:] = npy.where(npy.squeeze(field.variables['sivolu']) > 1e18 , 0, npy.squeeze(field.variables['siconc']))
		else:
	        	Ithdata_read[cur_month,:,:] = npy.nan
	        	Ifrdata_read[cur_month,:,:] = npy.nan
		cur_month = cur_month + 1 

	ice_volume = npy.sum( Ithdata_read[:,:,:] * dom_areatime[:,:,:] * tmask2Dtime[:,:,:], axis=(1,2) )
	ice_area   = npy.sum( Ifrdata_read[:,:,:] * dom_areatime[:,:,:] * tmask2Dtime[:,:,:], axis=(1,2) )

	return ice_volume, ice_area


def READ_OBS_LGTS_CRFFWC(lgTS_ys,lgTS_ye) :
	# Function dedicated to read observations data that don't change 
	t_months=(npy.arange(12)*30.+15.)/365.

	# CRF box mean FWC Obs.
        locpath='./'
	locfile = 'BeaufortGyreFWC-Obs-Proshutinsky_GRL2018_y2003-2017.nc'
	if chkfile(locpath+locfile) :
		field = Dataset(locpath+locfile)
		LongTS_OBS_FWC = npy.squeeze(field.variables['CRFBGFWC_mean'])
		time_axis_FWC = npy.squeeze(field.variables['time_obs'])
	else:
		LongTS_OBS_FWC = npy.arange((2017-2003+1))+npy.nan
		time_axis_FWC = npy.arange((2017-2003+1))+npy.nan

	return  LongTS_OBS_FWC, time_axis_FWC


def READ_OBS_LGTS_CRFEkm(lgTS_ys,lgTS_ye) :
        # Function dedicated to read observations data that don't change 
        
        # CRF box mean Ekman pumping Obs.
        locpath='./'
        locfile = 'ArcticEkmanPumping_MonthlyMean.nc'
        if chkfile(locpath+locfile) :
        	field = Dataset(locpath+locfile)
        	LongTS_OBS_Ekm = npy.squeeze(field.variables['weMooringMonth'])
        	time_axis_Ekm = npy.squeeze(field.variables['time'])
        else:
        	LongTS_OBS_Ekm = npy.arange((2015-2003+1))+npy.nan
        	time_axis_Ekm = npy.arange((2015-2003+1))+npy.nan
        
        return  LongTS_OBS_Ekm, time_axis_Ekm
