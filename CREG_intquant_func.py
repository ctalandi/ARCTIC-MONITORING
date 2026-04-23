import numpy as npy
import scipy.io as sio
from checkfile import *
import xarray as xr
import pandas as pd 
from fsspec.implementations.local import LocalFileSystem
fs = LocalFileSystem()

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
        locpath='./DATA/'
        locfile = 'PIOMAS_icevol_maskedBeringSea_interp'+CONFIG+'_1-12_1979-2024.nc'
        if chkfile(locpath+locfile) :
                ds_fld = xr.open_dataset(locpath+locfile)
                LongTS_OBS_icevol = ds_fld['icevol-BS']
        else:
                LongTS_OBS_icevol = npy.arange((2024-1979+1)*12)+npy.nan
        # Set the time axis for PIOMAS observations
        lgts_year=1979    ;      start = 1
        while  lgts_year <= 2024  :
                y_years=npy.tile(lgts_year,12)+t_months
                if start == 1:
                        time_axis_PIO=y_years
                        start=0
                else:
                        time_axis_PIO=npy.append(time_axis_PIO,y_years)
                lgts_year+=1

        # Ice extent obs
        # Data start in November 1978, so there is a 2 months shift
        locpath='./DATA/'
        locfile='NSIDC-G02202-V4_ice_area_and_extent_TiSe_y1978-11-2022-12_maskBeringSea_fullPoleGap.nc'
        if chkfile(locpath+locfile) :
                ds_fld = xr.open_dataset(locpath+locfile)
                LongTS_OBS_iceext = ds_fld['ice_extent'][2::]
                LongTS_OBS_iceare = ds_fld['ice_area'][2::]
        else:
                LongTS_OBS_iceext = npy.arange((2022-1979+1)*12)+npy.nan
                LongTS_OBS_iceare = npy.arange((2022-1979+1)*12)+npy.nan
	# Set the time axis for NSIDC observations
        lgts_year=1979    ;      start = 1
        while  lgts_year <= 2022  :
               y_years=npy.tile(lgts_year,12)+t_months
               if start == 1:
                       time_axis_NSIDC=y_years
                       start=0
               else:
                       time_axis_NSIDC=npy.append(time_axis_NSIDC,y_years)
               lgts_year+=1

        # Return also the September ice extent 
        LongTS_OBS_Septiceext = ds_fld['ice_extent'].sel(time=ds_fld.time.dt.month == 9)
        LongTS_OBS_Septiceext['time'] = pd.date_range(start='1979-01',end='2022-12',freq='YS') + pd.DateOffset(days=180)

        # Ice drift from IABP
        # Data start in 18/01/1979, with 784 Buoys and 2 smapling / day : 0 & 12 
        locpath='./DATA/'
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
        
        return LongTS_OBS_icevol, LongTS_OBS_iceext, LongTS_OBS_Septiceext, LongTS_OBS_iceare, IABPObservations, time_axis_obs, time_axis_PIO, time_axis_NSIDC

def CAL_ICE_VOL_AREA(CONFIG,CASE,lgts_year,data_dir,xiosfreq,dom_area,tmask2D) :

	# List files to be read
	locpath=data_dir+'/'+str(lgts_year)+'/'+xiosfreq+'/'
	locfile=CONFIG+'-'+CASE+'_y'+str(lgts_year)+'m??.'+xiosfreq+'_icemod.nc'
	ICE_files = [f for f in fs.glob(locpath+locfile)]
	
	drp_var=["time_centered_bounds","time_counter_bounds","siages"]
	if len(ICE_files) == 12 :
	   ds_Idata = xr.open_mfdataset(locpath+locfile, engine="netcdf4", concat_dim=["time_counter"], combine='nested', parallel=True, drop_variables=drp_var)
	   Ithdata_read = ds_Idata['sivolu']
	   Ifrdata_read = ds_Idata['siconc']

	   # Fill NaN with zero
	   Ithdata_read= Ithdata_read.fillna(0)
	   Ifrdata_read= Ifrdata_read.fillna(0)

	ice_volume = npy.sum( Ithdata_read * dom_area * tmask2D, axis=(1,2) )
	ice_area   = npy.sum( Ifrdata_read * dom_area * tmask2D, axis=(1,2) )
        
	return ice_volume, ice_area

def READ_OBS_LGTS_CRFFWC(lgTS_ys,lgTS_ye) :
        # Function dedicated to read observations data that don't change 

        # CRF box mean FWC Obs.
        locpath='./DATA/'
        locfile = 'BGFWC_OI_2023.nc'
        if chkfile(locpath+locfile) :
                ds_fld = xr.open_dataset(locpath+locfile)
                LongTS_OBS_FWC = ds_fld['fwc_total']
                time_axis_FWC = ds_fld['year']
        else:
                LongTS_OBS_FWC = npy.arange((2023-2003+1))+npy.nan
                time_axis_FWC = npy.arange((2023-2003+1))+npy.nan

        return  LongTS_OBS_FWC, time_axis_FWC


def READ_OBS_LGTS_CRFEkm(lgTS_ys,lgTS_ye) :
        # Function dedicated to read observations data that don't change 
        
        # CRF box mean Ekman pumping Obs.
        locpath='./DATA/'
        locfile = 'ArcticEkmanPumping_MonthlyMean.nc'
        if chkfile(locpath+locfile) :
                ds_fld = xr.open_dataset(locpath+locfile)
                LongTS_OBS_Ekm = ds_fld['weMooringMonth'].squeeze()
                time_axis_Ekm = ds_fld['time'].squeeze()
        else:
                LongTS_OBS_Ekm = npy.arange((2015-2003+1))+npy.nan
                time_axis_Ekm = npy.arange((2015-2003+1))+npy.nan
        
        return  LongTS_OBS_Ekm, time_axis_Ekm
