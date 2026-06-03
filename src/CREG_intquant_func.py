"""
CREG_intquant_func.py

Description:
This module defines a set of functions dedicated to :
- CAL_ICE_VOL_AREA     : compute model ice volume & area
- CONVERT_SA2PS        : convert salinity from Absolute to Practical units based on GSW package
- READ_OBS_LGTS_CRFEkm : read obs. Ekman pumping time-series computed in the Beaufort Gyre
- READ_OBS_LGTS_CRFFWC : read obs. Fresh Water Content time-series computed in the Beaufort Gyre
- READ_MOD_LGTS_DATA   : read model data time-series in a numpy npy file format (Default)
- READ_OBS_LGTS_DATA   : read obs. data from PIOMAS (ice thickness) & NSIDC v6 (ice concentration) & IABP (ice drift)
- P4Dtzyx              : transform a 1D vertical pressure into 4D field with 12 time records

Author:
Claude Talandier (claude.talandier@cnrs.fr)
"""
import numpy as npy
import scipy.io as sio
from checkfile import *
import xarray as xr
import pandas as pd 
import gsw as gsw 
from fsspec.implementations.local import LocalFileSystem
fs = LocalFileSystem()

################################################################################################################################
def READ_MOD_LGTS_DATA(locpath,locfile,zLongTS,npz=None,zvarname=None) :
################################################################################################################################
	"""
	Function dedicated to read model data time-series in a numpy npy file format (Default)
	
	Input:
	    locpath  : the full path directory where is stored the file name to read
	    locfile  : the file name to read
	    zLongTS  : a 1D array corresponding to the time-series to append
	    npz      : (optional) format of data files (default is npy format)
	    zvarname : (optional) the name of the variable corresponding to years 
	
	Output:
	    out : A long time-series for a given diagnostic
	"""
	# ----------------------------------------------------------------------
	# Function dedicated to read model data in numpy npy format (Default)
	if chkfile(locpath+locfile) :
		year_fld=npy.load(locpath+locfile,mmap_mode='r')
		zLongTS=npy.append(zLongTS,year_fld) if npz == None else npy.append(zLongTS,year_fld[zvarname])
	else:
		zLongTS=npy.append(zLongTS,npy.arange(12)+npy.nan)

	return  zLongTS


################################################################################################################################
def READ_OBS_LGTS_DATA(CONFIG,lgTS_ys,lgTS_ye) :
################################################################################################################################
        """
        Function dedicated to read obs. data from PIOMAS (ice thickness) & NSIDC v6 (ice concentration) & IABP (ice drift)
        Respective filenames are :
        	- PIOMAS_icevol_maskedBeringSea_interp+CONFIG+_1-12_1979-2024.nc (over the range 1979-2024)
        	- NSIDC-G02202-v6_ice_area_and_extent_TiSe_y1978-11-2026-03_maskBeringSea_fullPoleGap.nc (over the range 1979-2025)
        	- IABP_ice_drift_BG_1979-2016.mat (over the range 1979-2016)
        
        Compute also the ice extent and its September minimum as well
        
        Input:
            CONFIG  : the configuration name to know the grid on which ice thicness has been interpolated 
            lgTS_ys : the first year to read 
            lgTS_ye : the last year to read 
        
        Output:
            LongTS_OBS_icevol          : PIOMAS ice thickness time-series 
            LongTS_OBS_iceext          : NSIDC v6 ice extent time-series 
            LongTS_OBS_Septiceext      : NSIDC v6 september ice extent time-series 
            LongTS_OBS_iceare          : NSIDC v6 ice area time-series
            IABPObservations           : IABP ice drift time-series 
            time_axis_IABP             : IABP data time axis 
            time_axis_PIO              : PIOMAS data time axis 
            time_axis_NSIDC            : NSIDC data time axis 
            LongTS_OBS_Septiceext_time : NSIDC data time axis for September minimum ice extent 
        """
        # ----------------------------------------------------------------------
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
        # NSIDC-v6 starts in November 1978 ending in March 2026 so skip the first 2 months and last 3 months of the time-series
        locpath='./DATA/'
        locfile='NSIDC-G02202-v6_ice_area_and_extent_TiSe_y1978-11-2026-03_maskBeringSea_fullPoleGap.nc'
        if chkfile(locpath+locfile) :
                ds_fld = xr.open_dataset(locpath+locfile)
                LongTS_OBS_iceext = ds_fld['ice_extent'].sel(time=slice('1979','2025'))
                LongTS_OBS_iceare = ds_fld['ice_area'].sel(time=slice('1979','2025'))
        else:
                LongTS_OBS_iceext = npy.arange((2025-1979+1)*12)+npy.nan
                LongTS_OBS_iceare = npy.arange((2025-1979+1)*12)+npy.nan
	# Set the time axis for NSIDC observations
        lgts_year=1979    ;      start = 1
        while  lgts_year <= 2025  :
               y_years=npy.tile(lgts_year,12)+t_months
               if start == 1:
                       time_axis_NSIDC=y_years
                       start=0
               else:
                       time_axis_NSIDC=npy.append(time_axis_NSIDC,y_years)
               lgts_year+=1

        # Return also the September ice extent 
        LongTS_OBS_Septiceext = ds_fld['ice_extent'].sel(time=ds_fld.time.dt.month == 9)
        LongTS_OBS_Septiceext_time = pd.date_range(start='1979-01',end='2025-12',freq='YS') + pd.DateOffset(days=257)

        # Ice drift from IABP
        # Data start in 18/01/1979, with 784 Buoys and 2 smapling / day : 0 & 12 
        locpath='./DATA/'
        locfile='IABP_ice_drift_BG_1979-2016.mat'
        if chkfile(locpath+locfile) :
                IABPObservations_read = sio.loadmat(locpath+locfile,squeeze_me=True)
                IABPObservations = npy.array(IABPObservations_read['time_series'])
        else:
                IABPObservations = npy.arange(456)+npy.nan
        
        # Set the time axis for observations
        lgts_year=1979    ;      start = 1
        while  lgts_year <= 2016  :
                y_years=npy.tile(lgts_year,12)+t_months
                if start == 1:
                        time_axis_IABP=y_years
                        start=0
                else:
                        time_axis_IABP=npy.append(time_axis_IABP,y_years)
                lgts_year+=1
        
        return LongTS_OBS_icevol, LongTS_OBS_iceext, LongTS_OBS_Septiceext, LongTS_OBS_iceare, IABPObservations, time_axis_IABP, time_axis_PIO, time_axis_NSIDC, LongTS_OBS_Septiceext_time

################################################################################################################################
def CAL_ICE_VOL_AREA(CONFIG,CASE,lgts_year,data_dir,xiosfreq,dom_area,tmask2D) :
################################################################################################################################
	"""
	Function dedicated to compute model ice volume & area
	
	Input:
	    CONFIG    : the configuration name to know the grid on which ice thicness has been interpolated 
	    CASE      : the experiment name associated to the CONFIG 
	    lgts_year : the curent year to read 
	    data_dir  : the full path of model data 
	    xiosfreq  : the output frequency used in the model filename 
	    dom_area  : surface of each model grid cell
	    tmask2D   : ocean/land model mask at the surface (T-point)
	
	Output:
	    ice_volume : model ice volume time series 
            ice_area   : model ice area time series 
	"""
	# ----------------------------------------------------------------------

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

################################################################################################################################
def READ_OBS_LGTS_CRFFWC(lgTS_ys,lgTS_ye) :
################################################################################################################################
        """
        Function dedicated to read obs. Fresh Water Content computed in the Beaufort Gyre 
        Filename is : - BGFWC_OI_2023.nc (over the range 2003-2023)
        
        Input:
            lgTS_ys : the first year to read 
            lgTS_ye : the last year to read 
        
        Output:
            LongTS_OBS_FWC : FWC obs. time-series 
            time_axis_FWC  : Obs. time axis 
        """
        # ----------------------------------------------------------------------

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


################################################################################################################################
def READ_OBS_LGTS_CRFEkm(lgTS_ys,lgTS_ye) :
################################################################################################################################
        """
        Function dedicated to read obs. Ekman pumping time-series computed in the Beaufort Gyre 
        Filename is : - ArcticEkmanPumping_MonthlyMean.nc (over the range 2003-2015)
        
        Input:
            lgTS_ys : the first year to read 
            lgTS_ye : the last year to read 
        
        Output:
            LongTS_OBS_Ekm : Obs. Ekman pumping time-series 
            time_axis_Ekm  : Obs. data time axis 
        """
        # ----------------------------------------------------------------------
        
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

################################################################################################################################
def CONVERT_SA2PS( zSAL, zlon, zlat ) :
################################################################################################################################
        """
        Function dedicated to convert salinity from Absolute to Practical based on GSW package
        
        Input:
            zSAL : Absolute Salinity 3D field 
            zlon : longitude 2D field 
            zlat : latitude 2D field 
        
        Output:
            z_SP : Practical Salinity 3D field 
        """
        # ----------------------------------------------------------------------

        # Compute the pressure at each depth
        pressure = gsw.p_from_z( -zSAL.z.values.squeeze(), 77. )
        pressure4D = P4Dtzyx( pressure, zlat ) 

        # Apply the conversion
        z_SP = gsw.conversions.SP_from_SA( zSAL, pressure4D, zlon, zlat )

        return z_SP.astype('float32')

################################################################################################################################
def P4Dtzyx( zvector, zlat ) :
################################################################################################################################
        """
        Function dedicated to transform a 1D vertical pressure into 4D field with 12 time records 
        
        Input:
            zvector : Vertical pressure vector 
            zlat    : latitude 2D field 
        
        Output:
            zpressure4D : presure 4D field 
        """
        # ----------------------------------------------------------------------

        # Prepare this 1D field to be duplicated in 4D : time,z,y,x
        z2dt = npy.reshape( zvector, (1, len(zvector), 1, 1) )

        # Time axis 
        zplt = npy.repeat( z2dt, 12, axis=0 )

        # Horizontaly to fit the T/S on a global grid
        zpressure4D = npy.tile( zplt,( 1, 1, zlat.shape[0], zlat.shape[1] ) )

        return zpressure4D
