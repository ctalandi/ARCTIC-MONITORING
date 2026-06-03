"""
CREG_maps_obse.py

Description:
This module defines a set of functions dedicated to read Obs. or reanalysis data set :
- EKE_OBS          : EKE deduced from DOT field. File: EKE_DOT_based_2003-2014.nc
- ICE_THICK_OBS    : ice thickness from PIOMAS reanalysis data set. File: PIOMAS_icethic_interpCREGXX.L75_1-12_1979-2024.nc
- ICE_CONCE_OBS    : September/March ice concentration from NSIDC v6 release. File: NSIDC-G02202-v6_sic_psn25_197811-202603_v06r00.nc
- MLD_OBS          : Mixed Layer Depth computed from MIMOC obs. File: MLD_MIMOC_based_monthlyClim_rhocrit0.01.nc
- MLTS_OBS         : mean temperature/salinity in the MLD infered from MIMOC obs. File: MIMOC_ML_v2.2_PT_S_MLP_Clim.nc
- PHC3_OBS         : temperature/salinity from PHC 3.0 climatology. File: phc3.0_annual.nc
- SSH_OBS          : SSH data from DOT obs. File: EKE_DOT_based_2003-2014.nc [2003-2014] and Full_DOT_data_Arco_2025_09.nc [2015-2025]
- VONAPPEN_EKE_OBS : EKE data infered from direct in-situ obs. File: EKE_table_Pangaea_lon_sorted_zero_nan_depth.txt

Author:
Claude Talandier (claude.talandier@cnrs.fr)
"""
import sys
import subprocess
import xarray as xr 
import numpy as npy
from datetime import datetime
import pandas as pd 
from checkfile import *
import csv 

################################################################################################################################
def SSH_OBS( t_year=1959 ) :
################################################################################################################################
	"""
	Function dedicated to read obs. data from PIOMAS (ice thickness) & NSIDC v6 (ice concentration) & IABP (ice drift)
	Respective filenames are :
		- PIOMAS_icevol_maskedBeringSea_interp+CONFIG+_1-12_1979-2024.nc [1979-2024]
		- NSIDC-G02202-v6_ice_area_and_extent_TiSe_y1978-11-2026-03_maskBeringSea_fullPoleGap.nc [1979-2025]
	
	Input:
	    t_year : current year to read 
	
	Output:
	    out_ssh_OBS    : SSH from Obs. 
	    SSH_lon2D      : 2D longitude 
	    SSH_lat2D      : 2D latitude 
	    ssh_OBS_obsper : string to mention the considered period in plots title 
	"""
	# ----------------------------------------------------------------------

	# DOT data set his based on the following satelite:
	# From 2003 - 2011 : Envisat 
	# From 2012 - 2025 : Cryosat-2 

	# Focus on the ENVISAT observation period 
	##########################################
	if t_year < 2003 : 
		locpath='./DATA/'
		locfile='EKE_DOT_based_2003-2014.nc'
		if chkfile(locpath+locfile,zstop=True) :
			ds_sshobs = xr.open_dataset(locpath+locfile)
			lon = ds_sshobs['lon']
			lat = ds_sshobs['lat']
			ssh_init = ds_sshobs['DOT']
	
			SSH_lon2D=npy.tile(lon,(lat.size,1))
			SSH_lat2D=npy.tile(lat,(lon.size,1)).T
	
		# Compute the mean over the obs. monthly period 2003-2014
		out_ssh_OBS = ssh_init.mean(dim='date').squeeze()
		ssh_OBS_obsper = 'Envisat 2003-2014'
	
		# Remove the domain mean to get an anomaly
		out_ssh_OBS = out_ssh_OBS - ssh_init.mean()

	elif t_year >= 2003 and t_year <=2014 :
		locpath='./DATA/'
		locfile='EKE_DOT_based_2003-2014.nc'
		if chkfile(locpath+locfile,zstop=True) :
			ds_sshobs = xr.open_dataset(locpath+locfile)
			lon = ds_sshobs['lon']
			lat = ds_sshobs['lat']
			ssh_init = ds_sshobs['DOT']
	
			SSH_lon2D=npy.tile(lon,(lat.size,1))
			SSH_lat2D=npy.tile(lat,(lon.size,1)).T

		# Get the specific year
		s_ind=(t_year-2003)*12	 
		out_ssh_OBS = ssh_init.isel(date=slice(s_ind,s_ind+12)).mean(dim='date').squeeze()
		ssh_OBS_obsper = 'Envisat '+str(t_year)

		# Remove the domain mean to get an anomaly
		out_ssh_OBS = out_ssh_OBS - ssh_init.mean()

	# Focus on the CRYOSAT-2 observation period 
	############################################
	elif t_year > 2014 and t_year <=2024 : 

		# Focus on the Cryosat-2 observation period 
		locpath='./DATA/'
		locfile='Full_DOT_data_Arco_2025_09.nc'
		if chkfile(locpath+locfile,zstop=True) :
			ds_sshobs = xr.open_dataset(locpath+locfile)
			SSH_lon2D = ds_sshobs['lons']
			SSH_lat2D = ds_sshobs['lats']
			ssh_init = ds_sshobs['DOT_smoothed']
			# Redefine a new time axis to manage data easily
			time_values = ds_sshobs.time.values
	
			# Build a DatetimeIndex starting January 1st 2000 and adding days
			datetime_index = pd.to_datetime('2000-01-01') + pd.to_timedelta(time_values, unit='D')
			# Associate this new time axis to the data	
			ssh_init['time'] = datetime_index
	
		# Remove all data before 2015 and after 2024 (2025 is from January to September)
		ssh_20152024 = ssh_init.where(ssh_init.time.dt.year >= 2015, drop=True)
		ssh_20152024 = ssh_20152024.where(ssh_20152024.time.dt.year <= 2024, drop=True)
		out_ssh_OBS = ssh_20152024.sel(time=str(t_year)).mean(dim='time') - ssh_20152024.mean()
		#out_ssh_OBS = out_ssh_OBS.squeeze()
		ssh_OBS_obsper = 'Cryosat-2 '+str(t_year)

	return out_ssh_OBS, SSH_lon2D, SSH_lat2D, ssh_OBS_obsper

################################################################################################################################
def MLD_OBS() :
################################################################################################################################
	"""
	Function dedicated to read MIMOC MLD infered from climatological temperature/salinity 
	Filename : MLD_MIMOC_based_monthlyClim_rhocrit0.01.nc
	
	Input:
	    None 
	
	Output:
	    mld_m03   : March MLD climatology
	    mld_m09   : September MLD climatology
	    MLD_lon2D : 2D longitude 
	    MLD_lat2D : 2D latitude 
	"""
	# ----------------------------------------------------------------------

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

################################################################################################################################
def MLTS_OBS() :
################################################################################################################################
	"""
	Function dedicated to read MIMOC MLD mean temperature/salinity infered from climatological temperature/salinity 
	Filename : MIMOC_ML_v2.2_PT_S_MLP_Clim.nc
	
	Input:
	    None 
	
	Output:
	    mlT_init : MLD mean temperature 
	    mlS_init : MLD mean salinity 
	    lon2D    : 2D longitude 
	    lat2D    : 2D latitude 
	"""
	# ----------------------------------------------------------------------

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


################################################################################################################################
def EKE_OBS( t_year=1959 ) :
################################################################################################################################
	"""
	Function dedicated to read EKE infered from DOT obs. data 
	Filename : EKE_DOT_based_2003-2014.nc
	
	Input:
	    t_year : current year to read 
	
	Output:
	    out_EKE_OBS : EKE infered from Obs. 
	    EKE_lon2D   : 2D longitude 
	    EKE_lat2D   : 2D latitude 
	"""
	# ----------------------------------------------------------------------

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

################################################################################################################################
def ICE_THICK_OBS( zconfig='CREG025.L75', t_year=1959 ) :
################################################################################################################################
	"""
	Function dedicated to read ice thickness from PIOMAS reanalysis 
	Filename : PIOMAS_icevol_maskedBeringSea_interp+CONFIG+_1-12_1979-2024.nc [1979-2024]
	
	Input:
	    zconfig : (optional) the considered model configuration (default = 'CREG025.L75') 
	    t_year  : current year to read 
	
	Output:
	    out_ICE_thick : ice thickness from observation
	"""
	# ----------------------------------------------------------------------

	locpath='./DATA/'
	locfile='PIOMAS_icethic_interp'+zconfig+'_1-12_1979-2024.nc'
	if chkfile(locpath+locfile,zstop=True) :
		ds_icet = xr.open_dataset(locpath+locfile)
		ICE_thick_init = ds_icet['icethic'].squeeze()

	if t_year >= 1979 and t_year <= 2024 :
		s_ind = (t_year-1979)*12
		mean_ICE_thick = ICE_thick_init[s_ind:s_ind+12,:,:].copy()
	else:
		mean_ICE_thick = ICE_thick_init.copy()
	
	# Annual or climatological mean
	out_ICE_thick = npy.mean(mean_ICE_thick,axis=0).squeeze()

	return out_ICE_thick


################################################################################################################################
def ICE_CONCE_OBS( t_year=1959 ) :
################################################################################################################################
	"""
	Function dedicated to read ice concentration obs. data NSIDC v6 
	Filename : NSIDC-G02202-v6_ice_area_and_extent_TiSe_y1978-11-2026-03_maskBeringSea_fullPoleGap.nc [1979-2025]
	
	Input:
	    t_year : current year to read 
	
	Output:
	    mean_CONC_m03 : March ice concentration of the current year 
	    mean_CONC_m09 : September ice concentration of the current year
	    lon           : 2D longitude  
	    lat           : 2D latitude
	"""
	# ----------------------------------------------------------------------

	locpath='./DATA/'
	locfile='NSIDC-G02202-v6_sic_psn25_197811-202603_v06r00.nc'
	if chkfile(locpath+locfile,zstop=True) :
		ds_icec = xr.open_dataset(locpath+locfile)
		ds_coor = xr.open_dataset(locpath+locfile, group='cdr_supplementary')
		ds_coor = ds_coor.set_coords('longitude')
		ds_coor = ds_coor.set_coords('latitude')
		lon = ds_coor['longitude'].squeeze()
		lat = ds_coor['latitude'].squeeze()
		CONC_init = ds_icec['cdr_seaice_conc_monthly'].squeeze()

	if t_year >= 1979 and t_year <= 2025 : 
		# Select only March & September monthly mean
		mean_CONC_m03 = CONC_init.sel(time=str(t_year)+'-03').squeeze()
		mean_CONC_m09 = CONC_init.sel(time=str(t_year)+'-09').squeeze()
	else:
		# Compute a mean seasonal cycle and select March & September
		CONC_clim = COR_CONC_init.groupby('time.month').mean('time')
		mean_CONC_m03 = CONC_clim.isel(month=2)
		mean_CONC_m09 = CONC_clim.isel(month=8)

	return mean_CONC_m03, mean_CONC_m09, lon, lat


################################################################################################################################
def PHC3_OBS() :
################################################################################################################################
	"""
	Function dedicated to read temperature/salinity from PHC3.0 climatology
	Filename : phc3.0_annual.nc
	
	Input:
	    None
	
	Output:
	    My_varTinit : 3D temperature climatology 
	    My_varSinit : 3D salinity climatology 
	    PHC_lon2D   : 2D longitude 
	    PHC_lat2D   : 2D latitude 
	"""
	# ----------------------------------------------------------------------

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


################################################################################################################################
def VONAPPEN_EKE_OBS() :
################################################################################################################################
	"""
	Function dedicated to read EKE from in-situ obs. at given moorings
	Filename : EKE_table_Pangaea_lon_sorted_zero_nan_depth.txt
	
	Input:
	    None
	
	Output:
	    dsVAD : EKE dataset at specific location and depths
	"""
	# ----------------------------------------------------------------------

	print('				Read Von Appen et al. EKE infered from Obs. ')
	locpath='./DATA/'
	locfile='EKE_table_Pangaea_lon_sorted_zero_nan_depth.txt'
	if chkfile(locpath+locfile) : 
		with open(locpath+locfile, newline='\n') as csvfile:
		    rd = csv.reader(csvfile, delimiter='\t')
		    VAD = npy.zeros((29, 212))
		    Names_VAD = ['' for k in range(212)]
		    i = 0
		    for row in rd:
		        if i>=1:
		            VAD[:, i-1] = row[1:]
		            Names_VAD[i-1] = row[0]
		        else:
		            Headers = row
		        i+=1
		    print(i)
		dsVAD = xr.Dataset(coords = {'mooring_loc':npy.arange(len(Names_VAD))})
		dsVAD['Names'] = (('moorings_loc'), Names_VAD)
		for i in range(1, len(Headers)):
		    dsVAD[Headers[i]] = (('moorings_loc'), VAD[i-1, :])

	return dsVAD
