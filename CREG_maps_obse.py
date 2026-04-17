#!/usr/bin/env python

import sys
import subprocess
import xarray as xr 
import numpy as npy
from datetime import datetime
from checkfile import *
import csv 

################################################################################################################################
def SSH_OBS( t_year=1959 ) :
################################################################################################################################

	locpath='./DATA/'
	locfile='EKE_DOT_based_2003-2014.nc'
	if chkfile(locpath+locfile,zstop=True) :
		ds_sshobs = xr.open_dataset(locpath+locfile)
		lon = ds_sshobs['lon']
		lat = ds_sshobs['lat']
		ssh_init = ds_sshobs['DOT']

		SSH_lon2D=npy.tile(lon,(lat.size,1))
		SSH_lat2D=npy.tile(lat,(lon.size,1)).T

	if t_year >= 2003 and t_year <=2014 :
		# Get the specific year
		s_ind=(t_year-2003)*12	 
		out_ssh_OBS = ssh_init.isel(date=slice(s_ind,s_ind+12)).mean(dim='date').squeeze()
		ssh_OBS_obsper = str(t_year)
	else:
		# Compute the mean over the obs. monthly period 2003-2014
		out_ssh_OBS = ssh_init.mean(dim='date').squeeze()
		ssh_OBS_obsper = '2003-2014'

	# Remove the domain mean to get an anomaly
	out_ssh_OBS = out_ssh_OBS - ssh_init.mean()

	return out_ssh_OBS, SSH_lon2D, SSH_lat2D, ssh_OBS_obsper

################################################################################################################################
def MLD_OBS() :
################################################################################################################################

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

	locpath='./DATA/'
	if zconfig == 'CREG025.L75' : 
		locfile='PIOMAS_icethic_interpCREG025.L75_1-12_1979-2020.nc'
	elif zconfig == 'CREG12.L75' :
		locfile='PIOMAS_icethic_interpCREG12.L75_1-12_1979-2020.nc'
	if chkfile(locpath+locfile,zstop=True) :
		ds_icet = xr.open_dataset(locpath+locfile)
		ICE_thick_init = ds_icet['icethic'].squeeze()

	if t_year >= 1979 and t_year <= 2020 :
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

	locpath='./DATA/'
	locfile='NSIDC-0051_92585_monthly.nc'
	if chkfile(locpath+locfile,zstop=True) :
		ds_icec = xr.open_dataset(locpath+locfile)
		lon = ds_icec['longitude'].squeeze()
		lat = ds_icec['latitude'].squeeze()
		CONC_init = ds_icec['Average_Sea_Ice_Concentration_with_Final_Version'].squeeze()

	## Initial data are based on add_offset
	## 251 > missing pole
	## 252 > not used
	## 253 > coastline
	## 254 > land
	## 255 > missing value
	## data will be recovered in dividing it by 250
	CONC_init = xr.where( CONC_init == 255, npy.nan, CONC_init )
	CONC_init = xr.where( CONC_init == 254, npy.nan, CONC_init )
	CONC_init = xr.where( CONC_init == 253, npy.nan, CONC_init )
	CONC_init = xr.where( CONC_init == 251, npy.nan, CONC_init )
	COR_CONC_init = CONC_init/250.

	CONC_init_land = npy.squeeze(CONC_init[0,:,:].copy())

	if t_year >= 1979 and t_year <= 2015 : 
		print( " Simplification calculation")
		# Select only March & September monthly mean
		mean_CONC_m03 = COR_CONC_init.sel(time=str(t_year)+'-03').squeeze()
		mean_CONC_m09 = COR_CONC_init.sel(time=str(t_year)+'-09').squeeze()
	else:
		# Compute a mean seasonal cycle and select March & September
		CONC_clim = COR_CONC_init.groupby('time.month').mean('time')
		mean_CONC_m03 = CONC_clim.isel(month=2)
		mean_CONC_m09 = CONC_clim.isel(month=8)

	mean_CONC_m03 = xr.where( CONC_init_land == 254 , npy.nan, mean_CONC_m03 )
	mean_CONC_m09 = xr.where( CONC_init_land == 254 , npy.nan, mean_CONC_m09 )

	return mean_CONC_m03, mean_CONC_m09, lon, lat


################################################################################################################################
def PHC3_OBS() :
################################################################################################################################

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

	print('				Read Von Appen et al. EKE inferred from Obs. ')
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
