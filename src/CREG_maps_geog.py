#!/usr/bin/env python
"""
CREG_maps_geog.py

Description:
This module is dedicated to plot the moorings and sections location available in the current monitoring

Author:
Claude Talandier (claude.talandier@cnrs.fr)
"""

import sys 
import matplotlib
matplotlib.use('Agg')
import numpy as npy
from CREG_maps_func import *
from checkfile import *
import matplotlib.pylab as plt
import matplotlib as mpl
from cartopy import crs as ccrs
import xarray as xr

main_dir='./'
CONFIG='XXCONFXX'   ; CASE='XXCASEXX'	  
CONFCASE=CONFIG+'-'+CASE
grid_dir=main_dir+CONFIG+'/GRID/'

print()
print('				       Configuration :' + CONFCASE)
print()

########################################
# Read GRID 
########################################
#------------------------------------------------------------------------------------------------------------------------
# Read only once for the same section type
locpath=grid_dir
locfile=CONFCASE+'_mask.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	ds_msk = xr.open_dataset(locpath+locfile)[['glamt','gphit']]
	lon = ds_msk['glamt'].squeeze()
	lat = ds_msk['gphit'].squeeze()
#------------------------------------------------------------------------------------------------------------------------

	fig=plt.figure() ;  fram=111
	projection = ccrs.NorthPolarStereo(central_longitude=-60, true_scale_latitude=65)
	ax = fig.add_subplot(fram, projection=projection)

	BATHY_MAP( ztype='isol500', ax=ax )

	ax.set_aspect('equal')  
	ax.set_extent([-180, 180, 65, 90], crs=ccrs.PlateCarree())
	
	ax.add_feature(cartopy.feature.LAND, facecolor='dimgray')
	gridlines = ax.gridlines(draw_labels=True, linestyle=':', linewidth=0.4, alpha=0.7, xlocs=npy.arange(-180, 181, 20), y_inline=True, rotate_labels=False )
	gridlines.xlabel_style={'fontsize': 4}
	gridlines.ylabel_style={'fontsize': 4}
	gridlines.top_labels=False
	gridlines.left_labels=True
	gridlines.right_labels=True

	############################################################################################################
	############################################################################################################
	boxtoplot=1
	if boxtoplot == 1:
		tmskBFG = npy.ones((lon.shape[0],lon.shape[1]))
		tmskBFG = npy.ma.masked_where( lat[:,:] >  80.5,tmskBFG )
		tmskBFG = npy.ma.masked_where( lat[:,:] <  70.5,tmskBFG )
		tmskBFG = npy.ma.masked_where( lon[:,:] > -130.,tmskBFG )
		tmskBFG = npy.ma.masked_where( lon[:,:] < -170.,tmskBFG )

		norm = mpl.colors.Normalize(vmin=0., vmax=1.)
		pal = plt.get_cmap('cool')
		ax.contourf( lon,lat,tmskBFG,levels=[0.,1.],cmap=pal,alpha=0.4, transform=ccrs.PlateCarree() )
		props = dict(boxstyle='round', facecolor='w', alpha=1.0)
		fig.text(0.33, 0.61, 'CRF-Box', color='b',fontsize=7, bbox=props, alpha=0.5)
	############################################################################################################
	############################################################################################################

	#fig.title(" MONARC montioring moorings, sections & boxes location")
	############################################################################################################
	############################################################################################################
	moorplot=1
	if moorplot == 1 :
		bx_ARCB={'name':'B'  ,'lon_min':-150.,'lon_max':-150.,'lat_min':78.,'lat_max':78.}
		bx_EURA={'name':'EUR','lon_min':  60.,'lon_max':  60.,'lat_min':85.,'lat_max':85.}

		All_box=[bx_ARCB,bx_EURA]
		for box in All_box:
			lats = [box['lat_min'],box['lat_max']]
			lons = [box['lon_min'],box['lon_max']]
			ax.scatter(lons,lats,1,marker='o', color='r', transform=ccrs.PlateCarree())

	props = dict(boxstyle='round', facecolor='w', alpha=1.0)
	fig.text(0.37, 0.50, 'ARC-B', color='r',fontsize=7, bbox=props)
	fig.text(0.53, 0.57, 'EURA' , color='r',fontsize=7, bbox=props)
	############################################################################################################
	############################################################################################################
	############################################################################################################
	############################################################################################################
	secplot=1
	if secplot == 1 :
		sec_BEAU={'name':"Beaufort" ,'jmin': 395,'jmax':560 ,'imin':197,'labx':0.31 ,'laby':0.37}
		sec_ANNA={'name':"St Anna"  ,'imin': 334,'imax':358 ,'jmax':395,'labx':0.62 ,'laby':0.54}
		sec_KARA={'name':"Kara"    ,'imin': 291,'imax':362 ,'jmax':447,'labx':0.58 ,'laby':0.66}
		sec_FRAM={'name':"Fram"    ,'imin': 301,'imax':333 ,'jmax':338,'labx':0.63 ,'laby':0.40}
		sec_BERI={'name':"Bering"   ,'imin': 202,'imax':223 ,'jmax':601,'labx':0.20 ,'laby':0.63}
		props = dict(boxstyle='round', facecolor='w', alpha=1.0)

		All_sec=[sec_KARA,sec_FRAM,sec_BERI]
		for box in All_sec:
			#################
			# ZONAL SECTIONS
			#################
			ji=box['imin']
			while ji <= box['imax']-1 :
				lats = [lat[box['jmax'],ji], lat[box['jmax'],ji+1]]
				lons = [lon[box['jmax'],ji], lon[box['jmax'],ji+1]]
				ax.plot(lons,lats, color='g', transform=ccrs.PlateCarree())
				ji+=1
			fig.text(box['labx'],box['laby'] , box['name'], color='g',fontsize=7, bbox=props)

		All_sec=[sec_BEAU]
		for box in All_sec:
			#################
			# MERIDIONAL SECTIONS
			#################
			jj=box['jmin']
			while jj <= box['jmax']-1 :
				lats = [lat[jj,box['imin']], lat[jj+1,box['imin']]]
				lons = [lon[jj,box['imin']], lon[jj+1,box['imin']]]
				ax.plot(lons,lats, color='g', transform=ccrs.PlateCarree())
				jj+=1
			fig.text(box['labx'],box['laby'] , box['name'], color='g',fontsize=7, bbox=props)

	zfile_ext='_ARC-GEO_'
	plt.tight_layout()
	plt.savefig('MONARC_ARC-GEOLOC.png', dpi=300)
