#!/usr/bin/env python

import sys 
import matplotlib
matplotlib.use('Agg')
import numpy as npy
import CREG_maps_func
from checkfile import *
from netCDF4 import Dataset
import matplotlib.pylab as plt
import matplotlib as mpl

main_dir='./'
CONFIG='XXCONFXX'   ; CASE='XXCASEXX'     
CONFCASE=CONFIG+'-'+CASE
grid_dir=main_dir+CONFIG+'/GRID/'

print
print '                              Configuration :' , CONFCASE
print


########################################
# Read GRID 
########################################
#------------------------------------------------------------------------------------------------------------------------
# Read only once for the same section type
locpath=grid_dir
locfile=CONFCASE+'_mask.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	fieldmask=Dataset(locpath+locfile)
	lon = npy.array(fieldmask.variables['nav_lon'])
	lat = npy.array(fieldmask.variables['nav_lat'])
#------------------------------------------------------------------------------------------------------------------------

	fig=plt.figure()

	zoutmap,X,Y=CREG_maps_func.Arc_Bat(ztype='isomonarc')

        ############################################################################################################
	############################################################################################################
	boxtoplot=1
	if boxtoplot == 1:
		tmskBFG=npy.ones((lon.shape[0],lon.shape[1]))
		tmskBFG=npy.ma.masked_where(lat[:,:] >  80.5,tmskBFG)
		tmskBFG=npy.ma.masked_where(lat[:,:] <  70.5,tmskBFG)
		tmskBFG=npy.ma.masked_where(lon[:,:] > -130.,tmskBFG)
		tmskBFG=npy.ma.masked_where(lon[:,:] < -170.,tmskBFG)

		norm = mpl.colors.Normalize(vmin=0., vmax=1.)
		pal = plt.cm.get_cmap('cool')
		#pal = plt.cm.get_cmap('binary')
		C2= zoutmap.contourf(X,Y,tmskBFG,[0.,1.],cmap=pal,norm=norm,alpha=0.4)
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
        		x,y = zoutmap(lons,lats)
        		zoutmap.scatter(x,y,5,marker='o', color='r')

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
                	        x,y = zoutmap(lons,lats)
                	        zoutmap.plot(x,y,linewidth=1,color='g')
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
                	        x,y = zoutmap(lons,lats)
                	        zoutmap.plot(x,y,linewidth=1,color='g')
                	        jj+=1
			fig.text(box['labx'],box['laby'] , box['name'], color='g',fontsize=7, bbox=props)


	zfontsize=8.
	zoutmap.drawparallels(npy.arange(-90.,91.,5.),labels=[False,False,False,False], color='grey', size=zfontsize, linewidth=0.3, alpha=0.5)
	zoutmap.drawmeridians(npy.arange(-180.,181.,20.),labels=[True,True,False,True], color='grey', size=zfontsize, latmax=90.,linewidth=0.3, alpha=0.5)
	zoutmap.fillcontinents(color='grey',lake_color='white',alpha=0.8)


	zfile_ext='_ARC-GEO_'
        plt.tight_layout()
	plt.savefig('MONARC_ARC-GEOLOC.pdf')
