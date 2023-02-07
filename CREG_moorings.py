#!/usr/bin/env python
#### Enthought Python Distribution (free for Academic) recommended

import matplotlib
matplotlib.use('Agg')
import sys
import PyRaf
import matplotlib.pylab as plt
import numpy as npy
from netCDF4 import Dataset
from checkfile import *
from CREG_moorings_cont import *
from CREG_moorings_func import *

s_year=XXSYEAXX
e_year=XXEYEAXX
lgTS_ys=XXLGTSSXX
lgTS_ye=XXLGTSEXX
xiosfreq=XXXIOSFREQXX
NCDF_OUT=XXNCDFOUTXX
m_alpha=1.
main_dir='./'
CONFIG='XXCONFXX'   ; CASE='XXCASEXX'     ;  CASE2='XXCASE2XX'
CONFCASE=CONFIG+'-'+CASE
data_dir=main_dir+CONFIG+'/'+CONFCASE+'-MEAN/'
grid_dir=main_dir+CONFIG+'/GRID/'

DIR_FIG_OUT='./'

read_IAinputdata=1  # read Interannual data (0/1)
monthly_data = 1  # Read monthly mean field for temporal evolution 
depth_time_plot=True   # Do the plot or not 
do_Kprofile=True     # Read a previous profile to compare the current one

########################################
# Read GRID 
########################################
#------------------------------------------------------------------------------------------------------------------------
# Read only once for the same section type
print
print ''
print '    >>>>>>>>>>>>>>>>>>>>>>>>>>>>Perform the treatment of the '+CONFCASE+' experiment '
print ''
print
print '		Start to read geographical, scale factors & masks '
print

locpath=grid_dir
locfile=CONFCASE+'_mask.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	fieldmask=Dataset(locpath+locfile)
	tmask = npy.squeeze(fieldmask.variables['tmask'])

locpath=grid_dir
locfile=CONFCASE+'_mesh_hgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	fieldhmesh=Dataset(locpath+locfile)
	lon = fieldhmesh.variables['nav_lon']
	lat = fieldhmesh.variables['nav_lat']
	e1t= npy.squeeze(fieldhmesh.variables['e1t'])
	e2t= npy.squeeze(fieldhmesh.variables['e2t'])

locpath=grid_dir
locfile=CONFCASE+'_mesh_zgr.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	fieldzmesh=Dataset(locpath+locfile)
	e3t = npy.squeeze(fieldzmesh.variables['e3t_0'])

name_z='deptht'

locpath=data_dir+'/'+str(s_year)+'/'
locfile=CONFCASE+'_y'+str(s_year)+'.'+xiosfreq+'_gridT.nc'
if chkfile(locpath+locfile,zstop=True,zscript=sys.argv[0]) :
	z = PyRaf.readfullNC4(locpath+locfile,name_z)
z2dt=npy.reshape(z,(z.size,1))
print '		End to read geographical, scale factors & masks '

# Define time axis
if monthly_data == 1:
        t_months=(npy.arange(12)*30.+15.)/365.
        c_year=s_year   ; start = 1
        while c_year <= e_year:
                y_years=npy.tile(c_year,12)+t_months
                if start == 1:
                        time_axis=y_years
                        start=0
                else:
                        time_axis=npy.append(time_axis,y_years)
                c_year+=1
else:
        time_axis=npy.arange(e_year-s_year+1)+s_year+0.5

if monthly_data == 1 and e_year-s_year == 0:
	time_axis=npy.arange(12)+1.5

xplt= npy.tile(time_axis,(z.size,1))
zplt = npy.repeat(z2dt,time_axis.shape[0],axis=1)

e1t3D = npy.tile(e1t,(e3t.shape[0],1,1))
e2t3D = npy.tile(e2t,(e3t.shape[0],1,1))
del e1t
del e2t


if monthly_data == 1:
        time_dim=(e_year-s_year+1)*12
else:
        time_dim=e_year-s_year+1

if time_dim == 0:
	time_dim = 1

nbvar = 2 
Ar_size=(nbvar,time_dim,e3t.shape[0],lon.shape[0],lon.shape[1])
My_var=npy.empty(Ar_size)
Ar_size=(nbvar,e3t.shape[0],lon.shape[0],lon.shape[1])
My_varinit=npy.empty(Ar_size)

if read_IAinputdata != 1:
	Ar_size=(time_dim,e3t.shape[0],lon.shape[0],lon.shape[1])
	My_varTemp=npy.empty(Ar_size)
#del lat
#del lon


nvar=0
for var in DEF_INFO_VAR('gridT','votemper'):
	t_month=0

	if read_IAinputdata == 1:
		print 
		print 'Start to read ', var['name'],' data '
		print 
		# Read the whole time series
		c_year=s_year
		while c_year <= e_year:
               		if monthly_data == 1:
               		    c_month=0
               		    while c_month <= 11:
				smonth='m0'+str(c_month+1) if c_month <= 8 else 'm'+str(c_month+1)
        			locpath=data_dir+'/'+str(c_year)+'/'
        			locfile=CONFCASE+'_y'+str(c_year)+smonth+'.'+xiosfreq+'_gridT.nc'
				if chkfile(locpath+locfile) :
				 	My_var[nvar,t_month,:,:,:] = PyRaf.readfullNC4(locpath+locfile,var['name'])
                                t_month+=1
                                c_month+=1
				print 'Read year :', c_year, ' month:', c_month
			else:
        		    locpath=data_dir+'/'+str(c_year)+'/'
        		    locfile=CONFCASE+'_y'+str(c_year)+'.'+xiosfreq+'_gridT.nc'
			    if chkfile(locpath+locfile) :
			    	My_var[nvar,c_year-s_year,:,:,:] = PyRaf.readfullNC4(locpath+locfile,var['name'])
			    print 'Read year :', c_year

			c_year = c_year + 1 
	else:
        	locpath='./'
        	locfile=CONFCASE+'_gridT_'+var['shortname']+'_IA.npy'
		if chkfile(locpath+locfile) :
			print 'Reload data from:', locpath+locfile
			My_varTemp[:,:,:,:]=npy.load(locpath+locfile,mmap_mode='r')
                	My_var[nvar,:,:,:,:]=My_varTemp[:,:,:,:].copy()
			print 'Check shape of reloaded data array:', My_var.shape
			#print 'Reloaded data max value:', My_var.max()
			#print 'Reloaded data min value:', My_var.min()

	# Read the initial state 
        locpath=data_dir+'/'
        locfile=CONFCASE+'_init_gridT.nc'
	if chkfile(locpath+locfile) :
		My_varinit[0,:,:,:] = PyRaf.readfullNC4(locpath+locfile,'votemper')
		My_varinit[1,:,:,:] = PyRaf.readfullNC4(locpath+locfile,'vosaline')

	nvar+=1

####################################################################################################################
####################################################################################################################
###########################################   LONG TIME SERIES PROFILES   ##########################################
####################################################################################################################
####################################################################################################################
#if False: 
if lgTS_ye-lgTS_ys+1 > 1 :
        plt.clf()
	plt.subplots_adjust(bottom=0.1, left=0.1, right=0.95, top=0.9)
	num_fram=410
	All_box=XXALLBOXXX
	All_var=DEF_INFO_VAR('gridT','votemper')
	cumul_box=''
	moor_n=1
        for zmybox in All_box :
		zbox=DEF_MOOR_BOX(CONFIG,zmybox)
		hsct_lev,Red_My_varinit=DEF_REDVAR(CONFIG,CASE,My_var,My_varinit,zbox,nbvar,time_dim,e3t,All_var,s_year,e_year,tmask,lon,lat,e1t3D,e2t3D)
		fram=num_fram+moor_n
		DEF_ZTIME(CONFIG,CASE,lgTS_ys,lgTS_ye,zbox,z,z2dt,hsct_lev,zfram=fram,zoutNC=NCDF_OUT) 
		cumul_box=cumul_box+'-'+zbox['box']
		moor_n+=2

	# Save year time series to use for a longer time series
	timeadd='_month' if monthly_data == 1 else ''

	file_ext1='_TemSalEvol'+timeadd
	file_ext=file_ext1+cumul_box+'_LGTS_y'+str(lgTS_ys)+'LASTy'
	file_name=CONFCASE+file_ext

	plt.tight_layout()
        plt.savefig(DIR_FIG_OUT+file_name+'.pdf')


####################################################################################################################
####################################################################################################################
###########################################   VERTICAL PROFILES   ##################################################
####################################################################################################################
####################################################################################################################
if do_Kprofile :
	print
	print " 			##################################################################  " 
	print " 			##################################################################  " 
	print " 			######### PLOT VERTICAL PROFILES AT MOORINGS LOCATION ############  " 
	print " 			##################################################################  " 
	print " 			##################################################################  " 
	print
        plt.clf()
	num_fram=220
	All_box=XXALLBOXXX
	All_var=DEF_INFO_VAR('gridT','votemper')
	cumul_box=''
	moor_n=1
        for zmybox in All_box :
		zbox=DEF_MOOR_BOX(CONFIG,zmybox)
		fram=num_fram+moor_n
		hsct_lev,Red_My_varinit=DEF_REDVAR(CONFIG,CASE,My_var,My_varinit,zbox,nbvar,time_dim,e3t,All_var,s_year,e_year,tmask,lon,lat,e1t3D,e2t3D)
		DEF_ZPROFILE(CONFIG,CASE,str(s_year),hsct_lev,Red_My_varinit,zbox,zplt,zfram=fram,zs_year=s_year) 
		cumul_box=cumul_box+'-'+zbox['box']
		moor_n+=2

        pfile_ext='_TemSaliProfile'+cumul_box+'_y'+str(s_year)+'_z0-500'
        file_name=CONFIG+'-'+CASE+pfile_ext
	plt.tight_layout()
        plt.savefig(DIR_FIG_OUT+file_name+'.pdf')







