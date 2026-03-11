
import numpy as npy
import matplotlib.pylab as plt 
import scipy.io as sio
from checkfile import *
from CREG_moorings_cont import *
import subprocess
from datetime import datetime
import xarray as xr 
from fsspec.implementations.local import LocalFileSystem
fs = LocalFileSystem()

def DEF_INFO_VAR(zgtype,zMyvar):
	# The following box indices have been defined for the ARCTIC area
	D_temp={'name':'votemper','units':'DegC','longname':'Temperature','shortname':'T','minval':3,'maxval':5}
	D_sali={'name':'vosaline','units':'PSU' ,'longname':'Salinity','shortname':'S','minval':34.5,'maxval':35.5}
	D_sig0={'name':'rhop_sig0','units':'kg m-3' ,'longname':'potential density','shortname':'sig0','minval':26.5,'maxval':28.0}
	D_Kz={'name':'votkeavt','units':'m2 s-1' ,'longname':'Vertical diffusivity','shortname':'Kz','minval':-7.5,'maxval':.0}
	D_W={'name':'vovecrtz','units':'m s-1' ,'longname':'Vertical velocity','shortname':'W','minval':-1.,'maxval':1.}
	D_Umod={'name':'speed','units':'m s-1' ,'longname':'Velocity module','shortname':'Vel','minval':0.,'maxval':5.}
	
	if zgtype == "gridT" :
		if zMyvar == 'rhop_sig0' or zMyvar == 'vosigma0' :
			zAll_var=[D_sig0]
		else:
			zAll_var=[D_temp,D_sali]
	elif zgtype == "gridU":
		zAll_var=[D_Umod]
	elif zgtype == "gridW":
		zAll_var=[D_Kz,D_W]

	return zAll_var


def DEF_MOOR_BOX(CONFIG,box2sel):
	# The following box indices have been defined for the ARCTIC area

	if CONFIG == 'CREG12.L75':
		# Set the range of values for variables to plot Z-Time
		print()
		print('   >>>>>>>>>>	    Use the CREG12.L75 coordinates for moorings location')
		print()
		# Define a ARC_M1 box or point
		bx_ARC_M1={'name':'Arctic M1 box','imin':1060,'imax':1060,'jmin':1487,'jmax':1487,'depthlim':(-500.,0.),'box':'ARC-M1', 
			   'templim':(-2.,2.,0.2), 'sallim':(27.,35.,0.2),
			   'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(23.5,28.,0.2)
			  }
		#./cdffindij 125. 125. 78. 78. -c coordinates_CREG12_lbclnk_noz_vh20160930.nc
		#  dl_dis=    0.933 km
		#     1060	1060	  1487	    1487
		# 124.9641  124.9641   78.0039	 78.0039
		
		# Define a ARC_A box or point
		bx_ARC_A={'name':'ARC-A','imin':572,'imax':572,'jmin':1541,'jmax':1541,'depthlim':(-500.,0.),'box':'ARC-A',
			  'templim':(-2.,2.,0.2),'sallim':(29.,35,0.2),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(22.,27.5,0.5),
			  'N2lim':(0,5e-4,1e-5)
			 }
		# ./cdffindij -149.58 -149.58 75. 75. -c coordinates_CREG12_lbclnk_noz_vh20160930.nc
		#  dl_dis=    0.560 km
		#      573	 573	  1542	    1542
		# 210.4246  210.4246   75.0049	 75.0049
		
		# Define a ARC_B box or point
		bx_ARC_B={'name':'Arctic B box','imin':627,'imax':627,'jmin':1482,'jmax':1482,'depthlim':(-500.,0.),'box':'ARC-B',
			  'templim':(-2.,2.,0.25),'sallim':(28.,34.4,0.5),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(22.,27.5,0.5)
			 }
		#./cdffindij -150. -150. 78. 78. -c coordinates_CREG12_lbclnk_noz_vh20160930.nc
		#  dl_dis=    1.758 km
		#      627	 627	  1482	    1482
		# 210.0240  210.0240   78.0150	 78.0150
	
		# Define a ARC_C box or point
		bx_ARC_C={'name':'ARC-C','imin':550,'imax':550,'jmin':1466,'jmax':1466,'depthlim':(-500.,0.),'box':'ARC-C',
			  'templim':(-2.,2.,0.2),'sallim':(29.,35,0.2),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(22.,27.5,0.5),
			  'N2lim':(0,5e-4,1e-5)
			 }
		# ./cdffindij -139.54 -139.54 76.59 76.59 -c coordinates_CREG12_lbclnk_noz_vh20160930.nc
		#  dl_dis=    0.997 km
		#      551	 551	  1467	    1467
		# 220.4900  220.4900   76.5957	 76.5957

		# Define a MIK-B box or point (MIKE mooring in the GIN bassin)
		bx_MIK_B={'name':'MIKE box','imin':1067,'imax':1067,'jmin':707,'jmax':707,'depthlim':(-1000.,0.),'box':'MIK-B',
			  'templim':(0.,10.2,0.2),'sallim':(34.9,35.2,0.05),
			  'Wlim':(-10.,10.,2.),'Kzlim':(-6.25,1,0.5),'Sig0lim':(26.5,28.1,0.05)
			 }
		#./cdffindij 2. 2. 66. 66. -c coordinates_CREG12_lbclnk_noz_vh20160930.nc
		#  dl_dis=    2.978 km
		#      1067	 1067	    707       707
		#    2.0280    2.0280	65.9758   65.9758
	
		# Define a ALP-B box or point (ALPHA mooring in the Irminger bassin)
		bx_ALP_B={'name':'ALPHA box','imin':740,'imax':740,'jmin':595,'jmax':595,'depthlim':(-1000.,0.),'box':'ALP-B',
			  'templim':(2.8,8.2,0.2),'sallim':(34.7,35.,0.05),
			  'Wlim':(-10.,10.,2.),'Kzlim':(-5.25,1.,0.25),'Sig0lim':(26.5,27.8,0.05)
			 }
		#./cdffindij -33. -33. 62. 62. -c coordinates_CREG12_lbclnk_noz_vh20160930.nc
		#  dl_dis=    1.710 km
		#	740	  740	    595       595
		#  -33.0323  -33.0323	62.0024   62.0024
	
		# Define a NPEO-B box or point
		bx_NPE_B={'name':'NPEO box','imin':859,'imax':859,'jmin':1245,'jmax':1245,'depthlim':(-500.,0.),'box':'NPE-B',
			  'templim':(-2.2,2.,0.2),'sallim':(30.,34.4,0.2),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(25.5,28.,0.1)
			 }
		#./cdffindij 0. 0. 90 90 -c coordinates_CREG12_lbclnk_noz_vh20160930.nc
		#  dl_dis=    4.935 km
		#      859	 859	  1245	    1245
		# -16.4590  -16.4590   89.9556	 89.9556
	
		# Define a BRAVO-B box or point
		bx_BRA_B={'name':'BRAVO box','imin':521,'imax':521,'jmin':505,'jmax':505,'depthlim':(-2000.,0.),'box':'BRA-B',
			  'templim':(3.,9.,0.2),'sallim':(34.70,34.9,0.05),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(26.5,27.8,0.1)
			 }
		#./cdffindij -52.3 -52.3 56.45 56.45 -c coordinates_CREG12_lbclnk_noz_vh20160930.nc
		#  dl_dis=    3.075 km
		#	521	  521	    505       505
		#  -52.2610  -52.2610	56.4326   56.4326
	
		# Define a GIN-B box or point (Close to FRAM strait)
		bx_GIN_B={'name':'GIN seas box','imin':1061,'imax':1061,'jmin':925,'jmax':925,'depthlim':(-1000.,0.),'box':'GIN-B',
			  'templim':(-2.,8.5,0.5),'sallim':(34.7,34.95,0.025),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-6.,-3.,0.25),'Sig0lim':(26.5,28.5,0.05)
			 }
		#bx_GIN_B={'name':'GIN seas box','imin':1061,'imax':1061,'jmin':925,'jmax':925,'depthlim':(-1000.,0.),'box':'GIN-B','templim':(-2.,8.5,0.5),'sallim':(34.8,35.2,0.5),'Wlim':(-5.,5.,1.),'Kzlim':(-6.,-3.,0.25)}
		#./cdffindij 12.998   12.998   75.107	75.107 -c coordinates_CREG12_lbclnk_noz_vh20160930.nc
		#  dl_dis=    3.264 km
		#      1061	 1061	    925       925
		#   12.9736   12.9736	75.1357   75.1357
	
		# Define a STG-B box or point
		bx_STG_B={'name':'STG box','imin':476,'imax':476,'jmin':144,'jmax':144,'depthlim':(-2000.,0.),'box':'STG-B',
			  'templim':(3.,23.,1.),'sallim':(35.,36.6,0.2),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(26.5,28.5,0.05)
			 }
		#./cdffindij -54.053  -54.053	35.819	 35.819 -c coordinates_CREG12_lbclnk_noz_vh20160930.nc
		#  dl_dis=    5.246 km
		#	476	  476	    144       144
		#  -54.0941  -54.0941	35.7855   35.7855
	
		# Define a Eurasian mooring in the "middle" of this basin
		bx_EUR_B={'name':'EUR box','imin':992,'imax':992,'jmin':1214,'jmax':1214,'depthlim':(-500.,0.),'box':'EUR-B',
			  'templim':(-2.,3.,0.5),'sallim':(32.,35.,0.5),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(26.5,28.5,0.05)
			 }
		#./cdffindij 60. 60. 85.0 85.0 -c coordinates_CREG12_lbclnk_noz_vh20160930.nc -p T
		#   dl_dis=    1.091 km
		#	992	  992	   1214      1214
		#   60.0968   60.0968	85.0050   85.0050
	
	
	
	if CONFIG == 'CREG025.L75':
		print()
		print('   >>>>>>>>>>	    Use the CREG025.L75 coordinates for moorings location')
		print()
		# Define a ARC_M1 box or point
		bx_ARC_M1={'name':'Arctic M1 box','imin':354,'imax':354,'jmin':498,'jmax':498,'depthlim':(-500.,0.),'box':'ARC-M1',
			   'templim':(-2.,2.,0.2),'sallim':(27.,35.,0.2),
			   'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(23.5,28.,0.2)
			  }
		# ./cdffindij 125. 125. 78. 78. -c CREG025.L75_coordinates.nc -p T
		#    dl_dis=	8.699 km
		#	 354	   354	     498       498
		#   125.3725  125.3725	 77.9889   77.9889
		
		# Define a ARC_A box or point
		bx_ARC_A={'name':'ARC-A','imin':191,'imax':191,'jmin':515,'jmax':515,'depthlim':(-500.,0.),'box':'ARC-A',
			  'templim':(-2.,2.,0.2),'sallim':(29.,35,0.2),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(22.,27.5,0.5),
			  'N2lim':(0,5e-4,1e-5)
			 }
		# ./cdffindij -149.58 -149.58 75. 75. -c  CREG025.L75_coordinates.nc -p T 
		#  dl_dis=    2.406 km
		#      192	 192	   516	     516
		# 210.4358  210.4358   74.9787	 74.9787

		# Define a ARC_B box or point
		bx_ARC_B={'name':'Arctic B box','imin':210,'imax':210,'jmin':496,'jmax':496,'depthlim':(-500.,0.),'box':'ARC-B',
			  'templim':(-2.,2.,0.25),'sallim':(28.,34.4,0.5),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(22.,27.5,0.5)
			 }
		# ./cdffindij -150. -150. 78. 78. -c  CREG025.L75_coordinates.nc -p T 
		#    dl_dis=	1.435 km
		#	 210	   210	     496       496
		#   210.0315  210.0315	 77.9889   77.9889
	
		# Define a ARC_C box or point
		bx_ARC_C={'name':'ARC-C','imin':184,'imax':184,'jmin':490,'jmax':490,'depthlim':(-500.,0.),'box':'ARC-C',
			  'templim':(-2.,2.,0.2),'sallim':(29.,35,0.2),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(22.,27.5,0.5),
			  'N2lim':(0,5e-4,1e-5)
			 }
		# ./cdffindij -139.54 -139.54 76.59 76.59 -c  CREG025.L75_coordinates.nc -p T 
		#  dl_dis=    2.120 km
		#      185	 185	   491	     491
		# 220.3823  220.3823   76.5962	 76.5962

		# Define a MIK-B box or point (MIKE mooring in the GIN bassin)
		bx_MIK_B={'name':'MIKE box','imin':356,'imax':356,'jmin':237,'jmax':237,'depthlim':(-1000.,0.),'box':'MIK-B',
			  'templim':(0.,10.2,0.2),'sallim':(34.9,35.2,0.05),
			  'Wlim':(-10.,10.,2.),'Kzlim':(-6.25,1,0.5),'Sig0lim':(26.5,28.1,0.05)
			 }
		#./cdffindij 2. 2. 66. 66. -c CREG025.L75_coordinates.nc
		#  dl_dis=    3.657 km
		#	356	  356	    237       237
		#    1.9244    1.9244	65.9883   65.9883
	
		# Define a ALP-B box or point (ALPHA mooring in the Irminger bassin)
		bx_ALP_B={'name':'ALPHA box','imin':247,'imax':247,'jmin':200,'jmax':200,'depthlim':(-1000.,0.),'box':'ALP-B',
			  'templim':(2.8,8.2,0.2),'sallim':(34.7,35.,0.05),
			  'Wlim':(-10.,10.,2.),'Kzlim':(-5.25,1.,0.25),'Sig0lim':(26.5,27.8,0.05)
			 }
		#./cdffindij -33. -33. 62. 62. -c CREG025.L75_coordinates.nc
		#  dl_dis=    8.947 km
		#	247	  247	    200       200
		#  -33.1371  -33.1371	62.0484   62.0484
	
		# Define a NPEO-B box or point
		bx_NPE_B={'name':'NPEO box','imin':287,'imax':287,'jmin':416,'jmax':416,'depthlim':(-500.,0.),'box':'NPE-B',
			  'templim':(-2.2,2.,0.2),'sallim':(30.,34.4,0.2),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(25.5,28.,0.1)
			 }
		#./cdffindij 0. 0. 90 90 -c CREG025.L75_coordinates.nc
		#  dl_dis=    6.930 km
		#	287	  287	    416       416
		#  -17.6907  -17.6907	89.9376   89.9376
	
		# Define a BRAVO-B box or point
		bx_BRA_B={'name':'BRAVO box','imin':174,'imax':174,'jmin':170,'jmax':170,'depthlim':(-2000.,0.),'box':'BRA-B',
			  'templim':(2.,4.,0.1),'sallim':(34.60,35.0,0.05),
			  #'templim':(3.,9.,0.2),'sallim':(34.70,34.9,0.05),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(26.5,27.8,0.1)
			 }
		#./cdffindij -52.3 -52.3 56.45 56.45 -c CREG025.L75_coordinates.nc
		#  dl_dis=    4.580 km
		#	174	  174	    170       170
		#  -52.3611  -52.3611	56.4736   56.4736
	
		# Define a GIN-B box or point (Close to FRAM strait)
		bx_GIN_B={'name':'GIN seas box','imin':355,'imax':355,'jmin':310,'jmax':310,'depthlim':(-1000.,0.),'box':'GIN-B',
			  'templim':(0.,5.0,0.5),'sallim':(34.7,35.10,0.05),
			  #'templim':(-2.,8.5,0.5),'sallim':(34.7,34.95,0.025),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-6.,-3.,0.25),'Sig0lim':(26.5,28.5,0.05)
			 }
		#./cdfwhereij 355 355 310 310 -c CREG025.L75_coordinates.nc
		# Type of point   : T
		#   I J zoom	  :    355   355   310	 310
		#   LON LAT zoom  :    12.998	12.998	 75.107   75.107
	
		# Define a STG-B box or point
		bx_STG_B={'name':'STG box','imin':160,'imax':160,'jmin':50,'jmax':50,'depthlim':(-2000.,0.),'box':'STG-B',
			  'templim':(3.,23.,1.),'sallim':(35.,36.6,0.2),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(26.5,28.5,0.05)
			 }
		#./cdfwhereij 160 160 50 50 -c CREG025.L75_coordinates.nc
		# Type of point   : T
		#   I J zoom	  :    160   160    50	  50
		#   LON LAT zoom  :   -54.053  -54.053	 35.819   35.819
	
		# Define a Eurasian mooring in the "middle" of this basin
		bx_EUR_B={'name':'EUR box','imin':332,'imax':332,'jmin':406,'jmax':406,'depthlim':(-500.,0.),'box':'EUR-B',
			  'templim':(-2.,3.,0.5),'sallim':(32.,35.,0.5),
			  'Wlim':(-5.,5.,1.),'Kzlim':(-7.,-5.,0.25),'Sig0lim':(26.5,28.5,0.05)
			 }
		#./cdfwhereij 60 60 85 85 -c CREG025.L75_coordinates.nc
		# Type of point   : T
		#  dl_dis=    4.787 km
		#      332	 332	   406	     406
		#  59.7941   59.7941   84.9609	 84.9609
	
	if box2sel == 'ARCM1':
		sbox=bx_ARC_M1
	elif box2sel == 'ARCA':
		sbox=bx_ARC_A
	elif box2sel == 'ARCB':
		sbox=bx_ARC_B
	elif box2sel == 'ARCC':
		sbox=bx_ARC_C
	elif box2sel == 'GINB':
		sbox=bx_GIN_B
	elif box2sel == 'STGB':
		sbox=bx_STG_B
	elif box2sel == 'MIKB':
		sbox=bx_MIK_B
	elif box2sel == 'ALPB':
		sbox=bx_ALP_B
	elif box2sel == 'NPEB':
		sbox=bx_NPE_B
	elif box2sel == 'BRAB':
		sbox=bx_BRA_B
	elif box2sel == 'EURA':
		sbox=bx_EUR_B

	return sbox

def DEF_ZPROFILE(zCONFIG,zCASE,zclimyear,zhsct_lev,Red_My_varinit,box,zAll_var,zplt,zgtype=None,zMyvar=None,zfram=None,zs_year=None):

	nvar=0
	for var in zAll_var :
		print('   Vertical profile for '+ var['name']+' '+ box['box'] )
		
		ax=plt.subplot(zfram+nvar)
		if var['name'] == "votemper":
			if zfram == 221 : 
				plt.title(zclimyear+' mean temperature \n'+' @ moor '+box['name'],size=7)
			else:
				plt.title('\n @ moor '+box['name'],size=7)
				plt.xlabel(r' ($^\circ$C)', fontsize=7)
			if box['box'] == 'ARC-B':
				plt.axis([-2.,2.,-1000.,0.])
				plt.yticks(-100.*npy.arange(11),size=7)
			elif box['box'] == 'EUR-B':
				plt.axis([-2.,2.,-2000.,0.])
				plt.yticks(-200.*npy.arange(11),size=7)
			elif box['box'] == 'GIN-B':
				print(	box['box']+ ' set axis Ok')
				plt.axis([-2.,6.,-1000.,0.])
				plt.yticks(-100.*npy.arange(11),size=7)
			elif box['box'] == 'BRA-B':
				print(	box['box']+ ' set axis Ok')
				plt.axis([0.,5.,-1000.,0.])
				plt.yticks(-100.*npy.arange(11),size=7)
			elif box['box'] == 'MIK-B':
				print(	box['box']+ ' set axis Ok')
				plt.axis([-2.,8.,-1000.,0.])
				plt.yticks(-100.*npy.arange(11),size=7)
			elif box['box'] == 'ALP-B':
				print(	box['box']+ ' set axis Ok')
				plt.axis([2.,8.,-2000.,0.])
				plt.yticks(-200.*npy.arange(11),size=7)

			plt.xticks(size=7)
			plt.yticks(size=7)
			plt.ylabel('Depth (m)',size=7)

		elif var['name'] == "vosaline":
			if zfram == 221 : 
				plt.title(zclimyear+' mean Salinity \n'+' @ moor '+box['name'],size=7)
			else:
				plt.title('\n @ moor '+box['name'],size=7)
				plt.xlabel(r' (PSU)', fontsize=7)
			if box['box'] == 'ARC-B':
				plt.axis([28.,35.5,-1000.,0.])
				plt.yticks(-100.*npy.arange(11),size=7)
			elif box['box'] == 'EUR-B':
				plt.axis([28.,35.5,-2000.,0.])
				plt.yticks(-200.*npy.arange(11),size=7)
			elif box['box'] == 'GIN-B':
				print(	box['box']+ ' set axis Ok')
				plt.axis([34.8,35.2,-1000.,0.])
				plt.yticks(-100.*npy.arange(11),size=7)
			elif box['box'] == 'BRA-B':
				print(	box['box']+ ' set axis Ok')
				plt.axis([34.,35.,-1000.,0.])
				plt.yticks(-100.*npy.arange(11),size=7)
			elif box['box'] == 'MIK-B':
				print(	box['box']+ ' set axis Ok')
				plt.axis([34.8,35.4,-1000.,0.])
				plt.yticks(-100.*npy.arange(11),size=7)
			elif box['box'] == 'ALP-B':
				print(	box['box']+ ' set axis Ok')
				plt.axis([34.8,35.2,-2000.,0.])
				plt.yticks(-200.*npy.arange(11),size=7)

			plt.xticks(size=7)
			plt.yticks(size=7)
			plt.setp(ax.get_yticklabels(),visible=False)
	
		# Plot just the mean at each depth to plot the vertical profile
		mean_z=npy.mean(zhsct_lev,axis=2).squeeze()

		if var['name'] == "vovecrtz":
			plt.plot(1e6*mean_z[nvar,:]   ,-zplt[:,0],'k', linewidth=0.7, label=zCASE)
		elif var['name'] == "votkeavt":
			plt.semilogx(mean_z[nvar,:],-zplt[:,0],'k', linewidth=0.7, label=zCASE)
			#plt.plot(npy.log10(mean_z[nvar,:])	    ,-zplt[:,0],'k', linewidth=0.7, label=zCASE)
		elif var['name'] == "rhop_sig0":
			plt.plot(mean_z[nvar,:]-1.e+3	,-zplt[:,0],'k', linewidth=0.7, label=zCASE)
		else:
			plt.plot(mean_z[nvar,:]   ,-zplt[:,0],'k', linewidth=0.7, label=zCASE)
			plt.plot(Red_My_varinit[nvar,:] ,-zplt[:,0],'m', linewidth=0.7, label='woa09')

		if zfram == 221 :
			locpath='./DATA/'
			locfile='TimeMean_profiles_moorings.mat'
			if chkfile(locpath+locfile) :
				B_mooring = sio.loadmat(locpath+locfile,squeeze_me=True)
			else:
				B_mooring = npy.zeros(mean_z.shape)
			if nvar == 0 :
				if box['box'] == 'ARC-B' : plt.plot(B_mooring['TB_mean'],-1.*B_mooring['Zint'],'g', linewidth=0.7, label='Obs')
			else:
				if box['box'] == 'ARC-B' : plt.plot(B_mooring['SB_mean'],-1.*B_mooring['Zint'],'g', linewidth=0.7, label='Obs')
		plt.grid(True,linestyle='--',color='grey',alpha=0.7)
	
		if var['name'] == "vosaline" :
			plt.legend(loc='lower left')
			leg = plt.gca().get_legend()
			ltext  = leg.get_texts()
			plt.setp(ltext, fontsize=7.)
	
	
		nvar+=1
    
	# Save Vertical profiles in Python npy format
	file_extZ='_ZTSProfile'
	npy.save('./DATA/'+zCONFIG+'/'+zCASE+'_'+box['box']+file_extZ+'_y'+str(zs_year)+'.npy',mean_z)
	
	return 


def DEF_ZTIME(zCONFIG,zCASE,lgTS_ys,lgTS_ye,box,z,z2dt,hsct_lev,zgtype=None,zMyvar=None,zfram=None,zoutNC=False):

	LongTS_hsct_lev= []   

	# Start to read all yearly files
	################################
	lgts_year=lgTS_ys    ;	  t_months=(npy.arange(12)*30.+15.)/365.   ;   start = 1
	while  lgts_year <= lgTS_ye  :
	       # Read monthly Time-depth fields for each year
	       print('		>>>>		       Read LONG TIME SERIES PROFILES | year '+	str(lgts_year))
	       file_extZ='_ZTimeTS'
	       locpath='./DATA/'+zCONFIG+'/'
	       locfile=zCASE+'_'+box['box']+file_extZ+'_y'+str(lgts_year)+'.npy'
	       if chkfile(locpath+locfile) :
		        year_hsct_lev = npy.load(locpath+locfile,mmap_mode='r')
	       else:
		        year_hsct_lev = npy.zeros(hsct_lev.shape)+npy.nan

	       # Set the time axis
	       y_years=npy.tile(lgts_year,12)+t_months
	       if start == 1:
		       time_axis=y_years
		       start=0
		       LongTS_hsct_lev=year_hsct_lev.copy()
	       else:
		       time_axis=npy.append(time_axis,y_years)
		       LongTS_hsct_lev=npy.concatenate((LongTS_hsct_lev,year_hsct_lev),axis=2)

	       lgts_year+=1

	xplt= npy.tile(time_axis,(z.size,1))
	zplt = npy.repeat(z2dt,time_axis.shape[0],axis=1)
	lgtsclimyear=str(lgTS_ys)+str(lgTS_ye)

	time_grid=npy.arange(lgTS_ys,lgTS_ye+1,1.,dtype=float)
	newlocsx  = npy.array(time_grid,'f')
	newlabelsx = npy.array(time_grid,'i')

	# Plot the Time-depth time-series over SEVERAL YEARS
	########################################################
	lnvar=0
	for var in DEF_INFO_VAR('gridT','votemper') :
		ax=plt.subplot(zfram+lnvar)

		contours,norm,cmap,isoline = SET_CONT(box,var['name'])

		C = plt.contourf(xplt[:,:],-zplt[:,:],LongTS_hsct_lev[lnvar,:,:],contours,norm=norm,cmap=cmap,extend='both')
		if var['name'] == 'votemper':
			CS=plt.contour(xplt[:,:],-zplt[:,:],LongTS_hsct_lev[lnvar,:,:],colors='k',linewidths=0.5,levels=[0.,1.])
			plt.clabel(CS,fontsize=6,fmt='%1.0f')
		if var['name'] == "vosaline" :
			CS=plt.contour(xplt[:,:],-zplt[:,:],LongTS_hsct_lev[lnvar,:,:],colors='k',linewidths=0.5,levels=[31.,34.,34.8])
			plt.clabel(CS,fontsize=6,fmt='%4.1f')

		#plt.grid(True,linestyle='--',color='grey',alpha=0.8)
		plt.yticks(size=6)
		plt.xlim([lgTS_ys,lgTS_ye+1])
		plt.ylim([box['depthlim'][0],box['depthlim'][1]])
		plt.ylabel('Depth (m)',fontsize=6)

		if zfram == 411 and var['name'] == "votemper" :
			plt.xticks(newlocsx,newlabelsx,size=6)
			plt.setp(ax.get_xticklabels(),visible=False)
			plt.title(zCASE+' '+box['name'],fontsize=7)
			zmy_cblab=r'($^\circ$C)'
		elif zfram == 411 and var['name'] == "vosaline" :
			plt.xticks(newlocsx,newlabelsx,size=6)
			plt.setp(ax.get_xticklabels(),visible=False)
			zmy_cblab=r'(PSU)'
		elif zfram != 411 and var['name'] == "votemper" :
			plt.title(zCASE+' '+box['name'],fontsize=7)
			plt.setp(ax.get_xticklabels(),visible=False)
			plt.xticks(newlocsx,newlabelsx,size=6)
			zmy_cblab=r'($^\circ$C)'
		elif zfram != 411 and var['name'] == 'vosaline' :
			plt.setp(ax.get_xticklabels(),rotation=90)
			plt.xticks(newlocsx,newlabelsx,size=5)
			plt.xlabel('Years', fontsize=6)
			zmy_cblab=r'(PSU)'

		#divider = make_axes_locatable(ax)
		#cax = divider.append_axes("right", size="5%", pad=0.05)
		#cbar = plt.colorbar(C,cax=cax)

		cbar=plt.colorbar(C,format='%.1f',orientation='vertical',drawedges=True)
		#plt.subplots_adjust(hspace=0.5)
		#cbar = plt.colorbar(C,format='%.1f',shrink=1.5,orientation='vertical',drawedges=True)
		#posx=0.98  ; posy=0.25
		#cax=plt.axes([posx,posy,0.02,0.55])
		#cbar = plt.colorbar(C,cax,format='%.1f',orientation='vertical',drawedges=True)
		##cbar = plt.colorbar(C,format='%.1f',orientation='vertical',shrink=0.8,drawedges=True)
		cbar.set_label(zmy_cblab,fontsize=7)
		cl = plt.getp(cbar.ax, 'ymajorticklabels')
		plt.setp(cl, fontsize=6)

		lnvar+=1

	if zoutNC:
		ds_outZProf = xr.Dataset()
		# Define coordinates
		ds_outZProf.coords['time'] = (('time') , time_axis.astype('float32'))
		ds_outZProf.coords['z']    = (('z')    , z.values.astype('float32'))

		ds_outZProf['Temp']= (('z','time'), LongTS_hsct_lev[0,:,:].astype('float32')) 
		ds_outZProf['Temp'].attrs['long_name']='Monthly mean temperature profile @ mooring '+box['box']
		ds_outZProf['Temp'].attrs['units']='DegC'

		ds_outZProf['Sal']= (('z','time'), LongTS_hsct_lev[1,:,:].astype('float32')) 
		ds_outZProf['Sal'].attrs['long_name']='Monthly mean salinity profile @ mooring '+box['box']
		ds_outZProf['Sal'].attrs['units']='PSU'

		ds_outZProf['Time2D']= (('z','time'), xplt ) 
		ds_outZProf['Time2D'].attrs['long_name']='Depth versus time array for contour plot'
		ds_outZProf['Time2D'].attrs['units']='month'

		ds_outZProf['Depth2D']= (('z','time'), -1.*zplt ) 
		ds_outZProf['Depth2D'].attrs['long_name']='Depth versus time array for contour plot'
		ds_outZProf['Depth2D'].attrs['units']='m'

		# Write the NetCDF file 
		ds_outZProf.attrs['History'] = "Diagnostics have been calculated using the Arctic monitoring tool "
		ds_outZProf.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
		nc_f = './NETCDF/'+zCONFIG+'-'+zCASE+'_MOOR-'+box['box']+'_y'+str(lgTS_ys)+'LASTy.nc'
		ds_outZProf.to_netcdf(nc_f,engine='netcdf4')


	return 


def DEF_REDVAR(CONFIG,CASE,My_var,My_varinit,box,zAll_var,s_year,e_year,tmask,e1t,e2t,e3t):
	########################################
	# Start calculation
	########################################
	#------------------------------------------------------------------------------------------------------------------------
	# Number of variables temperature/salinity to process
	nbvar = 2 

	# Define a box area in which mean values will be computed
	boxmsk = tmask.copy()*.0
	boxmsk[:,box['jmin']:box['jmax']+1,box['imin']:box['imax']+1] = 1.

	# Compute the volume of each grid point
	My_box_vol = e3t * boxmsk * e1t * e2t
	Red_My_box_vol = My_box_vol.isel(y=slice(box['jmin'],box['jmax']+1),x=slice(box['imin'],box['imax']+1))
	box_vol_lev = Red_My_box_vol.sum(dim=['y','x'])
	print()
	print('			Reduce Array size step')
	print()

	Red_My_var = My_var.isel(y=slice(box['jmin'],box['jmax']+1),x=slice(box['imin'],box['imax']+1))
	My_varinitbox = My_varinit.isel(y=slice(box['jmin'],box['jmax']+1),x=slice(box['imin'],box['imax']+1)).squeeze()
	
	hsct_lev = xr.DataArray( dims=['nvar','z','time_counter'], coords = dict(nvar=npy.arange(nbvar),z=My_var.z,time_counter=My_var.time_counter) ) 
	Red_My_varinit = xr.DataArray( dims=['nvar','z'], coords = dict(nvar=npy.arange(nbvar),z=My_var.z) ) 
	
	nvar=0
	for var in zAll_var:
	
		print()
		print('		        Compute HSCT ZT for '+var['name'])
		print()
		# Compute the Temp/Salt mean
		hsct_lev[nvar,:,:] = (Red_My_var[var['name']] * Red_My_box_vol).sum(dim=['y','x'])/ box_vol_lev
		Red_My_varinit[nvar,:] = My_varinitbox[var['name']].values
		nvar+=1
	
	# Save year time series to use for a longer time series
	file_extZ='_ZTimeTS'
	npy.save('./DATA/'+CONFIG+'/'+CASE+'_'+box['box']+file_extZ+'_y'+str(s_year)+'.npy',hsct_lev)
	
	return hsct_lev, Red_My_varinit
