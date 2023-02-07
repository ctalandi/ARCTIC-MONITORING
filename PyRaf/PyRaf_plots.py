##########################################################################
### Copyright (c) 2011 by Raphael Dussin. For licensing, distribution 
### conditions, contact information, and additional documentation see
### in docs directory

##########################################################################
###                                                                    ###
###                          This is PyRaf (Plots)                     ###
###                                                                    ###
##########################################################################
### IMPORT the packages

import PyRaf
import PyRaf_colormaps

# Numpy
try:
        import numpy as npy
except:
        print 'numpy is not available on your machine'
	print 'check python path or install this package' ; exit()

# Matplotlib
try:
	import matplotlib.pylab as plt
	import matplotlib.cm as cm
	import matplotlib as mpl
	from matplotlib import rc
	from matplotlib import rcParams
except:
	print 'matplotlib is not available on your machine'
        print 'check python path or install this package' ; exit()

# Basemap
try:
#	from mpl_toolkits.basemap import Basemap
	from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
except:
	print 'Basemap is not available on your machine'
	print 'check python path or install this package' ; exit()


#######################################################
### PLOTS
#######################################################

def nemo_global_plot(lon,lat,tab,contours,limits,myticks=None,name=None,filename='test.pdf'):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
#	rcParams['text.usetex']=True
#	rcParams['text.latex.unicode']=True
	rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=PyRaf.limit_range(tab,limits[0],limits[1])
	fig = plt.figure(figsize=[12.8 , 6.])
	ax  = fig.add_subplot(111)
	# background
	m = Basemap(projection='cyl',llcrnrlat=-89,urcrnrlat=89,\
	            llcrnrlon=-180,urcrnrlon=180,resolution='c')
	m.drawcoastlines()
	m.fillcontinents(color='grey',lake_color='white')
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	if len(limits) > 2:
		norm = mpl.colors.Normalize(vmin=limits[2], vmax=limits[3])
	pal = PyRaf_colormaps.gen_pal_Testu()
	C = m.contourf(lon,lat,tab,25,cmap=pal,norm=norm)
	# colorbar	
	#if myticks is None:
	#cbar = plt.colorbar(C,format='%.2f',orientation='horizontal',shrink=0.8)
	#else:
	#cbar = plt.colorbar(C,format='%.2f',orientation='horizontal',shrink=0.8,ticks=myticks)
	plt.clim(limits[0],limits[1])
	if len(limits) > 2:
        	plt.clim(limits[2],limits[3])
	# contour (optional)
	if len(contours) > 0:
		CS2 = m.contour(lon, lat, tab, contours, colors='k')
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=8)
	#
	# x axis
	locs, labels = plt.xticks()
	newlocs   = npy.array([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180],'f')
	newlabels = npy.array([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180],'i')
	plt.xticks(newlocs,newlabels)
	plt.xlabel('Longitude',fontsize=16)
	#
	# y axis
	locsY,labelsy = plt.yticks()
	newlocsy   = npy.array([-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80],'f')
	newlabelsy = npy.array([-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80],'i')
	plt.yticks(newlocsy,newlabelsy)
	plt.ylabel('Latitude',fontsize=16)	
	# title
	plt.title(name,fontsize=18)
	plt.savefig(filename)
	return fig

def nemo_Arc_plot(lon,lat,tab,contours,limits,myticks=None,name=None,zmy_cblab=None,zmy_cmap=None,filename='test.pdf',zvar=None):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
#	rcParams['text.usetex']=True
#	rcParams['text.latex.unicode']=True
	rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=PyRaf.limit_range(tab,limits[0],limits[1])
	#fig = plt.figure(figsize=[12.8 , 6.])
	#ax  = fig.add_subplot(111)
	# background
	m = Basemap(projection='npstere',boundinglat=65,lon_0=-60, resolution='i')
	#m = Basemap(width=12000000,height=8000000, \
	#            resolution='i',projection='laea',\
	#            lat_ts=70,lat_0=60,lon_0=-45.)

	#m.drawcoastlines(linewidth=0.25)
	if zvar != 'Bathy' :
		m.drawparallels(npy.arange(-90.,91.,5.),labels=[False,False,False,False])
		m.drawmeridians(npy.arange(-180.,181.,20.),labels=[True,False,False,True])
		m.fillcontinents(color='grey',lake_color='white')
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])

	if zmy_cmap != None :
		pal = zmy_cmap
	else:
		pal = plt.cm.get_cmap('coolwarm')
		#pal = plt.cm.get_cmap('terrain')

	X,Y = m(lon,lat)
	C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm,extend='both')
        if zvar == 'vosaline':
                CS2 = m.contour(X,Y,tab,linewidths=1.,levels=[34.8], colors='g',alpha=0.6)
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=8, colors='g')
        if zvar == 'votemper':
                CS2 = m.contour(X,Y,tab,linewidths=1.,levels=[3.], colors='m',alpha=0.8)
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.0f', fontsize=8, colors='m')
                CS3 = m.contour(X,Y,tab,linewidths=1.,levels=[0.], colors='w',alpha=0.8)
		plt.clabel(CS3, CS3.levels, inline=True, fmt='%.0f', fontsize=8, colors='w')

	############################################################################################################
	############################################################################################################
	moorplot=1
	if moorplot == 1 :
        	bx_ARCB={'name':'B' ,'lon_min':-150.,'lon_max':-150.,'lat_min':78.,'lat_max':78.}
		bx_ARCM={'name':'M1','lon_min': 125.,'lon_max': 125.,'lat_min':78.,'lat_max':78.}

		All_box=[bx_ARCB,bx_ARCM]
		for box in All_box:
        		lats = [box['lat_min'],box['lat_max']]
        		lons = [box['lon_min'],box['lon_max']]
        		x,y = m(lons,lats)
        		m.scatter(x,y,6,marker='o', color='w')
        		#m.plot(x,y,linewidth=2, color='g')
	############################################################################################################
	############################################################################################################

	# colorbar	
	if myticks is None:
		cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
	else:
		cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8,drawedges=True)
		cbar.set_label(zmy_cblab,fontsize=9)
		cl = plt.getp(cbar.ax, 'ymajorticklabels')
		plt.setp(cl, fontsize=9)

	#plt.clim(limits[0],limits[1])
	# contour (optional)
	#if len(contours) > 0:
	#	CS2 = m.contour(X, Y, tab, myticks, colors='k',linewidths=0.3 )
	#	plt.clabel(CS2, CS2.levels, inline=True, fmt='%.0f', fontsize=8)

	#plt.grid(True)
	#
	# x axis
	#locs, labels = plt.xticks()
	#newlocs   = npy.array([-90,-80,-70,-60,-50,-40,-30,-20,-10,0],'f')
	#newlabels = npy.array([-90,-80,-70,-60,-50,-40,-30,-20,-10,0],'i')
	#Plt.xticks(newlocs,newlabels,size=8)
	#plt.xlabel('Longitude',fontsize=10)
	#
	# y axis
	#locsY,labelsy = plt.yticks()
	#newlocsy   = npy.array([20,30,40,50,60,70,80],'f')
	#newlabelsy = npy.array([20,30,40,50,60,70,80],'i')
	#plt.yticks(newlocsy,newlabelsy,size=8)
	#plt.ylabel('Latitude',fontsize=10)	
	# title
	plt.title(name,fontsize=11)
	#plt.savefig(filename)
	#plt.show()

	return 


def nemo_atl_plot(lon,lat,tab,contours,limits,myticks=None,name=None,filename='test.pdf'):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
#	rcParams['text.usetex']=True
#	rcParams['text.latex.unicode']=True
	rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=PyRaf.limit_range(tab,limits[0],limits[1])
	#fig = plt.figure(figsize=[12.8 , 6.])
	#ax  = fig.add_subplot(111)
	# background
	m = Basemap(projection='cyl',llcrnrlat=20,urcrnrlat=80,\
	            llcrnrlon=-90,urcrnrlon=0,resolution='c')
	#m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=-20,\
	#            llcrnrlon=-180,urcrnrlon=180,resolution='c')
	m.drawcoastlines(linewidth=0.25)
	m.fillcontinents(color='grey',lake_color='white')
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	if len(limits) > 2:
		norm = mpl.colors.Normalize(vmin=limits[2], vmax=limits[3])
	#pal = PyRaf_colormaps.gen_pal_Testu2()
	#pal = PyRaf_colormaps.gen_pal_blue2red3()
	#pal = plt.cm.get_cmap('jet')
	pal = plt.cm.get_cmap('seismic')
	pal = plt.cm.get_cmap('RdBu_r')
	# for EKE plot pal = plt.cm.get_cmap('BuPu')
	X,Y = m(lon,lat)
	#C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm)
	C = m.contourf(lon,lat,tab,25,cmap=pal,norm=norm)
        ## usefull for bathymetry contours #CC =m.contour(lon,lat,tab,colors='k',linewidths=.5,levels=[1000.,2000.,3000.,4000.])
        ## usefull for bathymetry contours #plt.clabel(CC,fontsize=7,fmt='%5.1f')
	# colorbar	
	#if myticks is None:
	#	cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
	#else:
	#	cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8,ticks=myticks)
	plt.clim(limits[0],limits[1])
	if len(limits) > 2:
        	plt.clim(limits[2],limits[3])
	# contour (optional)
	if len(contours) > 0:
		CS2 = m.contour(X, Y, tab, contours, colors='k',linewidths=0.3 )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=5)
		CS2 = m.contour(X, Y, tab, linewidths=0.6,levels=[0.,0.], colors='k' )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=5)
	plt.grid(True)
	#
	# x axis
	locs, labels = plt.xticks()
	newlocs   = npy.array([-90,-80,-70,-60,-50,-40,-30,-20,-10,0],'f')
	newlabels = npy.array([-90,-80,-70,-60,-50,-40,-30,-20,-10,0],'i')
	#newlocs   = npy.array([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180],'f')
	#newlabels = npy.array([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180],'i')
	plt.xticks(newlocs,newlabels,size=8)
	plt.xlabel('Longitude',fontsize=10)
	#
	# y axis
	locsY,labelsy = plt.yticks()
	newlocsy   = npy.array([20,30,40,50,60,70,80],'f')
	newlabelsy = npy.array([20,30,40,50,60,70,80],'i')
	plt.yticks(newlocsy,newlabelsy,size=8)
	plt.ylabel('Latitude',fontsize=10)	
	# title
	plt.title(name,fontsize=10)
	#plt.savefig(filename)
	#plt.show()
	#return fig

def nemo_Atl_Bat_plot(lon,lat,tab,contours,limits,myticks=None,name=None,area=None,filename='test.pdf',type=None):
	#
#	rcParams['text.usetex']=True
	rcParams['text.latex.unicode']=True
 #	rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=PyRaf.limit_range(tab,limits[0],limits[1])
	# background
	#m = Basemap(width=4000000,height=3500000,lat_1=30.,lat_2=70,lon_0=-40,lat_0=45,projection='aea',resolution='c')
	#m = Basemap(projection='ortho',lat_0=60.,lon_0=-40.,resolution='c')
	# m = Basemap(width=2000000,height=1150000,lat_1=50.,lat_2=70,lon_0=-20,lat_0=65,projection='aea',resolution='c')
	if type == 'spec1':
	    m = Basemap(width=800000,height=500000,lat_1=55.,lat_2=70,lon_0=-50,lat_0=61,projection='aea',resolution='c')
	elif type == 'labsea':
		# To plot the isopycne depth over Labrador sea erna_SigiDepth.py script
		m = Basemap(width=1500000,height=1600000,lat_1=50.,lat_2=65,lon_0=-50,lat_0=59.5,projection='aea',resolution='i')
		#m = Basemap(width=1400000,height=1400000,lat_1=54.,lat_2=65,lon_0=-50,lat_0=58.5,projection='aea',resolution='i')
	elif type == 'ar7zW':
		# To plot the isopycne depth over Labrador sea erna_SigiDepth.py script
		m = Basemap(width=400000,height=500000,lat_1=54.,lat_2=56,lon_0=-54,lat_0=55.0,projection='aea',resolution='c')
	elif type == 'ar7zE':
		# To plot the isopycne depth over Labrador sea erna_SigiDepth.py script
		m = Basemap(width=400000,height=600000,lat_1=59.,lat_2=62,lon_0=-49,lat_0=61.0,projection='aea',resolution='c')
        elif type == 'arc':
            # To plot the isopycne depth over Labrador sea erna_SigiDepth.py script
            m = Basemap(projection='npstere',boundinglat=65,lon_0=-60, resolution='i')
	elif type == 'rrex':
		# To plot the isopycne depth over Labrador sea erna_SigiDepth.py script
		m = Basemap(width=9000000,height=8000000,lat_1=20.,lat_2=80.,lon_0=-35,lat_0=50.0,projection='aea',resolution='c')
	elif type == 'atln':
		# To plot the North Atlantic 
		m = Basemap(width=6000000,height=5500000,lat_1=25.,lat_2=75.,\
                            lon_0=-50,lat_0=50,projection='aea',resolution='i')
	else:
	    m = Basemap(projection='cyl',llcrnrlat=area['latmin'],urcrnrlat=area['latmax'],\
	                llcrnrlon=area['lonmin'],urcrnrlon=area['lonmax'],resolution='c')
	#use this one for erna_VMOD.py m = Basemap(projection='cyl',llcrnrlat=area['latmin'],urcrnrlat=area['latmax'],\
	#use this one for erna_VMOD.py             llcrnrlon=area['lonmin'],urcrnrlon=area['lonmax'],resolution='c')
	#m.drawcoastlines(linewidth=0.25)
	#m.fillcontinents(color='grey',lake_color='white')
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	if len(limits) > 2:
		norm = mpl.colors.Normalize(vmin=limits[2], vmax=limits[3])
	pal = plt.cm.get_cmap('binary')
	#pal = plt.cm.get_cmap('cool_r')
	X,Y = m(lon,lat)
	if type == 'spec1': C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm,extend='both')

	############################################################################################################
	############################################################################################################
	secplot=0
	if secplot == 1 :
        	bx_FISS={'name':'FIS','lon_min':-51.6,'lon_max':-49.,'lat_min':52.7,'lat_max':53.7}
		bx_AR7S={'name':'AR7','lon_min':-55.5,'lon_max':-48.,'lat_min':53.4,'lat_max':61. }
		bx_GRES={'name':'GRE','lon_min':-44. ,'lon_max':-44.,'lat_min':55. ,'lat_max':60. }
		bx_IRMS={'name':'IRM','lon_min':-38. ,'lon_max':-32.,'lat_min':66. ,'lat_max':63. }
		bx_DSOS={'name':'DSO','lon_min':-30. ,'lon_max':-25.,'lat_min':67. ,'lat_max':65.5}
		bx_SCHS={'name':'SCH','lon_min':-49.5,'lon_max':-45.,'lat_min':43.1,'lat_max':42. }
		bx_LIWS={'name':'LIW','lon_min':-70. ,'lon_max':-68.,'lat_min':40. ,'lat_max':37.5}
		# # SEC LIW
		# ybox7 = (/ 40.0, 37.5, 37.5, 40., 40./)
		# xbox7 = (/ -70., -68., -68., -70., -70./)

		All_box=[bx_FISS,bx_AR7S,bx_GRES,bx_IRMS,bx_DSOS,bx_SCHS,bx_LIWS]
		#All_box=[bx_FISS]
		for box in All_box:
        		lats = [box['lat_min'],box['lat_max']]
        		lons = [box['lon_min'],box['lon_max']]
        		x,y = m(lons,lats)
        		m.plot(x,y,linewidth=2, color='r')
	############################################################################################################
	############################################################################################################

	# colorbar	
	#if myticks is None:
	#	cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
	#else:
	#	cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8,ticks=myticks)
	#plt.clim(limits[0],limits[1])
	if len(limits) > 2:
        	plt.clim(limits[2],limits[3])
	# contour (optional)
	if type == 'spec1': 
	   z_col='k'
	else: 
	   z_col='grey'
	if type == 'labsea' or type == 'ar7zW' or type == 'ar7zE':
		#CS2 = m.contour(X, Y, tab, linewidths=0.6,levels=[500.,1500.,2500.,3000.,3100.], colors='k')
		CS2 = m.contour(X, Y, tab, linewidths=1.2,levels=[1500.,2500.,3000.], colors='b', alpha=0.5 )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=7)
        elif type == 'arc' or type == 'atln' :
                CS2 = m.contour(X, Y, tab, linewidths=0.5,levels=contours, colors='k')
                #CS2 = m.contour(X, Y, tab, linewidths=0.5,levels=contours, colors='k',alpha=0.5)
                plt.clabel(CS2, CS2.levels, inline=True, fmt='%.0f', fontsize=4)
	else:
		CS2 = m.contour(X, Y, tab, myticks, colors=z_col,linewidths=0.6 , alpha=0.7)
		#CS2 = m.contour(X, Y, tab, myticks, colors=z_col,linewidths=1.3 )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.0f', fontsize=8, colors='k')
		#plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=5, manual=True)
		#CS2 = m.contour(X, Y, tab, linewidths=0.6,levels=[0.,0.], colors='k' )
		#plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=9)
	#plt.grid(True)
	
	xspace=5. ; yspace=5.
	my_size=9
	if type == 'spec1':
	  	m.drawmeridians(npy.arange(area['lonmin'],area['lonmax']+xspace,xspace),labels=[False,False,False,True],size=my_size,color='grey',alpha=0.70)
	
	if type == 'spec1':
	  	m.drawparallels(npy.arange(area['latmin'],area['latmax']+yspace,yspace),labels=[True,False,False,False],size=my_size,color='grey',alpha=0.70)
	#plt.title(name,fontsize=10)
	return m

def nemo_Atl_MLD_plot(lon,lat,tab,tab2,contours,limits,myticks=None,name=None,zmy_cblab=None,area=None,filename='test.pdf',zmy_cmap=None,type=None,vartype=None,isocont=0,znumplt=None,zobs=0, onepp=0):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
#	rcParams['text.usetex']=True
#	rcParams['text.latex.unicode']=True
	rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=PyRaf.limit_range(tab,limits[0],limits[1])
	if type == 'spec1':
	        if vartype == 'DENS':
		    # To plot the isopycne depth over GIN seas erna_SigiDepth.py script
		    m = Basemap(width=3000000,height=2500000,lat_1=60.,lat_2=90,lon_0=-5,lat_0=70,projection='aea',resolution='c')
	        else:
		    m = Basemap(width=2300000,height=2000000,lat_1=50.,lat_2=70,lon_0=-45,lat_0=60,projection='aea',resolution='c')
		    #m = Basemap(width=4000000,height=3500000,lat_1=50.,lat_2=80,lon_0=-30,lat_0=65,projection='aea',resolution='c')
		    type=None
	elif type == 'labsea':
		# To plot the isopycne depth over Labrador sea erna_SigiDepth.py script
		m = Basemap(width=1500000,height=1600000,lat_1=50.,lat_2=65,lon_0=-50,lat_0=59.5,projection='aea',resolution='i')
	elif type == 'ar7zW':
		# To plot the isopycne depth over Labrador sea erna_SigiDepth.py script
		m = Basemap(width=400000,height=500000,lat_1=54.,lat_2=56,lon_0=-54,lat_0=55.0,projection='aea',resolution='c')
	elif type == 'ar7zE':
		# To plot the isopycne depth over Labrador sea erna_SigiDepth.py script
		m = Basemap(width=400000,height=600000,lat_1=59.,lat_2=62,lon_0=-49,lat_0=61.0,projection='aea',resolution='c')
	elif type == 'atln':
		# To plot the North Atlantic 
		m = Basemap(width=6000000,height=5500000,lat_1=25.,lat_2=75.,\
                            lon_0=-50,lat_0=50,projection='aea',resolution='i')
	else:
		m = Basemap(projection='cyl',llcrnrlat=area['latmin'],urcrnrlat=area['latmax'],\
		            llcrnrlon=area['lonmin'],urcrnrlon=area['lonmax'],resolution='c')

	#m.drawcoastlines(linewidth=0.25)
	m.fillcontinents(color='grey',lake_color='white')

	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])

	if zmy_cmap != None :
                pal = zmy_cmap
        else:
		if vartype == 'VMOD' or  vartype == 'SSH' or vartype == 'SST' or vartype == 'SSS' or vartype == 'DENS': pal = plt.cm.get_cmap('RdBu_r')
		if vartype == 'EKE' or vartype == 'MLD': pal = plt.cm.get_cmap('Blues')
		if vartype == 'EKE' or vartype == 'MLD': pal = plt.cm.get_cmap('YlOrRd')
		if vartype == 'PSI': pal = plt.cm.get_cmap('YlOrRd_r')
		pal = plt.cm.get_cmap('RdBu_r')
		# Usefull to represent surface heat flux or density flux
		#pal = plt.cm.get_cmap('RdYlBu_r') # color used for heat budget maps
		#pal = plt.cm.get_cmap('jet_r')
		#pal = plt.cm.get_cmap('seismic')
		#if vartype == 'EKE' or vartype == 'MLD': pal = plt.cm.get_cmap('BuPu')
	X,Y = m(lon,lat)
	#if zobs == 0: C = m.contourf(X,Y,tab2,[9999.,9999.],colors='grey')
	#if isocont == 0: C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm,extend='both')
	C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm,extend='both')

	############################################################################################################
	############################################################################################################
	boxtoplot=0
	if boxtoplot == 1:
		lon  = PyRaf.readfull('/gpfs5r/workgpfs/rech/ost/rost832/PRE-POST/WORK/PSI/1_FER-BCTFE04_y1990_PSI.nc','nav_lon')
		lat  = PyRaf.readfull('/gpfs5r/workgpfs/rech/ost/rost832/PRE-POST/WORK/PSI/1_FER-BCTFE04_y1990_PSI.nc','nav_lat')
        	bx_LABS={'name':'Labrador box','imin':466,'imax':472,'jmin':388,'jmax':399}
		bx_IRMw={'name':'West Irminger box','imin':495,'imax':501,'jmin':392,'jmax':400}
		bx_IRMe={'name':'East Irminger box','imin':510,'imax':517,'jmin':395,'jmax':402}
		bx_GINS={'name':'GIN box','imin':548,'imax':553,'jmin':444,'jmax':455}
		# Boxes in the 1st zoom to compute Psi
		bx_LABC={'name':'Labrador box','imin':232,'imax':249,'jmin':328,'jmax':371,'depthlim':(-3.5,0.)}
		bx_WGCB={'name':'Greenland box','imin':265,'imax':288,'jmin':391,'jmax':424,'depthlim':(-3.5,0.)}
		bx_FARB={'name':'Farewell box','imin':298,'imax':321,'jmin':359,'jmax':403,'depthlim':(-3.5,0.)}
		All_box=[bx_LABC,bx_FARB,bx_WGCB]
		#All_box=[bx_LABS,bx_IRMw,bx_IRMe,bx_GINS]
		for box in All_box:
        		lats = [lat[box['jmin'],box['imin']], lat[box['jmin'],box['imax']], lat[box['jmax'],box['imax']], lat[box['jmax'],box['imin']] ]
        		lons = [lon[box['jmin'],box['imin']], lon[box['jmin'],box['imax']], lon[box['jmax'],box['imax']], lon[box['jmax'],box['imin']] ]
        		x,y = m(lons,lats)
        		# plot filled circles at the locations of the cities.
			if box['name'] == 'West Irminger box' or box['name'] == 'East Irminger box':
				m.plot(x,y,'go')
			else:
        			m.plot(x,y,linewidth=2,color='r')
				# pour la section west de la box
        			lats = [lat[box['jmin'],box['imin']], lat[box['jmax'],box['imin']] ]
        			lons = [lon[box['jmin'],box['imin']], lon[box['jmax'],box['imin']] ]
        			x,y = m(lons,lats)
        			m.plot(x,y,linewidth=2,color='r')
	############################################################################################################
	############################################################################################################
	############################################################################################################
	############################################################################################################
	secplot=0
	if secplot == 1 :
		bx_DSO7={'name':'DSO','lon_min':-29. ,'lon_max':-22,'lat_min':69.0,'lat_max':65.5}
		bx_AR7S={'name':'AR7','lon_min':-55.5,'lon_max':-48.,'lat_min':53.4,'lat_max':61. }
		All_box=[bx_AR7S]
		for box in All_box:
        		lats = [box['lat_min'],box['lat_max']]
        		lons = [box['lon_min'],box['lon_max']]
        		x,y = m(lons,lats)
        		m.plot(x,y,linewidth=1.5, color='green')
	############################################################################################################


	if isocont == 0: 
		plt.clim(limits[0],limits[1])
	# contour (optional)
	#if len(contours) > 0 :
	# The following line to be used with MLD only
	if len(contours) > 0 and isocont == 1:
	#if len(contours) > 0 :
	   	if vartype == 'voeke': 
		   my_fsize=6
		   mycol='k'
		   mywid1=0.7
		   mywid2=0.7
		   myform='%.1f'
		   local_cont=contours
		else:
		   my_fsize=6
		   if isocont == 1: 
			 mycol='k'
			 mywid1=0.7
			 mywid2=0.7
			 myform='%.1f'
		   else:
			 mycol='k'
			 mywid1=0.3
			 mywid2=0.6
			 myform='%.2f'
			 #myform='%.1f'
		   local_cont=npy.arange(limits[0],limits[1],(contours[1]-contours[0])*2.)
		   #local_cont=contours
		if vartype == 'sobarstf' or vartype == 'ssh' or vartype == 'voeke' :
			CS2 = m.contour(X, Y, tab, local_cont, linewidths=mywid2, colors=mycol )
			#CS2 = m.contour(X, Y, tab, linewidths=mywid2,levels=contours, colors='k' )
			plt.clabel(CS2, CS2.levels, inline=True, fmt=myform, fontsize=8)
			#CS3 = m.contour(X, Y, tab, linewidths=1.5,levels=[-35,-35.], colors='r' )
			#plt.clabel(CS3, CS3.levels, inline=True, fmt=myform, fontsize=10)
		else:
			CS2 = m.contour(X, Y, tab, local_cont, colors=mycol,linewidths=mywid1 )
			plt.clabel(CS2, CS2.levels, inline=True, fmt=myform, fontsize=my_fsize)
			#CS2 = m.contour(X, Y, tab, linewidths=mywid2,levels=[0.,0.], colors=mycol )
			#plt.clabel(CS2, CS2.levels, inline=True, fmt=myform, fontsize=my_fsize)
	plt.grid(True)
	myform='%.1f' ; my_fsize=8
	#CS2 = m.contour(X, Y, tab, linewidths=0.7,levels=[3000.,3000.], colors='r' )
	#CS2 = m.contour(X, Y, tab, linewidths=0.7,levels=[1200.,1200.], colors='r' )
	#plt.clabel(CS2, CS2.levels, inline=True, fmt=myform, fontsize=my_fsize)
	#
	# x axis
	my_fsize=8
	int_lon=10.
	m.drawmeridians(npy.arange(-180.,181.,int_lon),labels=[False,False,False,True],size=my_fsize,color='grey',alpha=0.70)

	#
	# y axis
	my_fsize=8
	if type == 'labsea':
		int_lat=2.
	else:
		int_lat=5.
	m.drawparallels(npy.arange(-80.,81.,int_lat),labels=[True,False,False,False],size=my_fsize,color='grey',alpha=0.70)

	# title
	plt.title(name,fontsize=10)

	# colorbar	
	if ( vartype == 'mldr10_1' or vartype == 'votemper' or vartype == 'vosaline' or vartype == 'sobarstf' or vartype == 'ssh' or vartype == 'voeke' ) and znumplt == 'last' :
		if myticks is None:
			cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
		else:
		   	if type == 'spec1':
				posx=0.67  ; posy=0.25
			elif type == 'spec2':
				posx=0.89  ; posy=0.25
			elif type == 'atln' or type == 'labsea' :
				posx=0.98  ; posy=0.25
			else:
				if onepp == 1 :
					posx=0.20  ; posy=0.05
					cax=plt.axes([posx,posy,0.65,0.02])
				else:
					posx=0.80  ; posy=0.25
					cax=plt.axes([posx,posy,0.02,0.55])
			#cax=plt.axes([offx,0.15,offx+0.01,offy+0.5])
			if onepp == 1 :
				cb=plt.colorbar(C,cax,format='%.1f',orientation='horizontal',shrink=0.8,drawedges=True)
			else:
				#cb=plt.colorbar(C,cax,format='%.1f',orientation='vertical',shrink=0.8,drawedges=True)
				cb=plt.colorbar(C,format='%.1f',orientation='vertical',shrink=0.8,drawedges=True)
				cb.set_label(zmy_cblab,fontsize=9)
				cl = plt.getp(cb.ax, 'ymajorticklabels')
				plt.setp(cl, fontsize=9)
			#cb.ax.set_yticklabels(myticks,fontsize=6)
	#plt.savefig(filename)
	#plt.show()
	return m

def nemo_Atl_SSTS_plot(lon,lat,tab,tab2,contours,limits,myticks=None,name=None,area=None,filename='test.pdf',type=None,vartype=None,isocont=0,znumplt=None,zuniq=0,tab_ref=None,tab_sec=None,zobs=0):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
	rcParams['text.usetex']=True
#	rcParams['text.latex.unicode']=True
	#rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=PyRaf.limit_range(tab,limits[0],limits[1])
	if type == 'spec1':
		# this line for Arctic m = Basemap(width=4000000,height=3500000,lat_1=65.,lat_2=90,lon_0=0.,lat_0=90,projection='aea',resolution='c')
		if vartype == 'sobotsig0':
			m = Basemap(width=2000000,height=1150000,lat_1=50.,lat_2=70,lon_0=-20,lat_0=65,projection='aea',resolution='c')
		else:
			m = Basemap(width=4000000,height=3500000,lat_1=50.,lat_2=90,lon_0=-30,lat_0=70,projection='aea',resolution='c')
		#m = Basemap(width=4000000,height=3500000,lat_1=50.,lat_2=80,lon_0=-30,lat_0=65,projection='aea',resolution='c')
	else:
		m = Basemap(projection='cyl',llcrnrlat=area['latmin'],urcrnrlat=area['latmax'],\
		            llcrnrlon=area['lonmin'],urcrnrlon=area['lonmax'],resolution='c')
	m.drawcoastlines(linewidth=0.25)
	m.fillcontinents(color='grey',lake_color='white')
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	if len(limits) > 2:
		norm = mpl.colors.Normalize(vmin=limits[2], vmax=limits[3])
	#pal = plt.cm.get_cmap('jet')
	#pal = plt.cm.get_cmap('seismic')
	if vartype == 'VMOD' or  vartype == 'SSH' or vartype == 'SST' or vartype == 'SSS' or \
	      vartype == 'icethick' : pal = plt.cm.get_cmap('RdBu_r')
	      # for Psi only vartype == 'icethick' : pal = plt.cm.get_cmap('spectral_r')
	#pal = plt.cm.get_cmap('BuPu')
	if vartype == 'EKE' or vartype == 'MLD': pal = plt.cm.get_cmap('BuPu')
	if vartype == 'sobotsig0' : pal = plt.cm.get_cmap('BuPu')
	if vartype == 'sobotsig0' : pal = plt.cm.get_cmap('seismic')
	X,Y = m(lon,lat)
	if zobs == 0 : C = m.contourf(X,Y,tab2,[9999.,9999.],colors='grey')
	C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm,extend='both')


	############################################################################################################
	############################################################################################################
	boxtoplot=1
	if boxtoplot == 1 :
	#if boxtoplot == 1 and znumplt == 'last':
		lon  = PyRaf.readfull('/Users/ctalandi/TOOLS/My_PYTHON/DATA-EXT/ERNA/ERNA-BCTGSP3-MEAN/1990/ERNA-BCTGSP3_y1990_somxl030.nc','nav_lon')
		lat  = PyRaf.readfull('/Users/ctalandi/TOOLS/My_PYTHON/DATA-EXT/ERNA/ERNA-BCTGSP3-MEAN/1990/ERNA-BCTGSP3_y1990_somxl030.nc','nav_lat')
		bx_LABS={'name':'Labrador box','imin':468,'imax':477,'jmin':386,'jmax':399}
		bx_IRMw={'name':'Irminger box','imin':497,'imax':503,'jmin':392,'jmax':403}
		bx_IRMe={'name':'Irminger box','imin':510,'imax':515,'jmin':395,'jmax':405}
		bx_GINS={'name':'GIN box','imin':548,'imax':553,'jmin':444,'jmax':455}
		bx_FLMC={'name':'Flemish box','imin':487,'imax':496,'jmin':351,'jmax':365}
		#All_box=[bx_LABS,bx_IRMw,bx_IRMe,bx_GINS]
		All_box=[bx_LABS,bx_IRMw,bx_IRMe]
		for box in All_box:
        		lats = [lat[box['jmin'],box['imin']], lat[box['jmin'],box['imax']], lat[box['jmax'],box['imax']], lat[box['jmax'],box['imin']], lat[box['jmin'],box['imin']]]
        		lons = [lon[box['jmin'],box['imin']], lon[box['jmin'],box['imax']], lon[box['jmax'],box['imax']], lon[box['jmax'],box['imin']], lon[box['jmin'],box['imin']]]
        		x,y = m(lons,lats)
        		# plot filled circles at the locations of the cities.
			if False:
			#if box['name'] == 'West Irminger box' or box['name'] == 'East Irminger box':
				m.plot(x,y,'go')
			else:
        			m.plot(x,y,linewidth=2, color='k')
        			#m.plot(x,y,linewidth=2, color='g', markerfacecolor='b')
        			#m.plot(x,y,'ro')
	############################################################################################################
	############################################################################################################
	############################################################################################################
	############################################################################################################
	zoomtoplot=0
	if zoomtoplot == 1 :
		z_lon  = PyRaf.readfull('/Users/ctalandi/TOOLS/My_PYTHON/DATA-EXT/ERNA/GRID/ERNA-BCTGSP3/1_ERNA-BCTGSP3_mesh_hgr.nc','glamt')
		z_lat  = PyRaf.readfull('/Users/ctalandi/TOOLS/My_PYTHON/DATA-EXT/ERNA/GRID/ERNA-BCTGSP3/1_ERNA-BCTGSP3_mesh_hgr.nc','gphit')
		lons = z_lon[0,:]   ; lats = z_lat[0,:]
        	x,y = m(lons,lats)
        	m.plot(x,y,linewidth=2, color='k')
		lons = z_lon[:,0]   ; lats = z_lat[:,0]
        	x,y = m(lons,lats)
        	m.plot(x,y,linewidth=2, color='k')
		lons = z_lon[z_lon.shape[0]-1,:]   ; lats = z_lat[z_lon.shape[0]-1,:]
        	x,y = m(lons,lats)
        	m.plot(x,y,linewidth=2, color='k')
		lons = z_lon[:,z_lon.shape[1]-1]   ; lats = z_lat[:,z_lon.shape[1]-1]
        	x,y = m(lons,lats)
        	m.plot(x,y,linewidth=2, color='k')
	############################################################################################################
	############################################################################################################

	############################################################################################################
	############################################################################################################
	secplot=1
	if secplot == 1 :
		bx_AR7S={'name':'AR7','lon_min':-55.5,'lon_max':-48.,'lat_min':53.4,'lat_max':61. }
		bx_FRAM1={'name':'NON','lon_min':-19. ,'lon_max':-0.5,'lat_min':77.5,'lat_max':77.0}
		bx_FRAM2={'name':'NON','lon_min':-0.5 ,'lon_max':15.,'lat_min':77. ,'lat_max':78.0}
		# # SEC LIW
		# ybox7 = (/ 40.0, 37.5, 37.5, 40., 40./)
		# xbox7 = (/ -70., -68., -68., -70., -70./)

		All_box=[bx_AR7S]
		#All_box=[bx_FRAM1,bx_FRAM2]
		#All_box=[bx_FISS]
		for box in All_box:
        		lats = [box['lat_min'],box['lat_max']]
        		lons = [box['lon_min'],box['lon_max']]
        		x,y = m(lons,lats)
        		m.plot(x,y,linewidth=2, color='r')
	############################################################################################################
	############################################################################################################


	#C = m.contourf(lon,lat,tab,25,cmap=pal,norm=norm)
	plt.clim(limits[0],limits[1])
	if len(limits) > 2:
        	plt.clim(limits[2],limits[3])
	# contour (optional)
	if len(contours) > 0 and isocont == 1:
	   	if vartype == 'EKE': 
		   my_fsize=7
		   mycol='k'
		   mywid1=0.4
		   mywid2=0.7
		else:
		   my_fsize=8
		   #my_fsize=6
		   if isocont == 1: 
			 mycol='k'
			 mywid1=0.7
			 mywid2=0.7
		   else:
			 mycol='k'
			 mywid1=0.3
			 mywid2=0.6
		if vartype == 'SSH': 
		   contours=npy.arange(limits[0],limits[1]+1.,1.)
		#CS2 = m.contour(X, Y, tab, contours, colors=mycol,linewidths=mywid1 )
		#plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=my_fsize)
		CS2 = m.contour(X, Y, tab, linewidths=mywid2,levels=[0.,0.], colors=mycol )
		#plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=my_fsize)
	plt.grid(True)
	###   if vartype == 'icethick' or vartype == 'SST' or vartype == 'SSS':
	###   	m.contour(X, Y, tab_ref, linewidths=1.,levels=[0.05,0.05], colors='m' )
	###   	m.contour(X, Y, tab_sec, linewidths=1.,levels=[0.05,0.05], colors='g' )
	if vartype == 'sobotsig0':
		CS4=m.contour(X, Y, tab    , linewidths=1.,levels=[27.8,27.8], colors='y' )
		plt.clabel(CS4, [27.8], inline=True, fmt='%.1f', fontsize=9)
		CS4=m.contour(X, Y, tab    , linewidths=1.,levels=[27.88,27.88], colors='b' )
		plt.clabel(CS4, [27.88], inline=True, fmt='%.1f', fontsize=9)
	#
	# x axis
	locs, labels = plt.xticks()
	if type == 'spec1' or vartype =='SSH' or vartype == 'EKE' or vartype == 'SST' or vartype == 'SSS' or vartype == 'MLD' or vartype == 'sobotsig0':
	   	if vartype == 'EKE': 
		   my_fsize=9
		else:
		   my_fsize=6
	   	if vartype == 'MLD': 
		   int_lon=10.
		else:
		   int_lon=10.
		m.drawmeridians(npy.arange(-180.,181.,int_lon),labels=[False,False,False,True],size=my_fsize,color='grey',alpha=0.70)
     	else:
        	lon_grid=npy.arange(area['lonmin'],area['lonmax']+5.,5., dtype=float)
		newlocs   = npy.array(lon_grid,'f')
		newlabels = npy.array(lon_grid,'i')
		plt.xticks(newlocs,newlabels,size=8)
		plt.xlabel('Longitude',fontsize=8)
	#
	# y axis
	locsY,labelsy = plt.yticks()
	if type == 'spec1' or vartype =='SSH' or vartype == 'EKE' or vartype == 'SST' or vartype == 'SSS' or vartype == 'MLD' or vartype == 'sobotsig0':
	   	if vartype == 'EKE': 
		   my_fsize=9
		else:
		   my_fsize=6
		if vartype == 'SST':
		    m.drawparallels(npy.arange(-80.,81.,5.),labels=[True,False,False,False],size=my_fsize,color='grey',alpha=0.70)
		if vartype == 'SSS':
		    m.drawparallels(npy.arange(-80.,81.,5.),labels=[False,True,False,False],size=my_fsize,color='grey',alpha=0.70)
		if vartype == 'icethick':
		    m.drawparallels(npy.arange(-80.,81.,5.),labels=[False,True,False,False],size=my_fsize,color='grey',alpha=0.70)
		if vartype == 'sobotsig0':
		    m.drawparallels(npy.arange(-80.,81.,5.),labels=[False,True,False,False],size=my_fsize,color='grey',alpha=0.70)
		if vartype == 'SSH':
		    m.drawparallels(npy.arange(-80.,81.,5.),labels=[False,True,False,False],size=my_fsize,color='grey',alpha=0.70)
     	else:
        	lat_grid=npy.arange(area['latmin'],area['latmax']+2.,2., dtype=float)
        	#lat_grid=npy.arange(area['latmin'],area['latmax']+5.,5., dtype=float)
		newlocsy   = npy.array(lat_grid,'f')
		newlabelsy = npy.array(lat_grid,'i')
		plt.yticks(newlocsy,newlabelsy,size=8)
		plt.ylabel('Latitude',fontsize=8)	
	# title
	plt.title(name,fontsize=10)

	# colorbar	
	if vartype == 'SST' or vartype == 'SSS' or vartype == 'SSH' or vartype == 'icethick' and znumplt == 'last' or vartype == 'sobotsig0':
		if myticks is None:
			cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
		else:
		   	if zuniq == 1 :
			   if type=='spec1':
		        	posx=0.35  ; posy=0.05
			   else:
		        	posx=0.35  ; posy=0.10
			else:
		        	#if vartype == 'SST':  posx=0.35  ; posy=0.05
		        	if vartype == 'SST':  posx=0.15  ; posy=0.50
				#if vartype == 'SSS':  posx=0.35  ; posy=0.05
				if vartype == 'SSS':  posx=0.56  ; posy=0.50
				#posx=0.35  ; posy=0.05
			cax=plt.axes([posx,posy,0.33,0.02])
			cb=plt.colorbar(C,cax,format='%.2f',orientation='horizontal')
			cl = plt.getp(cb.ax, 'xmajorticklabels')
			plt.setp(cl, fontsize=8)
			#cb=plt.colorbar(C,cax,format='%.1f',orientation='vertical',shrink=0.8,ticks=myticks)
		        if vartype == 'SST':  cb.set_label(r'SST ($^\circ$C)',fontsize=7)
		        if vartype == 'SSS':  cb.set_label(r'SSS (PSU)',fontsize=7)
			cb.ax.set_yticklabels(myticks,fontsize=7)
			#cbar = plt.colorbar(C,orientation='vertical',shrink=0.8,ticks=myticks)
			#cbar = plt.colorbar(C,format='%.1f',orientation='vertical',shrink=0.8,ticks=myticks)
	#plt.savefig(filename)
	#plt.show()
	#return fig

def nemo_Atl_Rdef_plot(lon,lat,tab,contours,limits,myticks=None,name=None,area=None,filename='test.pdf',type=None,vartype=None,isocont=0,znumplt=None):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
	rcParams['text.usetex']=True
#	rcParams['text.latex.unicode']=True
	#rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=PyRaf.limit_range(tab,limits[0],limits[1])
	if type == 'spec1':
		m = Basemap(width=4000000,height=3500000,lat_1=50.,lat_2=80,lon_0=-30,lat_0=65,projection='aea',resolution='c')
		#m = Basemap(width=4000000,height=3500000,lat_1=50.,lat_2=80,lon_0=-30,lat_0=65,projection='aea',resolution='c')
	else:
		m = Basemap(projection='cyl',llcrnrlat=area['latmin'],urcrnrlat=area['latmax'],\
		            llcrnrlon=area['lonmin'],urcrnrlon=area['lonmax'],resolution='c')
	m.drawcoastlines(linewidth=0.25)
	m.fillcontinents(color='grey',lake_color='white')
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	if len(limits) > 2:
		norm = mpl.colors.Normalize(vmin=limits[2], vmax=limits[3])
	#pal = plt.cm.get_cmap('jet')
	#pal = plt.cm.get_cmap('seismic')
	#if vartype == 'Rdef' : pal = plt.cm.get_cmap('cool')
	if vartype == 'Rdef' : pal = plt.cm.get_cmap('jet')
	X,Y = m(lon,lat)
	C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm,extend='both')


	############################################################################################################
	############################################################################################################
	boxtoplot=0
	if boxtoplot == 1:
		lon  = PyRaf.readfull('/Users/ctalandi/TOOLS/My_PYTHON/DATA-EXT/ERNA/ERNA-BCTGSP3-MEAN/1990/ERNA-BCTGSP3_y1990_somxl030.nc','nav_lon')
		lat  = PyRaf.readfull('/Users/ctalandi/TOOLS/My_PYTHON/DATA-EXT/ERNA/ERNA-BCTGSP3-MEAN/1990/ERNA-BCTGSP3_y1990_somxl030.nc','nav_lat')
        	bx_LABS={'name':'Labrador box','imin':466,'imax':472,'jmin':388,'jmax':399}
		bx_IRMw={'name':'West Irminger box','imin':495,'imax':501,'jmin':392,'jmax':400}
		bx_IRMe={'name':'East Irminger box','imin':510,'imax':517,'jmin':395,'jmax':402}
		bx_GINS={'name':'GIN box','imin':548,'imax':553,'jmin':444,'jmax':455}
		All_box=[bx_LABS,bx_IRMw,bx_IRMe,bx_GINS]
		for box in All_box:
        		lats = [lat[box['jmin'],box['imin']], lat[box['jmin'],box['imax']], lat[box['jmax'],box['imax']], lat[box['jmax'],box['imin']] ]
        		lons = [lon[box['jmin'],box['imin']], lon[box['jmin'],box['imax']], lon[box['jmax'],box['imax']], lon[box['jmax'],box['imin']] ]
        		x,y = m(lons,lats)
        		# plot filled circles at the locations of the cities.
			if box['name'] == 'West Irminger box' or box['name'] == 'East Irminger box':
				m.plot(x,y,'go')
			else:
        			m.plot(x,y,'ro')
	############################################################################################################
	############################################################################################################


	#C = m.contourf(lon,lat,tab,25,cmap=pal,norm=norm)
	plt.clim(limits[0],limits[1])
	if len(limits) > 2:
        	plt.clim(limits[2],limits[3])
	# contour (optional)
	if len(contours) > 0 and isocont == 1:
		my_fsize=8
		mycol='k'
		mywid1=0.4
		mywid2=0.7
		CS2 = m.contour(X, Y, tab, contours, colors=mycol,linewidths=mywid1 )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=my_fsize)
		CS2 = m.contour(X, Y, tab, linewidths=mywid2,levels=[0.,0.], colors=mycol )
		#plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=my_fsize)
	plt.grid(True)
	#
	# x axis
	locs, labels = plt.xticks()
	my_fsize=6   ; int_lon=10.
	m.drawmeridians(npy.arange(-180.,181.,int_lon),labels=[False,False,False,True],size=my_fsize,color='grey',alpha=0.70)
	#
	# y axis
	locsY,labelsy = plt.yticks()
	my_fsize=6
	m.drawparallels(npy.arange(-80.,81.,5.),labels=[True,False,False,False],size=my_fsize,color='grey',alpha=0.70)
	# title
	plt.title(name,fontsize=10)

	# colorbar	
	if vartype == 'Rdef' and znumplt == 'last':
		if myticks is None:
			cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
		else:
			posx=0.70  ; posy=0.25
			cax=plt.axes([posx,posy,0.02,0.55])
			if isocont == 1:
				cb=plt.colorbar(C,cax,format='%.1f',orientation='vertical')
			else:
			        vmin=limits[0] ; vmax=limits[1] ; vint=10.
				cb=plt.colorbar(C,cax,format='%.1f',orientation='vertical',ticks=npy.arange(vmin,vmax+vint,vint))
				#cb=plt.colorbar(C,cax,format='%.1f',orientation='vertical',ticks=myticks)
			cb.set_label(r'Rdef (km)',fontsize=9)
			cl = plt.getp(cb.ax, 'ymajorticklabels')
			plt.setp(cl, fontsize=9)
	#plt.savefig(filename)
	#plt.show()
	#return fig

def nemo_Atl_VMOD_plot(lon,lat,tab,tabU,tabV,tabUup,tabVup,contours,limits,myticks=None,name=None,area=None,filename='test.pdf',type=None,vartype=None,isocont=0,znumplt=None,pvect=0,znormref=None,zscal=None,zloc=None):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
	rcParams['text.usetex']=True
#	rcParams['text.latex.unicode']=True
	#rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=PyRaf.limit_range(tab,limits[0],limits[1])
	if type == 'spec1':
		m = Basemap(width=4000000,height=3500000,lat_1=50.,lat_2=80,lon_0=-30,lat_0=65,projection='aea',resolution='c')
		#m = Basemap(width=4000000,height=3500000,lat_1=50.,lat_2=80,lon_0=-30,lat_0=65,projection='aea',resolution='c')
	else:
		m = Basemap(projection='cyl',llcrnrlat=area['latmin'],urcrnrlat=area['latmax'],\
		            llcrnrlon=area['lonmin'],urcrnrlon=area['lonmax'],resolution='c')
	m.drawcoastlines(linewidth=0.25)
	m.fillcontinents(color='grey',lake_color='white')
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	if len(limits) > 2:
		norm = mpl.colors.Normalize(vmin=limits[2], vmax=limits[3])
	#pal = plt.cm.get_cmap('jet')
	#pal = plt.cm.get_cmap('seismic')
	if vartype == 'VMOD' : pal = plt.cm.get_cmap('YlOrBr')
	#if vartype == 'VMOD' : pal = plt.cm.get_cmap('jet')
	#if vartype == 'VMOD' or  vartype == 'SSH' or vartype == 'SST' or vartype == 'SSS': pal = plt.cm.get_cmap('RdBu_r')
	#pal = plt.cm.get_cmap('BuPu')
	if vartype == 'EKE' or vartype == 'MLD': pal = plt.cm.get_cmap('BuPu')
	X,Y = m(lon,lat)
	if isocont == 0: C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm,extend='both')
	if pvect == 1:
	    #uproj,vproj,xx,yy = m.transform_vector(ugrid,vgrid,newlons,latitudes,31,31,returnxy=True,masked=True)
	    slon=2   ; slat=2
	    qu=m.quiver(lon[::slon,::slat], lat[::slon,::slat], tabU[::slon,::slat], tabV[::slon,::slat], color='k', angles='xy',scale=zscal,edgecolors=('k'),linewidths=(0.1,))
	    qk=plt.quiverkey(qu,-45.4,62.,znormref,'',coordinates='data',labelsep=0.05,fontproperties={'size':8})

	    qu=m.quiver(lon[::slon,::slat], lat[::slon,::slat], tabUup[::slon,::slat], tabVup[::slon,::slat], color='g', angles='xy',scale=zscal,edgecolors=('g'),linewidths=(0.1,))
	    qk=plt.quiverkey(qu,-45.4,62.7,znormref,str(znormref)+r' cm.s$^{-1}$',coordinates='data',labelsep=0.05,fontproperties={'size':8})

	############################################################################################################
	############################################################################################################
	boxtoplot=0
	if boxtoplot == 1:
		lon  = PyRaf.readfull('/Users/ctalandi/TOOLS/My_PYTHON/DATA-EXT/ERNA/ERNA-BCTGSP3-MEAN/1990/ERNA-BCTGSP3_y1990_somxl030.nc','nav_lon')
		lat  = PyRaf.readfull('/Users/ctalandi/TOOLS/My_PYTHON/DATA-EXT/ERNA/ERNA-BCTGSP3-MEAN/1990/ERNA-BCTGSP3_y1990_somxl030.nc','nav_lat')
        	bx_LABS={'name':'Labrador box','imin':466,'imax':472,'jmin':388,'jmax':399}
		bx_IRMw={'name':'West Irminger box','imin':495,'imax':501,'jmin':392,'jmax':400}
		bx_IRMe={'name':'East Irminger box','imin':510,'imax':517,'jmin':395,'jmax':402}
		bx_GINS={'name':'GIN box','imin':548,'imax':553,'jmin':444,'jmax':455}
		All_box=[bx_LABS,bx_IRMw,bx_IRMe,bx_GINS]
		for box in All_box:
        		lats = [lat[box['jmin'],box['imin']], lat[box['jmin'],box['imax']], lat[box['jmax'],box['imax']], lat[box['jmax'],box['imin']] ]
        		lons = [lon[box['jmin'],box['imin']], lon[box['jmin'],box['imax']], lon[box['jmax'],box['imax']], lon[box['jmax'],box['imin']] ]
        		x,y = m(lons,lats)
        		# plot filled circles at the locations of the cities.
			if box['name'] == 'West Irminger box' or box['name'] == 'East Irminger box':
				m.plot(x,y,'go')
			else:
        			m.plot(x,y,'ro')
	############################################################################################################
	############################################################################################################


	#C = m.contourf(lon,lat,tab,25,cmap=pal,norm=norm)
	plt.clim(limits[0],limits[1])
	if len(limits) > 2:
        	plt.clim(limits[2],limits[3])
	# contour (optional)
	if len(contours) > 0 and vartype != 'VMOD':
	#if len(contours) > 0 and isocont ==1:
	   	if vartype == 'EKE': 
		   my_fsize=7
		   mycol='k'
		   mywid1=0.4
		   mywid2=0.7
		else:
		   my_fsize=8
		   #my_fsize=6
		   if isocont == 1: 
			 mycol='k'
			 mywid1=0.7
			 mywid2=0.7
		   else:
			 mycol='k'
			 mywid1=0.3
			 mywid2=0.6
		CS2 = m.contour(X, Y, tab, contours, colors=mycol,linewidths=mywid1 )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=my_fsize)
		CS2 = m.contour(X, Y, tab, linewidths=mywid2,levels=[0.,0.], colors=mycol )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=my_fsize)
	plt.grid(True)
	#
	# x axis
	locs, labels = plt.xticks()
	if type == 'spec1' or vartype =='VMOD' or vartype == 'EKE' or vartype == 'SST' or vartype == 'SSS' or vartype == 'MLD':
	   	if vartype == 'EKE': 
		   my_fsize=9
		else:
		   my_fsize=7
	   	if vartype == 'MLD': 
		   int_lon=10.
		else:
		   int_lon=10.
		m.drawmeridians(npy.arange(-180.,181.,int_lon),labels=[False,False,False,True],size=my_fsize,color='grey',alpha=0.70)
     	else:
        	lon_grid=npy.arange(area['lonmin'],area['lonmax']+5.,5., dtype=float)
		newlocs   = npy.array(lon_grid,'f')
		newlabels = npy.array(lon_grid,'i')
		plt.xticks(newlocs,newlabels,size=8)
		plt.xlabel('Longitude',fontsize=8)
	#
	# y axis
	locsY,labelsy = plt.yticks()
	if type == 'spec1' or vartype =='VMOD' or vartype == 'EKE' or vartype == 'SST' or vartype == 'SSS' or vartype == 'MLD':
	   	if vartype == 'EKE': 
		   my_fsize=9
		else:
		   my_fsize=7
		if zloc == 'uleft' or zloc == 'lleft':
		    m.drawparallels(npy.arange(-80.,81.,5.),labels=[True,False,False,False],size=my_fsize,color='grey',alpha=0.70)
		if zloc == 'lright' or  zloc == 'uright'  :
		    m.drawparallels(npy.arange(-80.,81.,5.),labels=[False,True,False,False],size=my_fsize,color='grey',alpha=0.70)
     	else:
        	lat_grid=npy.arange(area['latmin'],area['latmax']+2.,2., dtype=float)
        	#lat_grid=npy.arange(area['latmin'],area['latmax']+5.,5., dtype=float)
		newlocsy   = npy.array(lat_grid,'f')
		newlabelsy = npy.array(lat_grid,'i')
		plt.yticks(newlocsy,newlabelsy,size=8)
		plt.ylabel('Latitude',fontsize=8)	
	# title
	plt.title(name,fontsize=12)

	# colorbar	
	if vartype == 'MLD' or vartype == 'VMOD' and znumplt == 'last':
		if myticks is None:
			cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
		else:
		   	####     old way if type == 'spec1':
			####     old way 	posx=0.67  ; posy=0.25
			####     old way elif type == 'spec2':
			####     old way 	posx=0.89  ; posy=0.25
			####     old way else:
			####     old way 	posx=0.90  ; posy=0.25
			####     old way 	#posx=0.80  ; posy=0.25
			####     old way cax=plt.axes([posx,posy,0.02,0.55])
			####     old way #cax=plt.axes([offx,0.15,offx+0.01,offy+0.5])
			####     old way if isocont == 1:
			####     old way 	cb=plt.colorbar(CS2,cax,format='%.1f',orientation='vertical')
			####     old way else:
			####     old way 	cb=plt.colorbar(C,cax,format='%.1f',orientation='vertical')
			####     old way #cb=plt.colorbar(C,cax,format='%.1f',orientation='vertical',shrink=0.8,ticks=myticks)
			####     old way cb.set_label(r'(1 x 10$^{-2}$ m s$^{-1}$)',fontsize=9)
			# horizontal
		        #if zloc == 'lright':  posx=0.14  ; posy=0.50
			#if zloc == 'uright':  posx=0.55  ; posy=0.50
			#cax=plt.axes([posx,posy,0.33,0.02])
			#cb=plt.colorbar(C,cax,format='%.2f',orientation='horizontal')
			#cl = plt.getp(cb.ax, 'xmajorticklabels')
			#plt.setp(cl, fontsize=8)
		        #cb.set_label(r'velocity module (cm s$^{-1}$)',fontsize=7)
			#cb.ax.set_yticklabels(myticks,fontsize=7)
			# vertical
		        if zloc == 'uright':  posx=0.49  ; posy=0.54
			if zloc == 'lright':  posx=0.49  ; posy=0.14
			cax=plt.axes([posx,posy,0.01,0.33])
			cb=plt.colorbar(C,cax,format='%.0f',orientation='vertical')
			cl = plt.getp(cb.ax, 'ymajorticklabels')
			plt.setp(cl, fontsize=8)
		        cb.set_label(r'velocity module (cm s$^{-1}$)',fontsize=7)
			#cb.ax.set_yticklabels(myticks,fontsize=7)

	return None

def nemo_Atl_BVS_plot(lon,lat,tab,contours,limits,myticks=None,name=None,area=None,filename='test.pdf',type=None,vartype=None,isocont=0,znumplt=None,znormref=None,zscal=None):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
	#rcParams['text.usetex']=True
#	rcParams['text.latex.unicode']=True
	rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=PyRaf.limit_range(tab,limits[0],limits[1])
	if type == 'spec1':
		m = Basemap(width=4000000,height=3500000,lat_1=50.,lat_2=80,lon_0=-30,lat_0=65,projection='aea',resolution='c')
		#m = Basemap(width=4000000,height=3500000,lat_1=50.,lat_2=80,lon_0=-30,lat_0=65,projection='aea',resolution='c')
	else:
		m = Basemap(projection='cyl',llcrnrlat=area['latmin'],urcrnrlat=area['latmax'],\
		            llcrnrlon=area['lonmin'],urcrnrlon=area['lonmax'],resolution='c')
	m.drawcoastlines(linewidth=0.25)
	m.fillcontinents(color='grey',lake_color='white')
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	if len(limits) > 2:
		norm = mpl.colors.Normalize(vmin=limits[2], vmax=limits[3])
	#if vartype == 'BVS' or  vartype == 'WDIA' : pal=plt.cm.get_cmap('RdBu_r')
	if vartype == 'BVS' or  vartype == 'WDIA' : pal=plt.cm.get_cmap('seismic')

	X,Y = m(lon,lat)
	if isocont == 1: C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm,extend='both')

	############################################################################################################
	############################################################################################################
	boxtoplot=0
	if boxtoplot == 1:
		lon  = PyRaf.readfull('/Users/ctalandi/TOOLS/My_PYTHON/DATA-EXT/ERNA/ERNA-BCTGSP3-MEAN/1990/ERNA-BCTGSP3_y1990_somxl030.nc','nav_lon')
		lat  = PyRaf.readfull('/Users/ctalandi/TOOLS/My_PYTHON/DATA-EXT/ERNA/ERNA-BCTGSP3-MEAN/1990/ERNA-BCTGSP3_y1990_somxl030.nc','nav_lat')
        	bx_LABS={'name':'Labrador box','imin':466,'imax':472,'jmin':388,'jmax':399}
		bx_IRMw={'name':'West Irminger box','imin':495,'imax':501,'jmin':392,'jmax':400}
		bx_IRMe={'name':'East Irminger box','imin':510,'imax':517,'jmin':395,'jmax':402}
		bx_GINS={'name':'GIN box','imin':548,'imax':553,'jmin':444,'jmax':455}
		All_box=[bx_LABS,bx_IRMw,bx_IRMe,bx_GINS]
		for box in All_box:
        		lats = [lat[box['jmin'],box['imin']], lat[box['jmin'],box['imax']], lat[box['jmax'],box['imax']], lat[box['jmax'],box['imin']] ]
        		lons = [lon[box['jmin'],box['imin']], lon[box['jmin'],box['imax']], lon[box['jmax'],box['imax']], lon[box['jmax'],box['imin']] ]
        		x,y = m(lons,lats)
        		# plot filled circles at the locations of the cities.
			if box['name'] == 'West Irminger box' or box['name'] == 'East Irminger box':
				m.plot(x,y,'go')
			else:
        			m.plot(x,y,'ro')
	############################################################################################################
	############################################################################################################


	#C = m.contourf(lon,lat,tab,25,cmap=pal,norm=norm)
	plt.clim(limits[0],limits[1])
	if len(limits) > 2:
        	plt.clim(limits[2],limits[3])
	# contour (optional)
	if len(contours) > 0 and vartype != 'BVS' and vartype != 'WDIA':
	#if len(contours) > 0 and isocont ==1:
	   	if vartype == 'EKE': 
		   my_fsize=7
		   mycol='k'
		   mywid1=0.4
		   mywid2=0.7
		else:
		   my_fsize=6
		   #my_fsize=6
		   if isocont == 1: 
			 mycol='k'
			 mywid1=0.7
			 mywid2=0.7
		   else:
			 mycol='k'
			 mywid1=0.3
			 mywid2=0.6
		CS2 = m.contour(X, Y, tab, contours, colors=mycol,linewidths=mywid1 )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=my_fsize)
		CS2 = m.contour(X, Y, tab, linewidths=mywid2,levels=[0.,0.], colors=mycol )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=my_fsize)
	plt.grid(True)
	#
	# x axis
	locs, labels = plt.xticks()
	if type == 'spec1' or vartype =='BVS' or vartype == 'EKE' or vartype == 'WDIA' :
	   	if vartype == 'EKE': 
		   my_fsize=9
		else:
		   my_fsize=9
	   	if vartype == 'MLD': 
		   int_lon=10.
		else:
		   int_lon=10.
		m.drawmeridians(npy.arange(-180.,181.,int_lon),labels=[False,False,False,True],size=my_fsize,color='grey',alpha=0.70)
     	else:
        	lon_grid=npy.arange(area['lonmin'],area['lonmax']+5.,5., dtype=float)
		newlocs   = npy.array(lon_grid,'f')
		newlabels = npy.array(lon_grid,'i')
		plt.xticks(newlocs,newlabels,size=8)
		plt.xlabel('Longitude',fontsize=10)
	#
	# y axis
	locsY,labelsy = plt.yticks()
	if type == 'spec1' or vartype =='BVS' or vartype == 'EKE' or vartype == 'WDIA' :
	   	if vartype == 'EKE': 
		   my_fsize=9
		else:
		   my_fsize=9
		m.drawparallels(npy.arange(-80.,81.,5.),labels=[True,False,False,False],size=my_fsize,color='grey',alpha=0.70)
     	else:
        	lat_grid=npy.arange(area['latmin'],area['latmax']+2.,2., dtype=float)
        	#lat_grid=npy.arange(area['latmin'],area['latmax']+5.,5., dtype=float)
		newlocsy   = npy.array(lat_grid,'f')
		newlabelsy = npy.array(lat_grid,'i')
		plt.yticks(newlocsy,newlabelsy,size=8)
		plt.ylabel('Latitude',fontsize=10)	
	# title
        if vartype == 'WDIA':
	   plt.title(name,fontsize=9)
        else:
	   plt.title(name,fontsize=12)

	# colorbar	
	if vartype == 'WDIA' or vartype == 'BVS' and znumplt == 'last':
		if myticks is None:
			cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
		else:
		   	if type == 'spec1':
				posx=0.67  ; posy=0.25
			elif type == 'spec2':
				posx=0.89  ; posy=0.25
			else:
                                if vartype == 'WDIA' :
				   posx=0.91  ; posy=0.25
                                else:
				   posx=0.69  ; posy=0.25
				   #posx=0.83  ; posy=0.25
			cax=plt.axes([posx,posy,0.015,0.55])
			#cax=plt.axes([offx,0.15,offx+0.01,offy+0.5])
			if isocont == 1 and vartype != 'BVS' and vartype != 'WDIA':
				cb=plt.colorbar(CS2,cax,format='%.1f',orientation='vertical')
			else:
                                if vartype == 'WDIA' :
				    cb=plt.colorbar(C,cax,format='%.0f',orientation='vertical')
                                else:
				    cb=plt.colorbar(C,cax,format='%.1f',orientation='vertical')
			#cb=plt.colorbar(C,cax,format='%.1f',orientation='vertical',shrink=0.8,ticks=myticks)
			if vartype == 'BVS': cb.set_label(r'(1 x 10$^{-6} N m$^{-3}$)',fontsize=8)
			if vartype == 'WDIA': cb.set_label(r'( m year$^{-1}$)',fontsize=8)
			#cb.ax.set_yticks(myticks,fontsize=8)
			#cb.ax.set_yticklabels(myticks,fontsize=8)
			#cbar = plt.colorbar(C,orientation='vertical',shrink=0.8,ticks=myticks)
			#cbar = plt.colorbar(C,format='%.1f',orientation='vertical',shrink=0.8,ticks=myticks)
	#plt.savefig(filename)
	#plt.show()
	return None

def nemo_Atl_GS_plot(lon,lat,tab1,tab2,tab3,contours,limits,myticks=None,name=None,area=None,filename='test.pdf',type=None,snd_grid=None,snd_lat=None,snd_lon=None,):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
	rcParams['text.usetex']=True
#	rcParams['text.latex.unicode']=True
#	rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	#tab1=PyRaf.limit_range(tab1,limits[0],limits[1])
	#tab2=PyRaf.limit_range(tab2,limits[0],limits[1])
	# background
	## test   m = Basemap(projection='ortho',lat_0=60.,lon_0=-40.,resolution='c')
	m = Basemap(projection='cyl',llcrnrlat=area['latmin'],urcrnrlat=area['latmax'],\
		llcrnrlon=area['lonmin'],urcrnrlon=area['lonmax'],resolution='c')
	m.drawcoastlines(linewidth=0.25)
	m.fillcontinents(color='grey',lake_color='white')
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	if len(limits) > 2:
		norm = mpl.colors.Normalize(vmin=limits[2], vmax=limits[3])
	#pal = PyRaf_colormaps.gen_pal_Testu2()
	#pal = PyRaf_colormaps.gen_pal_blue2red3()
	#pal = plt.cm.get_cmap('jet')
	pal = plt.cm.get_cmap('seismic')
	pal = plt.cm.get_cmap('RdBu_r')
	pal = plt.cm.get_cmap('BuPu')
	# for EKE plot pal = plt.cm.get_cmap('BuPu')
	X,Y = m(lon,lat)
	if snd_grid == 'second': snd_X,snd_Y = m(snd_lon,snd_lat)
	#C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm,extend='both')
	#plt.clim(limits[0],limits[1])
	#if len(limits) > 2:
        #	plt.clim(limits[2],limits[3])
	# contour (optional)
	CS1 = m.contour(X, Y, tab1, linewidths=0.8,levels=[17.,17.], colors='k' )
	CS2 = m.contour(X, Y, tab2, linewidths=0.8,levels=[17.,17.], colors='r' )
	#CS2 = m.contour(X, Y, tab2, linewidths=1.4,levels=[17.,17.], colors='r' , linestyles='dashdot')
	if snd_grid == 'second':
	    CS3 = m.contour(snd_X, snd_Y, tab3, linewidths=1.1,levels=[17.,17.], colors='g' )
	else:
	    CS3 = m.contour(X, Y, tab3, linewidths=1.1,levels=[17.,17.], colors='g' )
	    #CS3 = m.contour(X, Y, tab3, linewidths=1.1,levels=[17.,17.], colors='g' , linestyles='dashed')
	#plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=5)
	plt.grid(True,color='grey',alpha=0.4)
	#
	# x axis
	locs, labels = plt.xticks()
	if type == 'spec1':
		m.drawmeridians(npy.arange(-180.,181.,5.),labels=[False,False,False,True],size=6,alpha=0.7)
     	else:
        	lon_grid=npy.arange(area['lonmin'],area['lonmax']+5.,5., dtype=float)
		newlocs   = npy.array(lon_grid,'f')
		newlabels = npy.array(lon_grid,'i')
		plt.xticks(newlocs,newlabels,size=8)
		plt.xlabel('Longitude',fontsize=8)
	#
	# y axis
	locsY,labelsy = plt.yticks()
	if type == 'spec1':
		m.drawparallels(npy.arange(-80.,81.,5.),labels=[True,False,False,False],size=6,alpha=0.7)
     	else:
        	lat_grid=npy.arange(area['latmin'],area['latmax']+5.,5., dtype=float)
		newlocsy   = npy.array(lat_grid,'f')
		newlabelsy = npy.array(lat_grid,'i')
		plt.yticks(newlocsy,newlabelsy,size=8)
		plt.ylabel('Latitude',fontsize=8)	
	# title
	plt.title(name,fontsize=10)
	#plt.savefig(filename)
	#plt.show()
	#return fig

def nemo_Atl_SSH_plot(lon,lat,tab,contours,limits,myticks=None,name=None,area=None,filename='test.pdf',type=None,vartype=None,isocont=0,znumplt=None):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
	rcParams['text.usetex']=True
#	rcParams['text.latex.unicode']=True
	#rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=PyRaf.limit_range(tab,limits[0],limits[1])
	if type == 'spec1':
		m = Basemap(width=4000000,height=3500000,lat_1=50.,lat_2=80,lon_0=-30,lat_0=65,projection='aea',resolution='c')
		#m = Basemap(width=4000000,height=3500000,lat_1=50.,lat_2=80,lon_0=-30,lat_0=65,projection='aea',resolution='c')
	else:
		m = Basemap(projection='cyl',llcrnrlat=area['latmin'],urcrnrlat=area['latmax'],\
		            llcrnrlon=area['lonmin'],urcrnrlon=area['lonmax'],resolution='c')
	m.drawcoastlines(linewidth=0.25)
	m.fillcontinents(color='grey',lake_color='white')
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	if len(limits) > 2:
		norm = mpl.colors.Normalize(vmin=limits[2], vmax=limits[3])
	#pal = plt.cm.get_cmap('jet')
	#pal = plt.cm.get_cmap('seismic')
	if vartype == 'SSH' or vartype == 'SST' or vartype == 'SSS': pal = plt.cm.get_cmap('RdBu_r')
	#pal = plt.cm.get_cmap('BuPu')
	if vartype == 'EKE' or vartype == 'MLD': pal = plt.cm.get_cmap('BuPu')
	X,Y = m(lon,lat)
	C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm,extend='both')


	############################################################################################################
	############################################################################################################
	boxtoplot=0
	if boxtoplot == 1:
		lon  = PyRaf.readfull('/Users/ctalandi/TOOLS/My_PYTHON/DATA-EXT/ERNA/ERNA-BCTGSP3-MEAN/1990/ERNA-BCTGSP3_y1990_somxl030.nc','nav_lon')
		lat  = PyRaf.readfull('/Users/ctalandi/TOOLS/My_PYTHON/DATA-EXT/ERNA/ERNA-BCTGSP3-MEAN/1990/ERNA-BCTGSP3_y1990_somxl030.nc','nav_lat')
        	bx_LABS={'name':'Labrador box','imin':466,'imax':472,'jmin':388,'jmax':399}
		bx_IRMw={'name':'West Irminger box','imin':495,'imax':501,'jmin':392,'jmax':400}
		bx_IRMe={'name':'East Irminger box','imin':510,'imax':517,'jmin':395,'jmax':402}
		bx_GINS={'name':'GIN box','imin':548,'imax':553,'jmin':444,'jmax':455}
		All_box=[bx_LABS,bx_IRMw,bx_IRMe,bx_GINS]
		for box in All_box:
        		lats = [lat[box['jmin'],box['imin']], lat[box['jmin'],box['imax']], lat[box['jmax'],box['imax']], lat[box['jmax'],box['imin']] ]
        		lons = [lon[box['jmin'],box['imin']], lon[box['jmin'],box['imax']], lon[box['jmax'],box['imax']], lon[box['jmax'],box['imin']] ]
        		x,y = m(lons,lats)
        		# plot filled circles at the locations of the cities.
			if box['name'] == 'West Irminger box' or box['name'] == 'East Irminger box':
				m.plot(x,y,'go')
			else:
        			m.plot(x,y,'ro')
	############################################################################################################
	############################################################################################################


	#C = m.contourf(lon,lat,tab,25,cmap=pal,norm=norm)
	plt.clim(limits[0],limits[1])
	if len(limits) > 2:
        	plt.clim(limits[2],limits[3])
	# contour (optional)
	if len(contours) > 0 :
	   	if vartype == 'EKE': 
		   my_fsize=7
		else:
		   my_fsize=6
		   if isocont == 1: 
			 mycol='k'
			 mywid1=0.7
			 mywid2=0.7
		   else:
			 mycol='k'
			 mywid1=0.3
			 mywid2=0.7
		CS2 = m.contour(X, Y, tab, contours, colors=mycol,linewidths=mywid1 )
		#plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=my_fsize)
		CS2 = m.contour(X, Y, tab, linewidths=mywid2,levels=[0.,0.], colors=mycol )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=my_fsize)
	#
	# x axis
	locs, labels = plt.xticks()
	if type == 'spec1' or vartype =='SSH' or vartype == 'EKE' or vartype == 'SST' or vartype == 'SSS' or vartype == 'MLD':
	   	if vartype == 'EKE': 
		   my_fsize=9
		else:
		   my_fsize=6
	   	if vartype == 'MLD': 
		   int_lon=10.
		else:
		   int_lon=10.
		m.drawmeridians(npy.arange(-180.,181.,int_lon),labels=[False,False,False,True],size=my_fsize,color='grey',alpha=0.70)
     	else:
        	lon_grid=npy.arange(area['lonmin'],area['lonmax']+5.,5., dtype=float)
		newlocs   = npy.array(lon_grid,'f')
		newlabels = npy.array(lon_grid,'i')
		plt.xticks(newlocs,newlabels,size=8)
		plt.xlabel('Longitude',fontsize=8)
	#
	# y axis
	locsY,labelsy = plt.yticks()
	if type == 'spec1' or vartype =='SSH' or vartype == 'EKE' or vartype == 'SST' or vartype == 'SSS' or vartype == 'MLD':
	   	if vartype == 'EKE': 
		   my_fsize=9
		else:
		   my_fsize=6
		m.drawparallels(npy.arange(-80.,81.,5.),labels=[True,False,False,False],size=my_fsize,color='grey',alpha=0.7)
     	else:
        	lat_grid=npy.arange(area['latmin'],area['latmax']+2.,2., dtype=float)
        	#lat_grid=npy.arange(area['latmin'],area['latmax']+5.,5., dtype=float)
		newlocsy   = npy.array(lat_grid,'f')
		newlabelsy = npy.array(lat_grid,'i')
		plt.yticks(newlocsy,newlabelsy,size=8)
		plt.ylabel('Latitude',fontsize=8)	
	# title
	plt.title(name,fontsize=10)

	# colorbar	
	if vartype == 'SSH' and znumplt == 'last':
		if myticks is None:
			cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
		else:
		   	if type == 'spec2':
				posx=0.66  ; posy=0.23
			else:
				posx=0.80  ; posy=0.25
			cax=plt.axes([posx,posy,0.02,0.55])
			if isocont == 1:
				cb=plt.colorbar(CS2,cax,format='%.1f',orientation='vertical')
			else:
			        vmin=limits[0] ; vmax=limits[1] ; vint=10.
				cb=plt.colorbar(C,cax,format='%.0f',orientation='vertical',ticks=npy.arange(vmin,vmax+vint,vint))
				#cb=plt.colorbar(C,cax,format='%.1f',orientation='vertical',ticks=myticks)
			cb.set_label(r'SSH (cm)',fontsize=7)
			cl = plt.getp(cb.ax, 'ymajorticklabels')
			plt.setp(cl, fontsize=7)

			#cb.ax.set_yticklabels(myticks,fontsize=7)
			#cbar = plt.colorbar(C,orientation='vertical',shrink=0.8,ticks=myticks)
			#cbar = plt.colorbar(C,format='%.1f',orientation='vertical',shrink=0.8,ticks=myticks)

	#plt.savefig(filename)
	#plt.show()
	#return fig


def nemo_Atl_std_plot(lon,lat,tab,std_tab,contours,limits,myticks=None,name=None,filename='test.pdf'):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
#	rcParams['text.usetex']=True
#	rcParams['text.latex.unicode']=True
	rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=PyRaf.limit_range(tab,limits[0],limits[1])
	#fig = plt.figure(figsize=[12.8 , 6.])
	#ax  = fig.add_subplot(111)
	# background
	m = Basemap(projection='cyl',llcrnrlat=20,urcrnrlat=80,\
	            llcrnrlon=-90,urcrnrlon=0,resolution='c')
	#m = Basemap(projection='cyl',llcrnrlat=-60,urcrnrlat=-20,\
	#            llcrnrlon=-180,urcrnrlon=180,resolution='c')
	m.drawcoastlines(linewidth=0.25)
	m.fillcontinents(color='grey',lake_color='white')
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	if len(limits) > 2:
		norm = mpl.colors.Normalize(vmin=limits[2], vmax=limits[3])
	#pal = PyRaf_colormaps.gen_pal_Testu2()
	#pal = PyRaf_colormaps.gen_pal_blue2red3()
	#pal = plt.cm.get_cmap('jet')
	pal = plt.cm.get_cmap('seismic')
	# for EKE plot pal = plt.cm.get_cmap('BuPu')
	X,Y = m(lon,lat)
	C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm)
	#C = m.contourf(lon,lat,tab,25,cmap=pal,norm=norm)
	# colorbar	
	#if myticks is None:
	#	cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
	#else:
	#	cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8,ticks=myticks)
	plt.clim(limits[0],limits[1])
	if len(limits) > 2:
        	plt.clim(limits[2],limits[3])
	# contour (optional)
	if len(contours) > 0:
		CS2 = m.contour(X, Y, tab, contours, colors='k',linewidths=0.3 )
		#CS2 = m.contour(lon, lat, tab, contours, colors='k',linewidths=0.3 )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=6)

        # Plot the standard deviation
        CS = plt.contour(X,Y,std_tab,colors='c',linewidths=1.,levels=npy.arange(-10.,12.,1.))
        plt.clabel(CS,fontsize=10,fmt='%2.1f')

	plt.grid(True)
	#
	# x axis
	locs, labels = plt.xticks()
	newlocs   = npy.array([-90,-80,-70,-60,-50,-40,-30,-20,-10,0],'f')
	newlabels = npy.array([-90,-80,-70,-60,-50,-40,-30,-20,-10,0],'i')
	#newlocs   = npy.array([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180],'f')
	#newlabels = npy.array([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180],'i')
	plt.xticks(newlocs,newlabels,size=8)
	plt.xlabel('Longitude',fontsize=10)
	#
	# y axis
	locsY,labelsy = plt.yticks()
	newlocsy   = npy.array([20,30,40,50,60,70,80],'f')
	newlabelsy = npy.array([20,30,40,50,60,70,80],'i')
	plt.yticks(newlocsy,newlabelsy,size=8)
	plt.ylabel('Latitude',fontsize=10)	
	# title
	plt.title(name,fontsize=10)
	#plt.savefig(filename)
	#plt.show()
	#return fig

def nemo_Atl_FLOzoom_plot(lon,lat,tab,contours,limits,myticks=None,name=None,filename='test.pdf'):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
#	rcParams['text.usetex']=True
#	rcParams['text.latex.unicode']=True
	rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=PyRaf.limit_range(tab,limits[0],limits[1])
	#fig = plt.figure(figsize=[12.8 , 6.])
	#ax  = fig.add_subplot(111)
	# background
	m = Basemap(projection='cyl',llcrnrlat=20,urcrnrlat=40,\
	            llcrnrlon=-90,urcrnrlon=-70,resolution='c')
	m.drawcoastlines(linewidth=0.25)
	m.fillcontinents(color='grey',lake_color='white')
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	if len(limits) > 2:
		norm = mpl.colors.Normalize(vmin=limits[2], vmax=limits[3])
	#pal = PyRaf_colormaps.gen_pal_Testu2()
	#pal = PyRaf_colormaps.gen_pal_blue2red3()
	#pal = plt.cm.get_cmap('jet')
	pal = plt.cm.get_cmap('seismic')
	pal = plt.cm.get_cmap('RdBu_r')
	# for EKE plot pal = plt.cm.get_cmap('BuPu')
	X,Y = m(lon,lat)
	C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm)
	#C = m.contourf(lon,lat,tab,25,cmap=pal,norm=norm)
	# colorbar	
	#if myticks is None:
	#	cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
	#else:
	#	cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8,ticks=myticks)
	plt.clim(limits[0],limits[1])
	if len(limits) > 2:
        	plt.clim(limits[2],limits[3])
	# contour (optional)
	if len(contours) > 0:
		CS2 = m.contour(X, Y, tab, contours, colors='k',linewidths=0.3 )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=5)
		CS2 = m.contour(X, Y, tab, linewidths=0.6,levels=[0.,0.], colors='k' )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=5)
	plt.grid(True)
	#
	# x axis
	locs, labels = plt.xticks()
	newlocs   = npy.array([-90,-88,-86,-84,-82,-80,-78,-76,-74,-72,-70],'f')
	newlabels = npy.array([-90,-88,-86,-84,-82,-80,-78,-76,-74,-72,-70],'i')
	#newlocs   = npy.array([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180],'f')
	#newlabels = npy.array([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180],'i')
	plt.xticks(newlocs,newlabels,size=8)
	plt.xlabel('Longitude',fontsize=10)
	#
	# y axis
	locsY,labelsy = plt.yticks()
	newlocsy   = npy.array([20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40],'f')
	newlabelsy = npy.array([20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40],'i')
	plt.yticks(newlocsy,newlabelsy,size=8)
	plt.ylabel('Latitude',fontsize=10)	
	# title
	plt.title(name,fontsize=10)
	#plt.savefig(filename)
	#plt.show()
	#return fig

def nemo_midAtl_plot(lon,lat,tab,contours,limits,myticks=None,name=None,filename='test.pdf'):
	#
# Les 2 lignes suivantes posent un probleme a l'execution lie a LaTex ....
#	rcParams['text.usetex']=True
#	rcParams['text.latex.unicode']=True
	rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=PyRaf.limit_range(tab,limits[0],limits[1])
	#fig = plt.figure(figsize=[12.8 , 6.])
	#ax  = fig.add_subplot(111)
	# background
	m = Basemap(projection='cyl',llcrnrlat=40,urcrnrlat=50,\
	            llcrnrlon=-30,urcrnrlon=-8,resolution='c')
	m.drawcoastlines(linewidth=0.25)
	m.fillcontinents(color='grey',lake_color='white')
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	if len(limits) > 2:
		norm = mpl.colors.Normalize(vmin=limits[2], vmax=limits[3])
	#pal = PyRaf_colormaps.gen_pal_Testu2()
	#pal = PyRaf_colormaps.gen_pal_blue2red3()
	#pal = plt.cm.get_cmap('jet')
	pal = plt.cm.get_cmap('seismic')
	pal = plt.cm.get_cmap('RdBu_r')
	# for EKE plot pal = plt.cm.get_cmap('BuPu')
	X,Y = m(lon,lat)
	C = m.contourf(X,Y,tab,contours,cmap=pal,norm=norm)
	#C = m.contourf(lon,lat,tab,25,cmap=pal,norm=norm)
	# colorbar	
	#if myticks is None:
	#	cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8)
	#else:
	#	cbar = plt.colorbar(C,format='%.2f',orientation='vertical',shrink=0.8,ticks=myticks)
	plt.clim(limits[0],limits[1])
	if len(limits) > 2:
        	plt.clim(limits[2],limits[3])
	# contour (optional)
	if len(contours) > 0:
		CS2 = m.contour(X, Y, tab, contours, colors='grey',linewidths=0.3 )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=5)
		CS2 = m.contour(X, Y, tab, linewidths=0.6,levels=[0.,0.], colors='grey' )
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=5)
	plt.grid(True)
	#
	# x axis
	locs, labels = plt.xticks()
	newlocs   = npy.array([-30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8],'f')
	newlabels = npy.array([-30,-28,-26,-24,-22,-20,-18,-16,-14,-12,-10,-8],'i')
	#newlocs   = npy.array([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180],'f')
	#newlabels = npy.array([-180,-150,-120,-90,-60,-30,0,30,60,90,120,150,180],'i')
	plt.xticks(newlocs,newlabels,size=8)
	plt.xlabel('Longitude',fontsize=10)
	#
	# y axis
	locsY,labelsy = plt.yticks()
	newlocsy   = npy.array([40,41,42,43,44,45,46,47,48,49,50],'f')
	newlabelsy = npy.array([40,41,42,43,44,45,46,47,48,49,50],'i')
	plt.yticks(newlocsy,newlabelsy,size=8)
	plt.ylabel('Latitude',fontsize=10)	
	# title
	plt.title(name,fontsize=10)
	#plt.savefig(filename)
	#plt.show()
	#return fig

def nemo_global_diffplot(lon,lat,tab,contours,limits,myticks=None,name=None,filename='test.pdf'):
	#
	rcParams['text.usetex']=True
	rcParams['text.latex.unicode']=True
	rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=PyRaf.limit_range(tab,limits[0],limits[1])
	fig = plt.figure(figsize=[12.8 , 6.])
	ax  = fig.add_subplot(111)
	# background
	m = Basemap(projection='cyl',llcrnrlat=-89,urcrnrlat=89,\
	            llcrnrlon=0,urcrnrlon=420,resolution='c')
	m.drawcoastlines()
	m.fillcontinents(color='grey',lake_color='white')
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	if len(limits) > 2:
		norm = mpl.colors.Normalize(vmin=limits[2], vmax=limits[3])
	pal = PyRaf_colormaps.gen_pal_diff()
	C = m.contourf(lon,lat,tab,25,cmap=pal,norm=norm)
	
	if myticks is None:
		cbar = plt.colorbar(C,format='%.2f',orientation='horizontal',shrink=0.8)
	else:
		cbar = plt.colorbar(C,format='%.2f',orientation='horizontal',shrink=0.8,ticks=myticks)
	plt.clim(limits[0],limits[1])
	if len(limits) > 2:
        	plt.clim(limits[2],limits[3])
	# contour
	if len(contours) > 0:
		CS2 = m.contour(lon, lat, tab, contours, colors='k')
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=8)
	#
	# x axis
	locs, labels = plt.xticks()
	newlocs   = npy.array([0,30,60,90,120,150,180,210,240,270,300,330,360,390,420],'f')
	newlabels = npy.array([0,30,60,90,120,150,180,210,240,270,300,330,0,30,60],'i')
	plt.xticks(newlocs,newlabels)
	plt.xlabel('Longitude',fontsize=16)
	#
	# y axis
	locsY,labelsy = plt.yticks()
	newlocsy   = npy.array([-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80],'f')
	newlabelsy = npy.array([-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80],'i')
	plt.yticks(newlocsy,newlabelsy)
	plt.ylabel('Latitude',fontsize=16)	
	# title
	plt.title(name,fontsize=18)
	plt.savefig(filename)
	return fig

### more plots soon
