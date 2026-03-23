#!/usr/bin/env python

import numpy as npy
import matplotlib.pylab as plt
import matplotlib as mpl
import subprocess
from netCDF4 import Dataset
import xarray as xr 
from datetime import datetime

################################################################################################################################
def DEF_LOC_SEC( CONFIG, zsect ) :
################################################################################################################################

	#------------------------------------------------------------------------------------------------------------------------
	# Localisation of straits >>>> FORTRAN indices
	# Indices are those from ocean namelists

	if CONFIG == "CREG025.L75" :
		if zsect['name'] == "Bering"   : loc_sec={'name':"Bering"   ,'is': 202,'ie':223 ,'jloc':601 ,'depthlim':(-500.,0.) ,'templim':(-2.,3.,0.5), 'salilim':(28.,35.,0.5) , 'velolim':(-5.,5.,1.)   ,'denslim':(24.0,27.9,0.2),'tkeklim':(-7.,1.5,1.)}
		if zsect['name'] == "SouthG"   : loc_sec={'name':"SouthG"   ,'is':  54,'ie':324 ,'jloc':  2 ,'depthlim':(-5000.,0.),'templim':( 4.,27.,2.), 'salilim':(33.,36.,0.2) , 'velolim':(-50.,50.,10.),'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
		if zsect['name'] == "FramS"    : loc_sec={'name':"FramS"    ,'is': 301,'ie':333 ,'jloc':338 ,'depthlim':(-1000.,0.),'templim':(-2.,6.,0.5), 'salilim':(33.,36.,0.2) , 'velolim':(-5.,5.,1.)   ,'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
		if zsect['name'] == "FramObs"  : loc_sec={'name':"FramObs"  ,'is': 326,'ie':333 ,'jloc':333 ,'depthlim':(-1000.,0.),'templim':(-2.,6.,0.5), 'salilim':(33.,36.,0.2) , 'velolim':(-5.,5.,1.)   ,'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
		if zsect['name'] == "Davis"    : loc_sec={'name':"Davis"    ,'is': 162,'ie':190 ,'jloc':252 ,'depthlim':(-1000.,0.),'templim':(-2.,6.,0.5), 'salilim':(32.,35.,0.25), 'velolim':(-10.,10.,2.) ,'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
		if zsect['name'] == "Beaufort" : loc_sec={'name':"Beaufort" ,'is': 395,'ie':560 ,'jloc':197 ,'depthlim':(-600.,0.) ,'templim':(-2.,2.,0.2), 'salilim':(28.,35.,0.5) , 'velolim':(-5.,5.,1.)   ,'denslim':(24.0,27.9,0.2),'tkeklim':(-7.,1.5,1.)}
		if zsect['name'] == "ArcAnna"  : loc_sec={'name':"ArcAnna"  ,'is': 334,'ie':358 ,'jloc':395 ,'depthlim':(-1000.,0.),'templim':(-2.,6.,0.5), 'salilim':(33.,35.,0.1) , 'velolim':(-10.,10.,1.) ,'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.,10.)}
		if zsect['name'] == "Kara"     : loc_sec={'name':"Kara"     ,'is': 291,'ie':362 ,'jloc':447 ,'depthlim':(-1000.,0.),'templim':(-2.,4.,0.5), 'salilim':(33.5,35.,0.1), 'velolim':(-10.,10.,1.) ,'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
		if zsect['name'] == "Barents"  : loc_sec={'name':"Barents"  ,'is': 280,'ie':470 ,'jloc':355 ,'depthlim':(-300. ,0.),'templim':(-2.,6.,0.5), 'salilim':(33.,36.,0.2) , 'velolim':(-5.,5.,1.)   ,'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
	elif CONFIG == "CREG12.L75" :
		if zsect['name'] == "Bering"   : loc_sec={'name':"Bering"   ,'is':  620,'ie':  646,'jloc':1799,'depthlim':(-500.,0.) ,'templim':(-2.,3.,0.5), 'salilim':(28.,35.,0.5) , 'velolim':(-5.,5.,1.)	, 'denslim':(24.0,27.9,0.2),'tkeklim':(-7.,1.5,1.)}
		if zsect['name'] == "SouthG"   : loc_sec={'name':"SouthG"   ,'is':  161,'ie':  964,'jloc':2   ,'depthlim':(-5000.,0.),'templim':( 4.,27.,2.), 'salilim':(33.,36.,0.2) , 'velolim':(-50.,50.,10.), 'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
		if zsect['name'] == "FramS"    : loc_sec={'name':"FramS"    ,'is':  880,'ie': 1000,'jloc':1011,'depthlim':(-1000.,0.),'templim':(-2.,6.,0.5), 'salilim':(33.,36.,0.2) , 'velolim':(-5.,5.,1.)	, 'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
		if zsect['name'] == "FramObs"  : loc_sec={'name':"FramObs"  ,'is':  976,'ie':  992,'jloc':990 ,'depthlim':(-700. ,0.),'templim':(-2.,6.,0.5), 'salilim':(33.,36.,0.2) , 'velolim':(-5.,5.,1.)	, 'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
		if zsect['name'] == "Davis"    : loc_sec={'name':"Davis"    ,'is':  483,'ie':  567,'jloc':751 ,'depthlim':(-1000.,0.),'templim':(-2.,6.,0.5), 'salilim':(32.,35.,0.25), 'velolim':(-10.,10.,2.) , 'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
		if zsect['name'] == "Beaufort" : loc_sec={'name':"Beaufort" ,'is': 1181,'ie': 1675,'jloc':589 ,'depthlim':(-600. ,0.),'templim':(-2.,2.,0.2), 'salilim':(28.,35.,0.5) , 'velolim':(-5.,5.,1.)	, 'denslim':(24.0,27.9,0.2),'tkeklim':(-7.,1.,10.)}
		if zsect['name'] == "ArcAnna"  : loc_sec={'name':"ArcAnna"  ,'is': 1000,'ie': 1070,'jloc':1180,'depthlim':(-1000.,0.),'templim':(-2.,6.,0.5), 'salilim':(33.,35.,0.1) , 'velolim':(-10.,10.,1.) , 'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.,10.)}
		if zsect['name'] == "Kara"     : loc_sec={'name':"Kara"     ,'is':  870,'ie': 1083,'jloc':1336,'depthlim':(-1000.,0.),'templim':(-2.,4.,0.5), 'salilim':(33.5,35.,0.1) , 'velolim':(-10.,10.,1.), 'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
		if zsect['name'] == "Barents"  : loc_sec={'name':"Barents"  ,'is':  838,'ie': 1408,'jloc':1060,'depthlim':(-500. ,0.),'templim':(-2.,6.,0.5), 'salilim':(33.,36.,0.2) , 'velolim':(-10.,10.,2.) , 'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}

	return loc_sec
	

################################################################################################################################
def PLOT_SEC( znum_plot, X, Y, tab_clim, zcont, znorm, zcol, zisol, dens_clim=None, zdens=None, data_type=None, zgrid=None, zstrait=None, ztitle=None, zfig=None ):
################################################################################################################################

	ax=plt.subplot(znum_plot,facecolor='darkslategray')
	zfmt='%4.2f'
	if ztitle != None : plt.title(ztitle,size=9)
	if zgrid != 'gridW': 
		C = plt.contourf(X,Y,tab_clim,zcont,norm=znorm,cmap=zcol,extend='both') 
	else:
		C = plt.contourf(X,Y,tab_clim,zcont,norm=znorm,cmap=zcol)

	if zstrait['name'] == 'Davis':	    loctxt=[-55.,-400.]
	elif zstrait['name'] == 'Beaufort': loctxt=[-110.,-570.]
	elif zstrait['name'] == 'FramS':    loctxt=[8.,-700.]
	elif zstrait['name'] == 'ArcAnna':  loctxt=[54.3,-950.]
	elif zstrait['name'] == 'Barents':  loctxt=[22.,-350.]
	elif zstrait['name'] == 'Kara':     loctxt=[81.2,-950.]
	else:	 loctxt=[50.0,-500.]

	if zgrid == 'gridT' : 
		citxt=str(zstrait['templim'][2])+r' $^\circ$C'
		CS=plt.contour(X,Y,tab_clim,colors='k',linewidths=0.5,levels=[-1.,0.])
		plt.clabel(CS,fontsize=7,fmt='%4.2f')
	elif zgrid == 'gridS' : 
		citxt=str(zstrait['salilim'][2])+' PSU'
		CS=plt.contour(X,Y,tab_clim,colors='k',linewidths=0.5,levels=[31.,33.,34.8])
		plt.clabel(CS,fontsize=7,fmt='%4.2f')
	elif zgrid == 'gridV' : 
		citxt=str(zstrait['velolim'][2])+r' cm $s^{-1}$'
		CS=plt.contour(X,Y,tab_clim,colors='k',linewidths=0.5,levels=[0.])
		plt.clabel(CS,fontsize=7,fmt='%4.2f')
	elif zgrid == 'gridW' : 
		citxt=str(zstrait['tkeklim'][2])+r' m$^{2}s$^{-1}$$'
		CS=plt.contour(X,Y,tab_clim,colors='k',linewidths=0.5,levels=[1e-5])
		plt.clabel(CS,fontsize=4,fmt='%.0e')


	if zgrid != 'gridW' : plt.text(loctxt[0],loctxt[1],'C.I: '+citxt,color='w',fontsize=6)

	plt.ylim((zstrait['depthlim'][0],zstrait['depthlim'][1]))
	if zgrid == 'gridV' or zgrid == 'gridU' : 
		if zstrait['name'] == 'Kara':
			plt.xlabel('Latitude ', fontsize=7)
		else:
			plt.xlabel('Longitude ', fontsize=7)
	#if zgrid == 'gridV' or zgrid == 'gridU' : plt.xlabel('Longitude ', fontsize=7)
	plt.yticks(fontsize=6)
	plt.xticks(fontsize=6)
	plt.grid(True,linestyle='--', color='grey',alpha=0.7)

	if data_type != None :
		if zgrid == 'gridT' : 
			zmy_cblab=r'T ($^\circ$C)'
			plt.ylabel('Temperature',size=7)
			posx=0.50  ; posy=0.70
		elif zgrid == 'gridS' : 
			zmy_cblab=r'S (PSU)'
			plt.ylabel('Salinity',size=7)
			posx=0.50  ; posy=0.40
		elif zgrid == 'gridV' : 
			zmy_cblab=r'V (cm $s^{-1}$)'
			plt.ylabel('Velocities (>0 northward)',size=7)
			posx=0.50  ; posy=0.1 
		elif zgrid == 'gridW' : 
			plt.setp(ax.get_yticklabels(),visible=False)
			plt.ylabel(r'Kz ($m^2 s^{-1}$)',size=7)
			posx=0.60  ; posy=0.04 

		# basic one 
		if zgrid == 'gridW' :
			cax=zfig.add_axes([posx,posy,0.36,0.01])
			plt.setp(ax.get_yticklabels(),visible=False)
			cb=zfig.colorbar(C,cax=cax,format='%.0e',orientation='horizontal',extend='both')
			cb.ax.tick_params(labelsize=5)
		else:
			cax=zfig.add_axes([posx,posy,0.01,0.20])
			cb=zfig.colorbar(C,cax=cax,format='%.1f',orientation='vertical')
			cb.ax.tick_params(labelsize=7)
	else:
		# basic one 
		#cb=plt.colorbar(C,ax=ax,format='%.1f',orientation='vertical')
		#cb.ax.tick_params(labelsize=7)

		plt.setp(ax.get_yticklabels(),visible=False)

	return

################################################################################################################################
def CAL_VOLHEATSALTICE( data_dir, CONFIG, CASE, s_year, zstrait, Tdata_read, Sdata_read, Vdata_read, ITdata_read, IVdata_read, ze3, ze1t, ze1v, vmask3Dtime, time_dim ) :
################################################################################################################################

        # Back to Pyhton arrays indices starting at zero
        jjloc = zstrait['jloc']-1
        iis = zstrait['is']-1
        iie = zstrait['ie']-1
        print() 
        print(' 					 >>>> Section indices  :'+ zstrait['name'])
        print('						       		jloc   :'+ str(zstrait['jloc']))
        print('						       		istart :'+ str(zstrait['is']))
        print('						       		iend   :'+ str(zstrait['ie']))
        
        if zstrait['name'] == 'Bering' :
            vmask3Dtime[:,:,jjloc,iie::]=0   # The vmask is changed a little bit to remove a small open ocean area eastward to Bering
        
        # Compute mean temp/sal field at the V-point 
        ############################################
        # The calculation is done using data on the boundary (external part) and the first row (internal part)
        Tdata_read_V = npy.empty((Tdata_read.shape[0],75,Tdata_read.shape[3]))
        Tdata_read_V[:,:,:] = 0.5 * ( Tdata_read[:,:,jjloc+1,:] + Tdata_read[:,:,jjloc,:] ) * vmask3Dtime[:,:,jjloc,:]
        
        Sdata_read_V = npy.empty((Sdata_read.shape[0],75,Sdata_read.shape[3]))
        Sdata_read_V[:,:,:] = 0.5 * ( Sdata_read[:,:,jjloc+1,:] + Sdata_read[:,:,jjloc,:] ) * vmask3Dtime[:,:,jjloc,:]
        
        ITdata_read_V = npy.empty((ITdata_read.shape[0],ITdata_read.shape[2]))
        ITdata_read_V[:,:] = 0.5 * ( ITdata_read[:,jjloc+1,:] + ITdata_read[:,jjloc,:] ) * vmask3Dtime[:,0,jjloc,:]
        
        # Compute net transport, heat and salt transport
        #################################################
        net_volu_trans2D=npy.zeros((Tdata_read.shape[0],75,Tdata_read.shape[3]))
        net_volu_trans=npy.zeros((Tdata_read.shape[0]))
        net_heat_trans2D=npy.zeros((Tdata_read.shape[0],75,Tdata_read.shape[3]))
        net_heat_trans=npy.zeros((Tdata_read.shape[0]))
        net_salt_trans2D=npy.zeros((Tdata_read.shape[0],75,Tdata_read.shape[3]))
        net_salt_trans=npy.zeros((Tdata_read.shape[0]))
        net_icet_trans2D=npy.zeros((Tdata_read.shape[0],Tdata_read.shape[3]))
        net_icet_trans=npy.zeros((Tdata_read.shape[0]))
        
        net_voluIN_trans2D=npy.zeros((Tdata_read.shape[0],75,Tdata_read.shape[3]))
        net_heatIN_trans2D=npy.zeros((Tdata_read.shape[0],75,Tdata_read.shape[3]))
        net_saltIN_trans2D=npy.zeros((Tdata_read.shape[0],75,Tdata_read.shape[3]))
        net_voluIN_trans=npy.zeros((Tdata_read.shape[0]))
        net_heatIN_trans=npy.zeros((Tdata_read.shape[0]))
        net_saltIN_trans=npy.zeros((Tdata_read.shape[0]))

        rhocp = 1029.*4160.
        Sref = 34.8   ;  Tref = -0.1
        
        if zstrait['name'] == 'Bering' :
        	nfact=-1.
        else:
        	nfact=+1.

        for ti in set(npy.arange(time_dim)):
            for ii in set(npy.arange(iie-iis+1)+iis-3):
            #for ii in set(npy.arange(Tdata_read.shape[3])):
                for ik in set(npy.arange(75)):
                    net_volu_trans2D[ti,ik,ii] =  nfact*Vdata_read[ti,ik,jjloc,ii] * ze3[ik,jjloc,ii] * ze1v[jjloc,ii] * vmask3Dtime[ti,ik,jjloc,ii]  # [ m3 s-1 ] 
                    net_heat_trans2D[ti,ik,ii] =  net_volu_trans2D[ti,ik,ii] * Tdata_read_V[ti,ik,ii] * rhocp		        # [ W ]
                    #net_heat_trans2D[ti,ik,ii] =  net_volu_trans2D[ti,ik,ii] * (Tdata_read_V[ti,ik,ii]-Tref) * rhocp		# [ W ]
                    net_salt_trans2D[ti,ik,ii] =  net_volu_trans2D[ti,ik,ii] * (Sref-Sdata_read_V[ti,ik,ii])/Sref	        # [ m3 s-1 ] 
                    if net_volu_trans2D[ti,ik,ii] > 0e0 : net_voluIN_trans2D[ti,ik,ii] = net_volu_trans2D[ti,ik,ii].copy()      # [ m3 s-1 ]
                    net_heatIN_trans2D[ti,ik,ii] =  net_voluIN_trans2D[ti,ik,ii] * Tdata_read_V[ti,ik,ii] * rhocp	        # [ W ]
                    #net_heatIN_trans2D[ti,ik,ii] =  net_voluIN_trans2D[ti,ik,ii] * (Tdata_read_V[ti,ik,ii]-Tref) * rhocp	# [ W ]
                    net_saltIN_trans2D[ti,ik,ii] =  net_voluIN_trans2D[ti,ik,ii] * (Sref-Sdata_read_V[ti,ik,ii])/Sref	        # [ m3 s-1 ] 
        
            net_icet_trans2D[ti,ii] = nfact*IVdata_read[ti,jjloc,ii] * ze1v[jjloc,ii] * ITdata_read_V[ti,ii] * vmask3Dtime[ti,0,jjloc,ii] 
        
            net_volu_trans[ti] =  npy.sum(net_volu_trans2D[ti,:,:].squeeze())
            net_heat_trans[ti] =  npy.sum(net_heat_trans2D[ti,:,:].squeeze())
            net_salt_trans[ti] =  npy.sum(net_salt_trans2D[ti,:,:].squeeze())
            net_icet_trans[ti] =  npy.sum(net_icet_trans2D[ti,:].squeeze())
            net_voluIN_trans[ti] =  npy.sum(net_voluIN_trans2D[ti,:,:].squeeze())
            net_heatIN_trans[ti] =  npy.sum(net_heatIN_trans2D[ti,:,:].squeeze())
            net_saltIN_trans[ti] =  npy.sum(net_saltIN_trans2D[ti,:,:].squeeze())
        
        # Save fields to be able to reload them later
        npy.savez(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_'+zstrait['name']+'_STRAITSTrans_y'+str(s_year), \
        	  net_volu_trans=npy.array(net_volu_trans), net_heat_trans=npy.array(net_heat_trans), \
        	  net_salt_trans=npy.array(net_salt_trans), net_icet_trans=npy.array(net_icet_trans), \
        	  net_voluIN_trans=npy.array(net_voluIN_trans), net_heatIN_trans=npy.array(net_heatIN_trans), \
        	  net_saltIN_trans=npy.array(net_saltIN_trans) )

        #------------------------------------------------------------------------------------------------------------------------
        ########################################
        # Compute specific diagnostic along a short section close to FRAM strait
        ########################################
        #------------------------------------------------------------------------------------------------------------------------
        FramObs={'name':"FramObs"}
        
        if zstrait['name'] == 'FramS' : 
                zstrait = DEF_LOC_SEC( CONFIG, FramObs )
                
                # Back to Pyhton arrays indices starting at zero
                jjloc=zstrait['jloc']-1
                iis=zstrait['is']-1
                iie=zstrait['ie']-1
                print() 
                print(' 							  >>>> The strait treated is   :'+ zstrait['name'])
                print('											jloc   :'+ str(zstrait['jloc']))
                print('											istart :'+ str(zstrait['is']))
                print('											iend   :'+ str(zstrait['ie']))
                
                # define 3D tmak for the corresponding section
                e1t3D = npy.tile(ze1t,(ze3.shape[0],1,1))
                e1v3D = npy.tile(ze1v,(ze3.shape[0],1,1))
                e1te3t=e1t3D[:,:,:]*ze3[:,:,:]
                e1ve3t=e1v3D[:,:,:]*ze3[:,:,:]
                e1te3tsec=npy.squeeze(e1te3t[:,jjloc,:].copy())
                e1ve3tsec=npy.squeeze(e1ve3t[:,jjloc,:].copy())
                e1te3tsec[43::,:] = 0.
                e1ve3tsec[43::,:] = 0.
                e1te3tsec[:,0:iis] = 0.
                e1ve3tsec[:,0:iis] = 0.
                e1te3tsec[:,iie+2::] = 0.
                e1ve3tsec[:,iie+2::] = 0.
                surfsec = npy.sum(e1te3tsec)
                surfsecV= npy.sum(e1ve3tsec)
                e1te3tsectime=npy.tile(e1te3tsec,(time_dim,1,1))
                e1ve3tsectime=npy.tile(e1te3tsec,(time_dim,1,1))

                # Compute mean temp & velocity along the selected section
                #########################################################
                mean_T_FraOb=npy.empty((time_dim))   ; mean_V_FraOb=npy.empty((time_dim))
                for ti in set(npy.arange(time_dim)):
                    mean_T_FraOb[ti] =  npy.sum(e1te3tsectime[ti,:,:]*Tdata_read[ti,:,jjloc,:])/surfsec
                    mean_V_FraOb[ti] =  npy.sum(e1ve3tsectime[ti,:,:]*Vdata_read[ti,:,jjloc,:])/surfsecV

                # Save fields to be able to reload them later
                npy.savez(data_dir+'/DATA/'+CONFIG+'-'+CASE+'_FramObs_STRAITSTrans_y'+str(s_year), mean_T_FraOb=mean_T_FraOb, mean_V_FraOb=mean_V_FraOb)
	
        return

################################################################################################################################
def PLOT_TRANS_TISE( data_dir, CONFIG, CASE, lgTS_ys, lgTS_ye, strait, zncout ) :
################################################################################################################################

	# Time series lenght 
	####################
	lgtstime_dim = (lgTS_ye-lgTS_ys+1)*12

	LongTS_volu_trans= []	  ;   LongTS_heat_trans=[]     ;   LongTS_salt_trans=[]     ;   LongTS_icet_trans=[]
	LongTS_voluIN_trans= []   ;   LongTS_heatIN_trans=[]   ;   LongTS_saltIN_trans=[]
	LongTS_FramObs_T=[]       ;   LongTS_FramObs_V=[]

	# Start to read all yearly transports previously computed
	##########################################################
	lgts_year=lgTS_ys    ;	  t_months=(npy.arange(12)*30.+15.)/365.   ;   start = 1

	while  lgts_year <= lgTS_ye  :
		print('   				>>>>	 Read year:'+ str(lgts_year))

		locpath = data_dir+'/DATA/'
		locfile = CONFIG+'-'+CASE+'_'+strait['name']+'_STRAITSTrans_y'+str(lgts_year)+'.npz'
		if chkfile(locpath+locfile) :
			open_npzfile = npy.load(locpath+locfile,mmap_mode='r')
			LongTS_volu_trans = npy.append( LongTS_volu_trans, open_npzfile['net_volu_trans'] )
			LongTS_heat_trans = npy.append( LongTS_heat_trans, open_npzfile['net_heat_trans'] )
			LongTS_salt_trans = npy.append( LongTS_salt_trans, open_npzfile['net_salt_trans'] )
			LongTS_icet_trans = npy.append( LongTS_icet_trans, open_npzfile['net_icet_trans'] )
			LongTS_voluIN_trans = npy.append( LongTS_voluIN_trans, open_npzfile['net_voluIN_trans'] )
			LongTS_heatIN_trans = npy.append( LongTS_heatIN_trans, open_npzfile['net_heatIN_trans'] )
			LongTS_saltIN_trans = npy.append( LongTS_saltIN_trans, open_npzfile['net_saltIN_trans'] )
		else:
			LongTS_volu_trans = npy.append( LongTS_volu_trans, npy.arange(12)+npy.nan )
			LongTS_heat_trans = npy.append( LongTS_heat_trans, npy.arange(12)+npy.nan )
			LongTS_salt_trans = npy.append( LongTS_salt_trans, npy.arange(12)+npy.nan )
			LongTS_icet_trans = npy.append( LongTS_icet_trans, npy.arange(12)+npy.nan )
			LongTS_voluIN_trans = npy.append( LongTS_voluIN_trans, npy.arange(12)+npy.nan )
			LongTS_heatIN_trans = npy.append( LongTS_heatIN_trans, npy.arange(12)+npy.nan )
			LongTS_saltIN_trans = npy.append( LongTS_saltIN_trans, npy.arange(12)+npy.nan )

		if strait['name'] == 'FramS' :
			# Read mean T and V to compare to obs.
			locpath = data_dir+'/DATA/'
			locfile = CONFIG+'-'+CASE+'_FramObs_STRAITSTrans_y'+str(lgts_year)+'.npz'
			if chkfile(locpath+locfile) :
				open_npzfile = npy.load(locpath+locfile,mmap_mode='r')
				LongTS_FramObs_T = npy.append( LongTS_FramObs_T, open_npzfile['mean_T_FraOb'] )
				LongTS_FramObs_V = npy.append( LongTS_FramObs_V, open_npzfile['mean_V_FraOb'] )
			else:
				LongTS_FramObs_T = npy.append( LongTS_FramObs_T, npy.arange(12)+npy.nan )
				LongTS_FramObs_V = npy.append( LongTS_FramObs_V, npy.arange(12)+npy.nan )

		# Set the time axis
		y_years = npy.tile(lgts_year,12)+t_months
		if start == 1:
			time_axis = y_years
			start=0
		else:
			time_axis = npy.append( time_axis, y_years )

		lgts_year+=1

	# Plot the time-series over SEVERAL YEARS
	########################################################

	time_grid = npy.arange( lgTS_ys, lgTS_ye+2, 1., dtype=int)
	newlocsx  = npy.array(time_grid,'f')
	newlabelsx = npy.array(time_grid,'i')
	
	lgtsclimyear=str(lgTS_ys)+str(lgTS_ye)

	if strait['name'] == 'SouthG': xwind = 310  
	if strait['name'] == 'FramS' or strait['name'] == 'Davis' or strait['name'] == 'Bering' : xwind = 410	 
	
	plt.clf()

	# Plot the net volume transport
	############################################################################################################################
	ax = plt.subplot(xwind+1)
	plt.plot(time_axis, LongTS_volu_trans*1e-6, 'r', linewidth=0.7)

	if strait['name'] == 'FramS' :
		plt.plot(time_axis, LongTS_voluIN_trans*1e-6   , 'g', linewidth=0.7)
		plt.plot(time_axis, (LongTS_volu_trans-LongTS_voluIN_trans)*1e-6   , 'k', linewidth=0.7)
		plt.text(lgTS_ye+1.,-2. ,str(npy.round(npy.nanmean(LongTS_volu_trans*1e-6),decimals=2)),color='r',size=8)
		plt.text(lgTS_ye+1.,3.	,str(npy.round(npy.nanmean(LongTS_voluIN_trans*1e-6),decimals=2)),color='g',size=8)
		plt.text(lgTS_ye+1.,-6. ,str(npy.round(npy.nanmean((LongTS_volu_trans-LongTS_voluIN_trans)*1e-6),decimals=2)),color='k',size=8)
		#plt.text(lgTS_ys,7. ,'Obs: -1.6$\pm$3.9',color='m',size=8)
		plt.ylim([-10.,10.])

	elif strait['name'] == 'Davis' :
		plt.plot(time_axis, LongTS_voluIN_trans*1e-6   , 'g', linewidth=0.7)
		plt.plot(time_axis, (LongTS_volu_trans-LongTS_voluIN_trans)*1e-6   , 'k', linewidth=0.7)
		plt.text(lgTS_ye+1.,-0.5 ,str(npy.round(npy.nanmean(LongTS_volu_trans*1e-6),decimals=2)),color='r',size=8)
		plt.text(lgTS_ye+1.,2.	 ,str(npy.round(npy.nanmean(LongTS_voluIN_trans*1e-6),decimals=2)),color='g',size=8)
		plt.text(lgTS_ye+1.,-2. ,str(npy.round(npy.nanmean((LongTS_volu_trans-LongTS_voluIN_trans)*1e-6),decimals=2)),color='k',size=8)
		#plt.text(lgTS_ys,2. ,'Obs: -3.1$\pm$0.7',color='m',size=8)
		plt.ylim([-4.,4.])

	elif strait['name'] == 'Bering' :
		plt.text(lgTS_ye+1.,0.75 ,str(npy.round(npy.nanmean(LongTS_volu_trans*1e-6),decimals=2)),color='r',size=8)
		#plt.text(lgTS_ys,0.5 ,r'Obs: 1.$\pm$0.2',color='m',size=8)
		plt.ylim([0.,2.])

	elif strait['name'] == 'SouthG' :
		plt.text(lgTS_ye+1.,-2 ,str(npy.round(npy.nanmean(LongTS_volu_trans*1e-6),decimals=2)),color='r',size=8)
		plt.ylim([-4.,2.])
	
	plt.title(CASE+' integrated quantities @ '+strait['name']+' \n '+str(lgtsclimyear),size=9)
	plt.xticks(newlocsx,newlabelsx,size=6)
	plt.ylabel('Volume flux \n'+'(Sv)',size=7)
	plt.yticks(size=6)
	plt.setp(ax.get_xticklabels(),visible=False)
	plt.grid(True,linestyle='--',color='grey',alpha=0.7)
	
	# Plot the heat transport
	############################################################################################################################
	if strait['name'] == 'SouthG' :
	   z_alpha=1e-15
	else:
	   z_alpha=1e-12
	ax=plt.subplot(xwind+2)
	curve_tot = plt.plot(time_axis, LongTS_heat_trans*z_alpha, 'r', linewidth=0.7)

	if strait['name'] == 'FramS' :
		curve_in = plt.plot(time_axis, LongTS_heatIN_trans*z_alpha   , 'g', linewidth=0.7)
		curve_ou = plt.plot(time_axis, (LongTS_heat_trans-LongTS_heatIN_trans)*z_alpha	 , 'k', linewidth=0.7)
		plt.text(lgTS_ye+1.,30	,str(npy.round(npy.nanmean(LongTS_heat_trans*z_alpha),decimals=2)),color='r',size=8)
		plt.text(lgTS_ye+1.,40	,str(npy.round(npy.nanmean(LongTS_heatIN_trans*z_alpha),decimals=2)),color='g',size=8)
		plt.text(lgTS_ye+1.,10	,str(npy.round(npy.nanmean((LongTS_heat_trans-LongTS_heatIN_trans)*z_alpha),decimals=2)),color='k',size=8)
		#plt.text(lgTS_ys,40 ,r'Obs: 62.$\pm$17',color='m',size=8)
		plt.ylim([-10,50])

	elif strait['name'] == 'Davis' :
		curve_in = plt.plot(time_axis, LongTS_heatIN_trans*z_alpha   , 'g', linewidth=0.7)
		curve_ou = plt.plot(time_axis, (LongTS_heat_trans-LongTS_heatIN_trans)*z_alpha	 , 'k', linewidth=0.7)
		plt.text(lgTS_ye+1.,0 ,str(npy.round(npy.nanmean(LongTS_heat_trans*z_alpha),decimals=2)),color='r',size=8)
		plt.text(lgTS_ye+1.,20	,str(npy.round(npy.nanmean(LongTS_heatIN_trans*z_alpha),decimals=2)),color='g',size=8)
		plt.text(lgTS_ye+1.,10	,str(npy.round(npy.nanmean((LongTS_heat_trans-LongTS_heatIN_trans)*z_alpha),decimals=2)),color='k',size=8)
		#plt.text(lgTS_ys,23 ,r'Obs: 28.$\pm$3',color='m',size=8)
		plt.ylim([-20,30])

	elif strait['name'] == 'Bering' :
		plt.text(lgTS_ye+1.,30 ,str(npy.round(npy.nanmean(LongTS_heat_trans*z_alpha),decimals=2)),color='r',size=8)
		#plt.text(lgTS_ys,40 ,r'Obs: 13.$\pm$2',color='m',size=8)
		plt.ylim([-10,50.])

	elif strait['name'] == 'SouthG' :
		plt.text(lgTS_ye+1.,0.6 ,str(npy.round(npy.nanmean(LongTS_heat_trans*z_alpha),decimals=2)),color='r',size=8)
		plt.ylim([0,1.])
	
	plt.xticks(newlocsx,newlabelsx,size=6)
	plt.ylabel('Heat flux \n'+'(TW)',size=7)
	plt.setp(ax.get_xticklabels(),rotation=45)
	plt.yticks(size=6)
	if strait['name'] == 'Bering' or strait['name'] == 'FramS' or strait['name'] == 'Davis' or strait['name'] == 'SouthG' : plt.setp(ax.get_xticklabels(),visible=False)
	plt.grid(True,linestyle='--',color='grey',alpha=0.7)
	
	# Plot the salt transport
	############################################################################################################################
	if strait['name'] == 'Bering' or strait['name'] == 'FramS' or strait['name'] == 'Davis' or strait['name'] == 'SouthG':
		ax=plt.subplot(xwind+3)
		plt.plot(time_axis, LongTS_salt_trans*1e-3, 'r', label='Net', linewidth=0.7)

		if strait['name'] == 'FramS' :
			plt.plot(time_axis, LongTS_saltIN_trans*1e-3   , 'g', label='Northward', linewidth=0.7)
			plt.plot(time_axis, (LongTS_salt_trans-LongTS_saltIN_trans)*1e-3   , 'k', label='Southward', linewidth=0.7)
			plt.ylim([-160.,50.])
			plt.text(lgTS_ye+1.,-80 ,str(npy.round(npy.nanmean(LongTS_salt_trans*1e-3),decimals=2)),color='r',size=8)
			plt.text(lgTS_ye+1.,-10 ,str(npy.round(npy.nanmean(LongTS_saltIN_trans*1e-3),decimals=2)),color='g',size=8)
			plt.text(lgTS_ye+1.,-120,str(npy.round(npy.nanmean((LongTS_salt_trans-LongTS_saltIN_trans)*1e-3),decimals=2)),color='k',size=8)
			#plt.text(lgTS_ys,10 ,r'Obs: 70.$\pm$37',color='m',size=8)

		elif strait['name'] == 'Davis' :
			plt.plot(time_axis, LongTS_saltIN_trans*1e-3   , 'g', label='Northward', linewidth=0.7)
			plt.plot(time_axis, (LongTS_salt_trans-LongTS_saltIN_trans)*1e-3   , 'k', label='Southward', linewidth=0.7)
			plt.ylim([-200.,80.])
			plt.text(lgTS_ye+1.,-60 ,str(npy.round(npy.nanmean(LongTS_salt_trans*1e-3),decimals=2)),color='r',size=8)
			plt.text(lgTS_ye+1., 0 ,str(npy.round(npy.nanmean(LongTS_saltIN_trans*1e-3),decimals=2)),color='g',size=8)
			plt.text(lgTS_ye+1.,-100,str(npy.round(npy.nanmean((LongTS_salt_trans-LongTS_saltIN_trans)*1e-3),decimals=2)),color='k',size=8)
			#plt.text(lgTS_ys,25 ,r'Obs: 119.$\pm$14',color='m',size=8)

		elif strait['name'] == 'Bering' :
			plt.ylim([0.,140.])
			plt.text(lgTS_ye+1.,50 ,str(npy.round(npy.nanmean(LongTS_salt_trans*1e-3),decimals=2)),color='r',size=8)
			#plt.text(lgTS_ys,120 ,r'Obs: -72.$\pm$14',color='m',size=8)

		elif strait['name'] == 'SouthG' :
			plt.ylim([-200.,50.])
			plt.text(lgTS_ye+1.,-100 ,str(npy.round(npy.nanmean(LongTS_salt_trans*1e-3),decimals=2)),color='r',size=8)
		
		plt.ylabel('Freshwater \n'+'(mSv)',size=7)
		plt.xticks(newlocsx,newlabelsx,size=5)
		#plt.setp(ax.get_xticklabels(),visible=False)
		plt.setp(ax.get_xticklabels(),rotation=90)
		plt.yticks(size=6)
		if strait['name'] == 'FramS' or strait['name'] == 'Davis' : plt.setp(ax.get_xticklabels(),visible=False)
		plt.grid(True,linestyle='--',color='grey',alpha=0.7)

		plt.legend(loc='lower left',ncol=3)
		leg = plt.gca().get_legend()
		ltext = leg.get_texts()
		plt.setp(ltext, fontsize=5.)

	# Plot the ice transport
	############################################################################################################################
	if strait['name'] == 'FramS' or strait['name'] == 'Davis' or strait['name'] == 'Bering' :
		ax=plt.subplot(xwind+4)
		z_ialpha = 1.e-4
		plt.plot(time_axis, LongTS_icet_trans*z_ialpha, 'r', linewidth=0.7)

		if strait['name'] == 'FramS' :
		   plt.text(lgTS_ye+1.,-3.  ,str(npy.round(npy.nanmean(LongTS_icet_trans*z_ialpha),decimals=2)),color='r',size=8)
		   #plt.text(lgTS_ys,-10. ,r'Obs: -5.$\pm$2',color='m',size=8)
		   plt.ylim([-25.,1])

		elif strait['name'] == 'Davis':
		   plt.text(lgTS_ye+1.,-1.5  ,str(npy.round(npy.nanmean(LongTS_icet_trans*z_ialpha),decimals=2)),color='r',size=8)
		   plt.ylim([-7.,1])

		elif strait['name'] == 'Bering':
		   plt.text(lgTS_ye+1.,1.5  ,str(npy.round(npy.nanmean(LongTS_icet_trans*z_ialpha),decimals=2)),color='r',size=8)
		   plt.ylim([-1.,4.])

		plt.ylabel('Ice export \n '+r'(x10$^{-2}$Sv)',size=7)
		plt.xticks(newlocsx,newlabelsx,size=5)
		plt.setp(ax.get_xticklabels(),rotation=90)
		plt.yticks(size=6)
		plt.grid(True,linestyle='--',color='grey',alpha=0.7)

	plt.savefig(CONFIG+'-'+CASE+'_STRAITS_'+strait['name']+'_NetVoluHeatSalt_LGTS_y'+str(lgTS_ys)+'LASTy.png',dpi=300)

	if zncout:
		# FWC field 
		#######################
		cmd_ddate="date"  ;  get_output = subprocess.check_output(cmd_ddate)
		nc_f = './NETCDF/'+CONFIG+'-'+CASE+'_STRAITS_'+strait['name']+'_NetVoluHeatSalt_LGTS_y'+str(lgTS_ys)+'LASTy.nc'
		w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
		w_nc_fid.description = "Diagnostics have been calculated using the Arctic monitoring tool "
		w_nc_fid.date=get_output.decode("utf-8")
		w_nc_fid.createDimension('time', time_axis.shape[0])
		#w_nc_fid.createDimension('time_counter', None)

		w_nc_var = w_nc_fid.createVariable('LongTS_volu_trans', 'f4', ('time'))
		w_nc_var.long_name='Model net volume flux '
		w_nc_var.units="m3/s"
		w_nc_fid.variables['LongTS_volu_trans'][:] = LongTS_volu_trans

		w_nc_var = w_nc_fid.createVariable('LongTS_voluIN_trans', 'f4', ('time'))
		w_nc_var.long_name='Model northward volume flux '
		w_nc_var.units="m3/s"
		w_nc_fid.variables['LongTS_voluIN_trans'][:] = LongTS_voluIN_trans

		w_nc_var = w_nc_fid.createVariable('LongTS_voluOU_trans', 'f4', ('time'))
		w_nc_var.long_name='Model southward volume flux '
		w_nc_var.units="m3/s"
		w_nc_fid.variables['LongTS_voluOU_trans'][:] = LongTS_volu_trans-LongTS_voluIN_trans

		w_nc_var = w_nc_fid.createVariable('LongTS_heat_trans', 'f4', ('time'))
		w_nc_var.long_name='Model net heat flux '
		w_nc_var.units="W"
		w_nc_fid.variables['LongTS_heat_trans'][:] = LongTS_heat_trans

		w_nc_var = w_nc_fid.createVariable('LongTS_heatIN_trans', 'f4', ('time'))
		w_nc_var.long_name='Model northward heat flux '
		w_nc_var.units="W"
		w_nc_fid.variables['LongTS_heatIN_trans'][:] = LongTS_heatIN_trans

		w_nc_var = w_nc_fid.createVariable('LongTS_heatOU_trans', 'f4', ('time'))
		w_nc_var.long_name='Model southward flux '
		w_nc_var.units="W"
		w_nc_fid.variables['LongTS_heatOU_trans'][:] = LongTS_heat_trans-LongTS_heatIN_trans

		w_nc_var = w_nc_fid.createVariable('LongTS_salt_trans', 'f4', ('time'))
		w_nc_var.long_name='Model net salt flux using a reference salinity of 34.8 PSU'
		w_nc_var.units="Sv"
		w_nc_fid.variables['LongTS_salt_trans'][:] = LongTS_salt_trans

		w_nc_var = w_nc_fid.createVariable('LongTS_saltIN_trans', 'f4', ('time'))
		w_nc_var.long_name='Model northward salt flux using a reference salinity of 34.8 PSU'
		w_nc_var.units="Sv"
		w_nc_fid.variables['LongTS_saltIN_trans'][:] = LongTS_saltIN_trans

		w_nc_var = w_nc_fid.createVariable('LongTS_saltOU_trans', 'f4', ('time'))
		w_nc_var.long_name='Model southward salt flux using a reference salinity of 34.8 PSU'
		w_nc_var.units="Sv"
		w_nc_fid.variables['LongTS_saltOU_trans'][:] = LongTS_salt_trans-LongTS_saltIN_trans

		w_nc_var = w_nc_fid.createVariable('LongTS_icet_trans', 'f4', ('time'))
		w_nc_var.long_name='Model southward ice export flux '
		w_nc_var.units="Sv"
		w_nc_fid.variables['LongTS_icet_trans'][:] = LongTS_icet_trans


		w_nc_fid.close()  # close the file


	# Plot the mean Temperature through FRAM section
	################################################
	if strait['name'] == 'FramS':
		# Read observations 
		locpath='./'
		locfile='FRAM_inflow.mat'
		if chkfile(locpath+locfile) :
			FramObservations_read = sio.loadmat(locpath+locfile,squeeze_me=True)
			select_FramObs_T = npy.array(FramObservations_read['T_5degE_700m_obs'])
			select_FramObs_V = npy.array(FramObservations_read['V_5degE_700m_obs'])
		else:
			select_FramObs_T = npy.arange((2004-1997+1)*12)+npy.nan
			select_FramObs_V = npy.arange((2004-1997+1)*12)+npy.nan

		# Set the time axis for observations
		cyear=1997  ;  start=1
		while cyear <= 2004 :
			y_years=npy.tile(cyear,12)+t_months
			if start == 1:
				time_axis_obs=y_years
				start=0
			else:
				time_axis_obs=npy.append(time_axis_obs,y_years)
			cyear+=1


		xwind = 210
		
		plt.clf()
		ax=plt.subplot(xwind+1)
		plt.plot(time_axis    , LongTS_FramObs_T, 'k', linewidth=0.7)
		plt.plot(time_axis_obs, select_FramObs_T, 'g', label='Obs', linewidth=0.7)
		plt.text(lgTS_ye+1.,2.	,str(npy.round(npy.nanmean(LongTS_FramObs_T),decimals=2)),color='k',size=8)
		plt.text(2005.,3.5	,str(npy.round(npy.nanmean(select_FramObs_T),decimals=2)),color='g',size=8)
		plt.title(CASE+' Mean Temp & Velocity above 700m & 5degE @ FramObs \n '+str(lgtsclimyear),size=9)
		plt.ylabel('Mean temperature \n'+r'($^{\circ}$C)', size=7)
		plt.ylim([0.,4.])
		plt.xticks(newlocsx,newlabelsx,size=5)
		plt.setp(ax.get_xticklabels(),rotation=90)
		plt.yticks(size=6)
		plt.grid(True)
		plt.legend(loc='upper left')
		leg = plt.gca().get_legend()
		ltext = leg.get_texts()
		plt.setp(ltext, fontsize=5.)
		
		ax=plt.subplot(xwind+2)
		plt.plot(time_axis    , LongTS_FramObs_V*100 , 'k', linewidth=0.7)
		plt.plot(time_axis_obs, select_FramObs_V, 'g', label='Obs', linewidth=0.7)
		plt.text(lgTS_ye+1,10.	 ,str(npy.round(npy.nanmean(LongTS_FramObs_V)*100.,decimals=1)),color='k',size=8)
		plt.text(2005.,10.   ,str(npy.round(npy.nanmean(select_FramObs_V)     ,decimals=1)),color='g',size=8)
		plt.ylabel('Mean velocity \n '+r'(cm $s^{-1}$)', size=7)
		plt.ylim([-2.,12.])
		plt.xticks(newlocsx,newlabelsx,size=5)
		plt.setp(ax.get_xticklabels(),rotation=90)
		plt.yticks(size=6)
		plt.grid(True)

		plt.savefig(CONFIG+'-'+CASE+'_STRAITS_FramObs_meanVT_LGTS_y'+str(lgTS_ys)+'LASTy.png',dpi=300)

		if zncout:
			# FWC field 
			#######################
			cmd_ddate="date"  ;  get_output = subprocess.check_output(cmd_ddate)
			nc_f = './NETCDF/'+CONFIG+'-'+CASE+'_STRAITS_FramObs_meanVT_LGTS_y'+str(lgTS_ys)+'LASTy.nc'
			w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
			w_nc_fid.description = "Diagnostics have been calculated using the Arctic monitoring tool "
			w_nc_fid.data = "Observations time series have been retrieved from the following URL: http://www.whoi.edu/page.do?pid=30914"
			w_nc_fid.date=get_output.decode("utf-8")
			w_nc_fid.createDimension('time_mod', time_axis.shape[0])
			w_nc_fid.createDimension('time_obs', time_axis_obs.shape[0])
			#w_nc_fid.createDimension('time_counter', None)

			w_nc_var = w_nc_fid.createVariable('LongTS_FramObs_T', 'f4', ('time_mod'))
			w_nc_var.long_name='Model mean temperature above 700m & in the range [ 5degE- 8.4degE] along Fram Strait @ 79DegN'
			w_nc_var.units="DegC"
			w_nc_fid.variables['LongTS_FramObs_T'][:] = LongTS_FramObs_T

			w_nc_var = w_nc_fid.createVariable('select_FramObs_T', 'f4', ('time_obs'))
			w_nc_var.long_name='Observed monthly means temperature above 700m & in the range [ 5degE- 8.4degE] along Fram Strait @ 79DegN'
			w_nc_var.units="DegC"
			w_nc_fid.variables['select_FramObs_T'][:] = select_FramObs_T

			w_nc_var = w_nc_fid.createVariable('LongTS_FramObs_V', 'f4', ('time_mod'))
			w_nc_var.long_name='Model mean velocity above 700m & in the range [ 5degE- 8.4degE] along Fram Strait @ 79DegN'
			w_nc_var.units="m/s"
			w_nc_fid.variables['LongTS_FramObs_V'][:] = LongTS_FramObs_V

			w_nc_var = w_nc_fid.createVariable('select_FramObs_V', 'f4', ('time_obs'))
			w_nc_var.long_name='Observed monthly means velocities above 700m & in the range [ 5degE- 8.4degE] along Fram Strait @ 79DegN'
			w_nc_var.units="cm/s"
			w_nc_fid.variables['select_FramObs_V'][:] = select_FramObs_V

			w_nc_fid.close()  # close the file

	return

################################################################################################################################
def PLOT_SECTION( zCONF, zCASE, strait, Tsec, Ssec, Usec, Tsec_init, Ssec_init, Ksec, z2D, lon2D, zclimyear, zncout ) :
################################################################################################################################


	plt.clf()
	fig = plt.figure(320)
	ltitle = zCASE+' /  '+strait['name']+' Section '+' \n '+str(zclimyear)
	ititle = ' Initial state /  '+strait['name']+' Section '

	# Plot the Temperature section
	########################################################################################################################
	vmin=strait['templim'][0]    ;	 vmax=strait['templim'][1]     ;   vint=strait['templim'][2]
	zcont = npy.arange(vmin,vmax+vint,vint)
	znorm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
	zisol = npy.arange(vmin,vmax+vint,vint*2.)
	zcol = plt.get_cmap('jet')

	PLOT_SEC(321,lon2D,-z2D,Tsec ,zcont,znorm,zcol,zisol,data_type='T',zgrid='gridT',zstrait=strait,ztitle=ltitle,zfig=fig)
	PLOT_SEC(322,lon2D,-z2D,Tsec_init,zcont,znorm,zcol,zisol,zgrid='gridT',zstrait=strait,ztitle=ititle)

	# Plot the Salinity section
	########################################################################################################################
	vmin=strait['salilim'][0]    ;	 vmax=strait['salilim'][1]     ;   vint=strait['salilim'][2]
	zcont = npy.arange(vmin,vmax+vint,vint)
	znorm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
	zisol = npy.arange(vmin,vmax+vint,vint*2.)
	zcol = plt.get_cmap('Spectral_r')

	PLOT_SEC(323,lon2D,-z2D,Ssec ,zcont,znorm,zcol,zisol,data_type='S',zgrid='gridS',zstrait=strait,zfig=fig)
	PLOT_SEC(324,lon2D,-z2D,Ssec_init,zcont,znorm,zcol,zisol,zgrid='gridS',zstrait=strait,zfig=fig)

	# Plot the velocity section
	########################################################################################################################
	vmin=strait['velolim'][0]    ;	 vmax=strait['velolim'][1]     ;   vint=strait['velolim'][2]
	zcont = npy.arange(vmin,vmax+vint,vint)
	znorm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
	zisol = npy.arange(vmin,vmax+vint,vint)
	zcol = plt.get_cmap('bwr')

	PLOT_SEC(325,lon2D,-z2D,Usec*100.,zcont,znorm,zcol,zisol, data_type='V',zgrid='gridV',zstrait=strait,zfig=fig)

	# Plot the vertical diffusivity 
	########################################################################################################################
	vmin=strait['tkeklim'][0]    ;	 vmax=strait['tkeklim'][1]     ;   vint=strait['tkeklim'][2]
	zcont = npy.power(10,npy.arange(vmin,vmax,1))  ;    zcont1=zcont*5.
	zcont = npy.sort(npy.append(zcont,zcont1))
	znorm = mpl.colors.LogNorm(vmin=npy.min(zcont), vmax=npy.max(zcont))
	zisol = npy.power(10,npy.arange(vmin,vmax,1))
	zcol = plt.get_cmap('Spectral')

	PLOT_SEC(326,lon2D,-z2D,Ksec,zcont,znorm,zcol,zisol,data_type='W',zgrid='gridW',zstrait=strait,zfig=fig)

	plt.tight_layout(w_pad=4.)
	plt.savefig(zCONF+'-'+zCASE+'_SECTION_'+strait['name']+'_VTS_y'+str(zclimyear)+'.png',dpi=300)

	if zncout:
		ds_out = xr.Dataset()

		# T, S, U & Kz fields along a given section
		###########################################
		ds_out['Temp']= (('z','len'), Tsec.values.astype('float32')) 
		ds_out['Temp'].attrs['long_name']='Time mean temperature @ section '+strait['name']
		ds_out['Temp'].attrs['units']='DegC'
		
		ds_out['Sal']= (('z','len'), Ssec.values.astype('float32')) 
		ds_out['Sal'].attrs['long_name']='Time mean salinity @ section '+strait['name']
		ds_out['Sal'].attrs['units']='PSU'
		
		ds_out['Temp_init']= (('z','len'), Tsec_init.values.astype('float32')) 
		ds_out['Temp_init'].attrs['long_name']='Time mean temperature @ section '+strait['name']
		ds_out['Temp_init'].attrs['units']='DegC'
		
		ds_out['Sal_init']= (('z','len'), Ssec_init.values.astype('float32')) 
		ds_out['Sal_init'].attrs['long_name']='Time mean salinity @ section '+strait['name']
		ds_out['Sal_init'].attrs['units']='PSU'
		
		ds_out['Vel']= (('z','len'), Usec.values.astype('float32')) 
		ds_out['Vel'].attrs['long_name']='Time mean velocity across section '+strait['name']
		ds_out['Vel'].attrs['units']='m/s'
		
		ds_out['Kz']= (('z','len'), Ksec.values.astype('float32')) 
		ds_out['Kz'].attrs['long_name']='Time mean vertical diffusivity @ section '+strait['name']
		ds_out['Kz'].attrs['units']='m2/s'
		
		ds_out['Geoloc2D']= (('z','len'), lon2D.astype('float32')) 
		ds_out['Geoloc2D'].attrs['long_name']='2D Latitude or Longitude'
		ds_out['Geoloc2D'].attrs['units']='DegN or DegE'
		
		ds_out['Depth2D']= (('z','len'), -z2D.astype('float32')) 
		ds_out['Depth2D'].attrs['long_nam']='2D Depth'
		ds_out['Depth2D'].attrs['units']='m'
                
		ds_out = ds_out.set_coords(['Geoloc2D','Depth2D'])

		# Write the NetCDF file 
		ds_out.attrs['History'] = 'Diagnostics have been calculated using the Arctic monitoring tool '
		ds_out.attrs['Date'] = datetime.now().strftime("%a %b %e %H:%M:%S GMT %Y")
		nc_f = './NETCDF/'+zCONF+'-'+zCASE+'_SECTION-'+strait['name']+'_y'+zclimyear+'.nc'
		ds_out.to_netcdf(nc_f,engine='netcdf4')

	return
