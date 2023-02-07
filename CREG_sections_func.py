
import numpy as npy
import matplotlib.pylab as plt

def DEF_LOC_SEC(CONFIG,zsect) :

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
	        if zsect['name'] == "Bering"   : loc_sec={'name':"Bering"   ,'is':  620,'ie':  646,'jloc':1799,'depthlim':(-500.,0.) ,'templim':(-2.,3.,0.5), 'salilim':(28.,35.,0.5) , 'velolim':(-5.,5.,1.)   , 'denslim':(24.0,27.9,0.2),'tkeklim':(-7.,1.5,1.)}
	        if zsect['name'] == "SouthG"   : loc_sec={'name':"SouthG"   ,'is':  161,'ie':  964,'jloc':2   ,'depthlim':(-5000.,0.),'templim':( 4.,27.,2.), 'salilim':(33.,36.,0.2) , 'velolim':(-50.,50.,10.), 'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
	        if zsect['name'] == "FramS"    : loc_sec={'name':"FramS"    ,'is':  880,'ie': 1000,'jloc':1011,'depthlim':(-1000.,0.),'templim':(-2.,6.,0.5), 'salilim':(33.,36.,0.2) , 'velolim':(-5.,5.,1.)   , 'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
	        if zsect['name'] == "FramObs"  : loc_sec={'name':"FramObs"  ,'is':  976,'ie':  992,'jloc':990 ,'depthlim':(-700. ,0.),'templim':(-2.,6.,0.5), 'salilim':(33.,36.,0.2) , 'velolim':(-5.,5.,1.)   , 'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
	        if zsect['name'] == "Davis"    : loc_sec={'name':"Davis"    ,'is':  483,'ie':  567,'jloc':751 ,'depthlim':(-1000.,0.),'templim':(-2.,6.,0.5), 'salilim':(32.,35.,0.25), 'velolim':(-10.,10.,2.) , 'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
	        if zsect['name'] == "Beaufort" : loc_sec={'name':"Beaufort" ,'is': 1181,'ie': 1675,'jloc':589 ,'depthlim':(-600. ,0.),'templim':(-2.,2.,0.2), 'salilim':(28.,35.,0.5) , 'velolim':(-5.,5.,1.)   , 'denslim':(24.0,27.9,0.2),'tkeklim':(-7.,1.,10.)}
	        if zsect['name'] == "ArcAnna"  : loc_sec={'name':"ArcAnna"  ,'is': 1000,'ie': 1070,'jloc':1180,'depthlim':(-1000.,0.),'templim':(-2.,6.,0.5), 'salilim':(33.,35.,0.1) , 'velolim':(-10.,10.,1.) , 'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.,10.)}
	        if zsect['name'] == "Kara"     : loc_sec={'name':"Kara"     ,'is':  870,'ie': 1083,'jloc':1336,'depthlim':(-1000.,0.),'templim':(-2.,4.,0.5), 'salilim':(33.5,35.,0.1) , 'velolim':(-10.,10.,1.), 'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}
	        if zsect['name'] == "Barents"  : loc_sec={'name':"Barents"  ,'is':  838,'ie': 1408,'jloc':1060,'depthlim':(-500. ,0.),'templim':(-2.,6.,0.5), 'salilim':(33.,36.,0.2) , 'velolim':(-10.,10.,2.) , 'denslim':(27.6,28.9,0.1),'tkeklim':(-7.,1.5,1.)}

	return loc_sec
	

def PLOT_SEC(znum_plot,X,Y,tab_clim,zcont,znorm,zcol,zisol,dens_clim=None,zdens=None,data_type=None,zgrid=None,zstrait=None,ztitle=None,zfig=None):

	ax=plt.subplot(znum_plot,facecolor='darkslategray')
	zfmt='%4.2f'
	if ztitle != None : plt.title(ztitle,size=9)
	if zgrid != 'gridW': 
		C = plt.contourf(X,Y,tab_clim,zcont,norm=znorm,cmap=zcol,extend='both') 
	else:
		C = plt.contourf(X,Y,tab_clim,zcont,norm=znorm,cmap=zcol)

	if zstrait['name'] == 'Davis':      loctxt=[-55.,-400.]
	elif zstrait['name'] == 'Beaufort': loctxt=[-110.,-570.]
	elif zstrait['name'] == 'FramS':    loctxt=[8.,-700.]
	elif zstrait['name'] == 'ArcAnna':  loctxt=[54.3,-950.]
	elif zstrait['name'] == 'Barents':  loctxt=[22.,-350.]
	elif zstrait['name'] == 'Kara':     loctxt=[81.2,-950.]
	else:    loctxt=[50.0,-500.]

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

def CAL_VOLHEATSALTICE(data_dir,CONFIG,CASE,s_year,zstrait,Tdata_read,Sdata_read,Vdata_read,ITdata_read,IVdata_read,ze3,ze1t,ze1v,vmask3Dtime,time_dim) :

        # Back to Pyhton arrays indices starting at zero
        jjloc=zstrait['jloc']-1
        iis=zstrait['is']-1
        iie=zstrait['ie']-1
	print 
        print " >>>> The strait treated is   :", zstrait['name']
        print "                       jloc   :", zstrait['jloc']
        print "                       istart :", zstrait['is']
        print "                       iend   :", zstrait['ie']

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

        rhocp=1029.*4160.
        Sref=34.8   ;  Tref = -0.1
        
        if zstrait['name'] == 'Bering' :
                nfact=-1.
        else:
                nfact=+1.

        for ti in set(npy.arange(time_dim)):
            for ii in set(npy.arange(iie-iis+1)+iis-3):
            #for ii in set(npy.arange(Tdata_read.shape[3])):
                for ik in set(npy.arange(75)):
                    net_volu_trans2D[ti,ik,ii] =  nfact*Vdata_read[ti,ik,jjloc,ii] * ze3[ik,jjloc,ii] * ze1v[jjloc,ii] * vmask3Dtime[ti,ik,jjloc,ii]  # [ m3 s-1 ] 
                    net_heat_trans2D[ti,ik,ii] =  net_volu_trans2D[ti,ik,ii] * Tdata_read_V[ti,ik,ii] * rhocp           # [ W ]
                    #net_heat_trans2D[ti,ik,ii] =  net_volu_trans2D[ti,ik,ii] * (Tdata_read_V[ti,ik,ii]-Tref) * rhocp           # [ W ]
                    net_salt_trans2D[ti,ik,ii] =  net_volu_trans2D[ti,ik,ii] * (Sref-Sdata_read_V[ti,ik,ii])/Sref       # [ m3 s-1 ] 
		    if net_volu_trans2D[ti,ik,ii] > 0e0 : net_voluIN_trans2D[ti,ik,ii] = net_volu_trans2D[ti,ik,ii].copy()  # [ m3 s-1 ]
                    net_heatIN_trans2D[ti,ik,ii] =  net_voluIN_trans2D[ti,ik,ii] * Tdata_read_V[ti,ik,ii] * rhocp           # [ W ]
                    #net_heatIN_trans2D[ti,ik,ii] =  net_voluIN_trans2D[ti,ik,ii] * (Tdata_read_V[ti,ik,ii]-Tref) * rhocp           # [ W ]
                    net_saltIN_trans2D[ti,ik,ii] =  net_voluIN_trans2D[ti,ik,ii] * (Sref-Sdata_read_V[ti,ik,ii])/Sref       # [ m3 s-1 ] 

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
		zstrait = DEF_LOC_SEC(CONFIG,FramObs)

        	# Back to Pyhton arrays indices starting at zero
        	jjloc=zstrait['jloc']-1
        	iis=zstrait['is']-1
        	iie=zstrait['ie']-1
        	print " >>>> The strait treated is   :", zstrait['name']
        	print "                       jloc   :", zstrait['jloc']
        	print "                       istart :", zstrait['is']
        	print "                       iend   :", zstrait['ie']
		
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
