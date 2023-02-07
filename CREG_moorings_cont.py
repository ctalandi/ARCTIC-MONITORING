import numpy as npy 
import matplotlib as mpl
import matplotlib.pylab as plt

def SET_CONT(lbox, lvarname, lano_cmp=0):

	if lvarname == "votemper":
		# Temperature
                if lano_cmp == 1 :
		        norm = mpl.colors.Normalize(vmin=-2., vmax=2.)
        	        contours = npy.arange(-2.,2.+1,0.2)
        	        isoline  = npy.arange(-2.,2.+1,0.2)
        	        cmap=plt.cm.get_cmap('seismic')
                else:
			vmin=lbox['templim'][0]    ;  vmax=lbox['templim'][1]    ;  vint=lbox['templim'][2]
		        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        	        contours = npy.arange(vmin,vmax+vint,vint)
        	        isoline  = npy.arange(vmin,vmax+vint,vint*2)
        	        cmap=plt.cm.get_cmap('rainbow')
        	        #cmap=plt.cm.get_cmap('jet')
        	        #cmap=plt.cm.get_cmap('RdBu_r')
        	        #cmap=plt.cm.get_cmap('inferno')
        	        #cmap=plt.cm.get_cmap('seismic')
	elif lvarname == "vosaline":
		# Salinity
                if lano_cmp == 1 :
		        norm = mpl.colors.Normalize(vmin=-0.5, vmax=0.5)
        	        contours = npy.arange(-0.5,0.5+0.1,0.05)
        	        isoline  = npy.arange(-0.5,0.5+0.1,0.05)
        	        cmap=plt.cm.get_cmap('seismic')
                else:
			vmin=lbox['sallim'][0]    ;  vmax=lbox['sallim'][1]    ;  vint=lbox['sallim'][2]
			#vmin=27.8   ;  vmax=35.    ;  vint=0.2
		        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        	        contours = npy.arange(vmin,vmax+vint,vint)
        	        isoline  = npy.arange(vmin,vmax+vint,vint*2)
        	        cmap=plt.cm.get_cmap('Spectral_r')
        	        #cmap=plt.cm.get_cmap('jet')
        	        #cmap=plt.cm.get_cmap('seismic')
        	        #cmap=plt.cm.get_cmap('YlOrRd')
        	        #cmap=plt.cm.get_cmap('viridis')
        	        #cmap=plt.cm.get_cmap('magma')
	elif lvarname == "vovecrtz":
		# Vertical velocity
                if lano_cmp == 1 :
		        norm = mpl.colors.Normalize(vmin=-0.5, vmax=0.5)
        	        contours = npy.arange(-0.5,0.5+0.1,0.05)
        	        isoline  = npy.arange(-0.5,0.5+0.1,0.05)
        	        cmap=plt.cm.get_cmap('seismic')
                else:
			vmin=lbox['Wlim'][0]    ;  vmax=lbox['Wlim'][1]    ;  vint=lbox['Wlim'][2]
		        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        	        contours = npy.arange(vmin,vmax+vint,vint)
        	        isoline  = npy.arange(vmin,vmax+vint,vint*2)
        	        cmap=plt.cm.get_cmap('cool')
        	        #cmap=plt.cm.get_cmap('jet')
        	        #cmap=plt.cm.get_cmap('seismic')
        	        #cmap=plt.cm.get_cmap('YlOrRd')
        	        #cmap=plt.cm.get_cmap('viridis')
        	        #cmap=plt.cm.get_cmap('magma')
	elif lvarname == "votkeavt":
		# Vertical diffusivity
                if lano_cmp == 1 :
		        norm = mpl.colors.Normalize(vmin=-0.5, vmax=0.5)
        	        contours = npy.arange(-0.5,0.5+0.1,0.05)
        	        isoline  = npy.arange(-0.5,0.5+0.1,0.05)
        	        cmap=plt.cm.get_cmap('seismic')
                else:
			vmin=lbox['Kzlim'][0]    ;  vmax=lbox['Kzlim'][1]    ;  vint=lbox['Kzlim'][2]
		        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        	        contours = npy.arange(vmin,vmax+vint,vint)
        	        isoline  = npy.arange(vmin,vmax+vint,vint*2)
        	        #cmap=plt.cm.get_cmap('cool')
        	        #cmap=plt.cm.get_cmap('jet')
        	        #cmap=plt.cm.get_cmap('seismic')
        	        #cmap=plt.cm.get_cmap('YlOrRd')
        	        cmap=plt.cm.get_cmap('viridis')
        	        #cmap=plt.cm.get_cmap('magma')
	elif lvarname == "rhop_sig0":
		# sigma0
                if lano_cmp == 1 :
			vmin=-0.2    ;  vmax=0.2    ;  vint=0.01
		        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        	        contours = npy.arange(vmin,vmax+vint,vint)
        	        isoline  = npy.arange(vmin,vmax+vint,vint)
        	        cmap=plt.cm.get_cmap('seismic')
                else:
			vmin=lbox['Sig0lim'][0]    ;  vmax=lbox['Sig0lim'][1]  ;  vint=lbox['Sig0lim'][2]
		        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        	        contours = npy.arange(vmin,vmax+vint,vint)
        	        isoline  = npy.arange(vmin,vmax+vint,vint)
			if lbox['box'] == 'ARC-M1' or lbox['box'] == 'NPE-B': 
        	        	contours = [23.,24.,25.0,26.,26.5, 27., 27.3, 27.5, 27.7, 28.] 
        	        	isoline  = [23.,24.,25.0,26.,26.5, 27., 27.3, 27.5, 27.7, 28.] 
			elif lbox['box'] == 'ALP-B' :
        	        	contours = [26.4,26.7,27.,27.2, 27.4, 27.6, 27.65, 27.7, 27.8] 
        	        	isoline  = [26.4,26.7,27.,27.2, 27.4, 27.6, 27.65, 27.7, 27.8] 
			elif lbox['box'] == 'GIN-B' :
        	        	contours = [27.,27.3,27.5, 27.7, 27.9, 27.95, 28., 28.05] 
        	        	isoline  = [27.,27.3,27.5, 27.7, 27.9, 27.95, 28., 28.05] 
			elif lbox['box'] == 'STG-B' :
        	        	contours = [26.,26.3,26.4,26.5,26.6,26.8,27.,27.2,27.4,27.6, 27.7] 
        	        	isoline  = [26.,26.3,26.4,26.5,26.6,26.8,27.,27.2,27.4,27.6, 27.7] 
			elif lbox['box'] == 'BRA-B' :
        	        	contours = [27.2,27.3,27.4,27.5,27.6,27.65,27.7,27.75,27.8] 
        	        	isoline  = [27.2,27.3,27.4,27.5,27.6,27.65,27.7,27.75,27.8] 
        	        cmap=plt.cm.get_cmap('jet')

	return contours, norm, cmap, isoline 
