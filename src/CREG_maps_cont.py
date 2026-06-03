"""
CREG_maps_cont.py

Description:
This module defines a function dedicated to set properties of plots

Author:
Claude Talandier (claude.talandier@cnrs.fr)
"""
import numpy as npy
import matplotlib.pylab as plt

################################################################################################################################
def SET_ARC_CNT( zCASE, zc_year, seas, zMyvar, zslev=0, zplot_obs=0, zdiff=0 ) :
################################################################################################################################
	"""
	Function dedicated to set contours, colors, min., max for each variable
	
	Input:
	    zCASE     : the configuration name to know the grid on which ice thicness has been interpolated 
	    zc_year   : current year
	    seas      : specify the month to plot either 'm03' or 'm09' (default = '')
	    zMyvar    : variable name 
	    zslev     : (optional) vertical index level to plot
	    zplot_obs : (optional) to speficy if the data to plot is from the model or Obs. (default = 0)
	    zdiff     : (optional) to speficy if the data to plot is an anomaly or the field (default = 0)
	
	Output:
	    contours  : contours to plot
	    limits    : minimum a miximum values of the variables
	    myticks   : ticks of the colorbar 
	    ztitle    : plot title 
	    zfile_ext : filename extension of the figure
	    my_cblab  : colorbar label
	    my_cmap   : colormap 
	    m_alpha   : multiplicative factor of a given variable
	"""
	# ----------------------------------------------------------------------

	m_alpha=1.

	# Plot information 
	######### SSH #########
	if zMyvar == 'ssh' :
	   my_cblab=r'(cm)'
	   m_alpha=100.
	   my_cmap=plt.get_cmap('coolwarm')
	   my_cmap=plt.get_cmap('Spectral_r')
	   
	   zfile_ext='_SSHClim_'
	   if zplot_obs == 1 :
	        if zc_year < 2003 :
		        ztitle=' Mean DOT from Armitage et al. 2017 \n 2003-2014'
	        else : 
		        ztitle=' Mean DOT from Armitage et al. 2017 \n '+str(zc_year)
	        vmin=-40. ; vmax=40. ; vint=2.
	   else:
	        ztitle=zCASE +' mean SSH anomaly '
	        vmin=-36. ; vmax=36. ; vint=2.

	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]			  # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	######### MIXED LAYER DEPTH #########
	if zMyvar == 'mldr10_1' :
	   my_cblab=r'(m)'
	   my_cmap=plt.get_cmap('Blues')
	   
	   zfile_ext='_MLD01Clim_'
	   if seas == 'm03' : 
	        if zplot_obs == 1 :
		        ztitle=' MIMOC climatology MLD01 \n '+seas
	        else:
		        ztitle=zCASE +' mean MLD01 over \n'+str(zc_year)+' '+seas
	        vmin=0. ; vmax=80. ; vint=10.
	   if seas == 'm09' : 
	        if zplot_obs == 1 :
		        ztitle=' MIMOC climatology MLD01 '+seas
	        else:
		        ztitle=' mean MLD01 '+seas
	        vmin=0. ; vmax=40. ; vint=5.
	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]			  # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	######### SEA-ICE VOLUME #########
	if zMyvar == 'sivolu' :
	   my_cblab=r'(m)'
	   my_cmap=plt.get_cmap('Spectral')
	   #my_cmap=plt.get_cmap('Blues')
	  
	   zfile_ext='_SITHICKClim_'
	   if zplot_obs == 1 :
	        if zc_year >= 1979 and zc_year <= 2024: 
		        ztitle=' PIOMAS mean SITHICK over \n '+str(zc_year)
	        else:
		        ztitle=' PIOMAS mean SITHICK over \n 1979-2024'
	   else:
	        ztitle=zCASE +' mean SITHICK over \n'+str(zc_year)
	   vmin=0. ; vmax=5. ; vint=0.5
	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]			  # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	######### SEA-ICE CONCENTRATION #########
	if zMyvar == 'siconc' :
	   my_cblab=r'(%)'
	   my_cmap=plt.get_cmap('Blues')
	   m_alpha=100.
	   
	   zfile_ext='_SICONClim_'
	   if zplot_obs == 1 :
	        if zc_year >= 1979 and zc_year <= 2025: 
		        ztitle=' NSIDC-v6 mean SICON '+seas+' \n '+str(zc_year)
	        else:
		        ztitle=' NSIDC-v6 mean SICON '+seas+' \n 1979-2025'
	   else:
	        ztitle=' mean SICON '+seas
	   vmin=0. ; vmax=100. ; vint=10.
	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]			  # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	######### BARTROPIC STREAMFUNCTION #########
	if zMyvar == 'sobarstf' :
	   my_cblab=r'(Sv)'
	   m_alpha=1.e-6
	   my_cmap=plt.get_cmap('coolwarm')
	   
	   zfile_ext='_PSIClim_'
	   ztitle=' mean PSI over \n'+str(zc_year)
	   vmin=-5. ; vmax=5. ; vint=.25
	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]			  # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,2.*vint)	# optional colorbar ticks (None)

	######### EKE #########
	if zMyvar == 'voeke' :
	   my_cblab=r'(log10 $m^2s^{-2}$)'
	   m_alpha=1.
	   my_cmap=plt.get_cmap('RdYlBu_r')
	   
	   zfile_ext='_EKEClim_'
	   if zslev == '0' : 
	        vmin=-6. ; vmax=-2. ; vint=0.2
	        if zplot_obs == 0 :
		        ztitle=zCASE +' - '+str(zc_year)+' mean EKE @ '+zslev+' m '
	        else:
		        if zc_year >= 2003 and zc_year <= 2014 : 
			        ztitle=' Mean EKE from Armitage et al. 2017 \n '+str(zc_year)
		        else:
			        ztitle=' Mean EKE from Armitage et al. 2017 \n 2003-2014'
	   else : 
	        vmin=-6. ; vmax=-2. ; vint=0.2
	        ztitle=zCASE +' - '+str(zc_year)+' mean EKE @ '+zslev+' m \n patches from Von Appen et al. 2022'

	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]			  # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	######### TEMPERATURE #########
	if zMyvar == 'votemper' :
	   my_cblab=r'($^\circ$C)'
	   #my_cmap=plt.get_cmap('coolwarm')
	   my_cmap=plt.get_cmap('jet')
	   
	   if zdiff == 1 :
	        zfile_ext='_TSDiffClim_z'+zslev+'m_'
	        ztitle=zCASE +' mean Temp diff with init. state over '+str(zc_year)+'\n @ depth '+zslev+' m'
	        vmin=-3. ; vmax=3. ; vint=0.5
	        my_cmap=plt.get_cmap('RdBu_r')
	        #my_cmap=plt.get_cmap('coolwarm')
	   else:
	        zfile_ext='_TClim_z'+zslev+'m_'
	        ztitle=zCASE +' mean Temp over'+str(zc_year)+'\n @ depth '+zslev+' m'
	        vmin=-2. ; vmax=8. ; vint=0.5
	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]			  # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	######### SALINITY #########
	if zMyvar == 'vosaline' :
	   my_cblab=r'(PSU)'
	   my_cmap=plt.get_cmap('jet')
	   #my_cmap=plt.get_cmap('YlOrRd')
	 
	   if zdiff == 1 :
	        zfile_ext='_TSDiffClim_z'+zslev+'m_'
	        ztitle=zCASE +' mean Sal diff with init. state over '+str(zc_year)+'\n @ depth '+zslev+' m'
	        vmin=-4. ; vmax=4. ; vint=0.5
	        if zslev == '97' or zslev == '199' : vmin=-2. ; vmax=2. ; vint=0.5
	        my_cmap=plt.get_cmap('RdBu_r')
	   else:
	        zfile_ext='_SClim_z'+zslev+'m_'
	        ztitle=zCASE +' mean Sal over'+str(zc_year)+'\n @ depth '+zslev+' m'
	        vmin=15. ; vmax=34. ; vint=1.
	        if zslev == '97' : vmin=32. ; vmax=36. ; vint=0.5
	        if zslev == '508': vmin=34. ; vmax=36. ; vint=0.2
	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]			  # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	return contours, limits, myticks, ztitle, zfile_ext, my_cblab, my_cmap, m_alpha
