
import numpy as npy
import matplotlib.pylab as plt

# Set the contours of fields over the Arctic area

def SET_ARC_CNT(zCASE,zclimyear,seas,zMyvar,zslev=0,zplot_obs=0,zdiff=0):

	m_alpha=1.

	# Plot information 
	######### SSH #########
	if zMyvar == 'ssh' :
	   my_cblab=r'(cm)'
	   m_alpha=100.
	   my_cmap=plt.cm.get_cmap('coolwarm')
	   my_cmap=plt.cm.get_cmap('Spectral_r')
	   
	   zfile_ext='_SSHClim_'
	   if zplot_obs == 1 :
	   	if npy.int(zclimyear[0:4]) >= 2003 and npy.int(zclimyear[0:4]) <= 2014 : 
	   		ztitle=' Mean DOT from Armitage et al. 2017 \n '+str(zclimyear[0:4])
	   	else:
	   		ztitle=' Mean DOT from Armitage et al. 2017 \n 2003-2014'
	   	vmin=-40. ; vmax=40. ; vint=2.
	   else:
	   	ztitle=zCASE +' mean SSH anomaly '
	   	vmin=-36. ; vmax=36. ; vint=2.

	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]                        # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	######### MIXED LAYER DEPTH #########
	if zMyvar == 'mldr10_1' :
	   my_cblab=r'(m)'
	   my_cmap=plt.cm.get_cmap('Blues')
	   
	   zfile_ext='_MLD01Clim_'
	   if seas == 'm03' : 
	   	if zplot_obs == 1 :
	   		ztitle=' MIMOC climatology MLD01 \n '+seas
	   	else:
	   		ztitle=zCASE +' mean MLD01 over \n'+zclimyear+' '+seas
	   	vmin=0. ; vmax=80. ; vint=10.
	   if seas == 'm09' : 
	   	if zplot_obs == 1 :
	   		ztitle=' MIMOC climatology MLD01 '+seas
	   	else:
	   		ztitle=' mean MLD01 '+seas
	   	vmin=0. ; vmax=40. ; vint=5.
	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]                        # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	######### SEA-ICE VOLUME #########
	if zMyvar == 'sivolu' :
	   my_cblab=r'(m)'
	   my_cmap=plt.cm.get_cmap('Spectral')
	   #my_cmap=plt.cm.get_cmap('Blues')
	  
	   zfile_ext='_SITHICKClim_'
	   if zplot_obs == 1 :
	   	if zclimyear >= 1979 and zclimyear <= 2018: 
	   		ztitle=' PIOMAS mean SITHICK over \n '+str(zclimyear)
	   	else:
	   		ztitle=' PIOMAS mean SITHICK over \n 1979-2018'
	   else:
	   	ztitle=zCASE +' mean SITHICK over \n'+zclimyear
	   vmin=0. ; vmax=5. ; vint=0.5
	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]                        # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	######### SEA-ICE CONCENTRATION #########
	if zMyvar == 'siconc' :
	   my_cblab=r'(%)'
	   my_cmap=plt.cm.get_cmap('Blues')
	   m_alpha=100.
	   
	   zfile_ext='_SICONClim_'
	   if zplot_obs == 1 :
	   	if zclimyear >= 1979 and zclimyear <= 2015: 
	   		ztitle=' NSIDC mean SICON '+seas+' \n '+str(zclimyear)
	   	else:
	   		ztitle=' NSIDC mean SICON '+seas+' \n 1979-2015'
	   else:
	   	ztitle=' mean SICON '+seas
	   vmin=0. ; vmax=100. ; vint=10.
	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]                        # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	######### BARTROPIC STREAMFUNCTION #########
	if zMyvar == 'sobarstf' :
	   my_cblab=r'(Sv)'
	   m_alpha=1.e-6
	   my_cmap=plt.cm.get_cmap('coolwarm')
	   
	   zfile_ext='_PSIClim_'
	   ztitle=' mean PSI over \n'+zclimyear
	   vmin=-5. ; vmax=5. ; vint=1.
	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]                        # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	######### EKE #########
	if zMyvar == 'voeke' :
	   my_cblab=r'($cm^2s^{-2}$)'
	   m_alpha=1.e+4
	   my_cmap=plt.cm.get_cmap('YlOrRd')
	   
	   zfile_ext='_EKEClim_'
           if zslev == '0' : 
	   	vmin=0. ; vmax=10. ; vint=1.
	   	if zplot_obs == 0 :
	   		ztitle=zCASE +' mean EKE @ '+zslev+' m over \n'+zclimyear
	   	else:
	   		if npy.int(zclimyear[0:4]) >= 2003 and npy.int(zclimyear[0:4]) <= 2014 : 
	   			ztitle=' Mean EKE from Armitage et al. 2017 \n '+zclimyear[0:4]
	   		else:
	   			ztitle=' Mean EKE from Armitage et al. 2017 \n 2003-2014'
           if zslev == '97' : 
	   	vmin=0. ; vmax=10. ; vint=1.
	   	ztitle=' mean EKE @ '+zslev+' m '

	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]                        # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	######### TEMPERATURE #########
	if zMyvar == 'votemper' :
	   my_cblab=r'($^\circ$C)'
	   #my_cmap=plt.cm.get_cmap('coolwarm')
	   my_cmap=plt.cm.get_cmap('jet')
	   
	   if zdiff == 1 :
	   	zfile_ext='_TSDiffClim_z'+zslev+'m_'
	   	ztitle=zCASE +' mean Temp diff with init. state over '+zclimyear+'\n @ depth '+zslev+' m'
	   	vmin=-3. ; vmax=3. ; vint=0.5
	   	my_cmap=plt.cm.get_cmap('coolwarm')
	   else:
	   	zfile_ext='_TClim_z'+zslev+'m_'
	   	ztitle=zCASE +' mean Temp over'+zclimyear+'\n @ depth '+zslev+' m'
	   	vmin=-2. ; vmax=8. ; vint=0.5
	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]                        # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	######### SALINITY #########
	if zMyvar == 'vosaline' :
	   my_cblab=r'(PSU)'
	   my_cmap=plt.cm.get_cmap('jet')
	   #my_cmap=plt.cm.get_cmap('YlOrRd')
	 
	   if zdiff == 1 :
	   	zfile_ext='_TSDiffClim_z'+zslev+'m_'
	   	ztitle=zCASE +' mean Sal diff with init. state over '+zclimyear+'\n @ depth '+zslev+' m'
	   	vmin=-4. ; vmax=4. ; vint=0.5
           	if zslev == '97' : vmin=-2. ; vmax=2. ; vint=0.5
	   	my_cmap=plt.cm.get_cmap('coolwarm')
	   else:
	   	zfile_ext='_SClim_z'+zslev+'m_'
	   	ztitle=zCASE +' mean Sal over'+zclimyear+'\n @ depth '+zslev+' m'
	   	vmin=15. ; vmax=34. ; vint=1.
           	if zslev == '97' : vmin=32. ; vmax=36. ; vint=0.5
           	if zslev == '508': vmin=34. ; vmax=36. ; vint=0.2
	   contours=npy.arange(vmin,vmax+vint,vint)  # optional contours
	   limits=[vmin,vmax,vint]                        # limits for eke
	   myticks=npy.arange(vmin,vmax+vint,vint)   # optional colorbar ticks (None)

	return contours, limits, myticks, ztitle, zfile_ext, my_cblab, my_cmap, m_alpha


