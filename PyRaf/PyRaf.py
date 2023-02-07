##########################################################################
### Copyright (c) 2011 by Raphael Dussin. For licensing, distribution 
### conditions, contact information, and additional documentation see
### in docs directory

##########################################################################
###                                                                    ###
###                          This is PyRaf                             ###
###                                                                    ###
##########################################################################
### IMPORT the packages

# Numpy
try:
        import numpy as npy
except:
        print 'numpy is not available on your machine'
	print 'check python path or install this package' ; exit()

# Scipy
try:
        from scipy.io import netcdf as nc
except:
        print 'scipy is not available on your machine'
	print 'check python path or install this package' ; exit()

# NetCDF4
try:
	from netCDF4 import Dataset
except:
	print 'NetCDF4 is not available on your machine'
	print 'NetCDF writing functions will not work '
	print 'Advice : use Enthought Python Distrib'

#######################################################
### NetCDF reading functions
#######################################################

## Read the whole array, suitable for 2d arrays
def readfullNC4(file,varname):
        fid = Dataset(file, 'r')
        out = npy.array(fid.variables[varname][:]).squeeze()
        fid.close()
        return out

## Read the whole array, suitable for 2d arrays
def readfull(file,varname):
	fid = nc.netcdf_file(file, 'r')
	out = npy.array(fid.variables[varname][:]).squeeze()
	fid.close()
	return out

# read the kt step of a (3d + time) field
def readnc_3d_tNC4(ncfile,varname,kt):
        fid = Dataset(ncfile, 'r')
        out = npy.array(fid.variables[varname][kt,:,:,:]).squeeze()
        fid.close()
        return out

# read the kt step of a (3d + time) field
def readnc_3d_t(ncfile,varname,kt):
        fid = nc.netcdf_file(ncfile, 'r')
        out = npy.array(fid.variables[varname][kt,:,:,:]).squeeze()
        fid.close()
        return out

# read the kt step of a (2d + time) field
def readnc_2d_t(ncfile,varname,kt):
        fid = nc.netcdf_file(ncfile, 'r')
        out = npy.array(fid.variables[varname][kt,:,:]).squeeze()
        fid.close()
        return out

# read the kt step of a (2d + time) field
def readnc_2d_tNC4(ncfile,varname,kt):
        fid = Dataset(ncfile, 'r')
        out = npy.array(fid.variables[varname][kt,:,:]).squeeze()
        fid.close()
        return out

# read a single 3d field
def readnc_3d(ncfile,varname):
        fid = nc.netcdf_file(ncfile, 'r')
        out = npy.array(fid.variables[varname][:,:,:]).squeeze()
        fid.close()
        return out

# read a single 2d field
def readnc_2d(ncfile,varname):
        fid = nc.netcdf_file(ncfile, 'r')
        out = npy.array(fid.variables[varname][:,:]).squeeze()
        fid.close()
        return out

# read a simple 1d field (timserie for example)
def readnc_1d(ncfile,varname):
        fid = nc.netcdf_file(ncfile, 'r')
        out = npy.array(fid.variables[varname][:]).squeeze()
        fid.close()
        return out

# read a simple section field (timserie for example)
def readnc_sec1d(ncfile,varname,sec):
        fid = nc.netcdf_file(ncfile, 'r')
        out = npy.array(fid.variables[varname][:,sec]).squeeze()
        fid.close()
        return out

## read only one level of a 3d or 4d array
def readnc_3d_lev(file,varname,level):
	fid = nc.netcdf_file(file, 'r')
	var = fid.variables[varname]
	if len(var.shape) == 2 :
		print 'This array is only 2D, use PyRaf.readnc_2d instead'
		sys.exit()
	elif len(var.shape) == 3 :
		out = npy.array(fid.variables[varname][level,:,:]).squeeze()
	elif len(var.shape) == 4 :
		out = npy.array(fid.variables[varname][:,level,:,:]).squeeze()
	else :
		print 'Unsupported number of dimensions'
		sys.exit()
	fid.close()
	return out

## read only one level of a 3d or 4d array
def readnc_3d_tlev(file,varname,kt,level):
	fid = nc.netcdf_file(file, 'r')
	var = fid.variables[varname]
	if len(var.shape) == 2 :
		print 'This array is only 2D, use PyRaf.readnc_2d instead'
		sys.exit()
	elif len(var.shape) == 3 :
		print 'This array is only 3D, use PyRaf.readnc_3d instead'
		sys.exit()
	elif len(var.shape) == 4 :
		out = npy.array(fid.variables[varname][kt,level,:,:]).squeeze()
	else :
		print 'Unsupported number of dimensions'
		sys.exit()
	fid.close()
	return out

## read only one level of a 3d or 4d array
def readnc_3d_tlevNC4(file,varname,kt,level):
        fid = Dataset(file, 'r')
        var = fid.variables[varname]
        if len(var.shape) == 2 :
                print 'This array is only 2D, use PyRaf.readnc_2d instead'
                sys.exit()
        elif len(var.shape) == 3 :
                print 'This array is only 3D, use PyRaf.readnc_3d instead'
                sys.exit()
        elif len(var.shape) == 4 :
                out = npy.array(fid.variables[varname][kt,level,:,:]).squeeze()
        else :
                print 'Unsupported number of dimensions'
                sys.exit()
        fid.close()
        return out

#######################################################
### Masking, saturation,...
#######################################################

## Mask one value (e.g. zeros)
def mask_value(tab, value):
	masque = npy.equal(tab, value)
	out = npy.ma.array(tab, mask=masque)
	return out

## Mask all values greater than given value
def mask_value_up(tab, value):
	masque = npy.greater_equal(tab, value)
	out = npy.ma.array(tab, mask=masque)
	return out

## Mask all values less than given value
def mask_value_down(tab, value):
	masque = npy.less_equal(tab, value)
	out = npy.ma.array(tab, mask=masque)
	return out

## Set all values greater than vmax to vmax
## and all values less than vmin to vmin
def limit_range(tab,vmin,vmax):
        masque = tab.mask
        tab=npy.array((tab),'f')
        tab[npy.where(tab>=vmax)] = vmax
        tab[npy.where(tab<=vmin)] = vmin
        out = npy.ma.array(tab, mask=masque)
        return out

## The same for a non-masked array
def limit_range_nomask(tab,vmin,vmax):
        tab=npy.array((tab),'f')
        tab[npy.where(tab>=vmax)] = vmax
        tab[npy.where(tab<=vmin)] = vmin
        out = tab 
        return out

## replace special value
def change_spval(tab,old,new):
	tab[npy.where(tab==old)] = new
	return tab

#######################################################
### Basic computations
#######################################################

### time operations

## compute the time mean of an 1d or 2d array
## works with non-uniformly distributed values along the time axis
def temp_mean_multiscales(array, time):
        if array.ndim == 1:
                out = temp_mean_multiscales_1d(array, time)
        elif array.ndim == 2:
                out = temp_mean_multiscales_2d(array, time)
        else:
                print 'Unable to compute mean with array w/ more than 2 dims' ; pass
        return out

def temp_mean_multiscales_1d(array, time):
        # create a dt array from time
        dt = (time[1:]-time[:-1])
        # create a extended array for the last value
        dte = npy.zeros((dt.shape[0]+1))
        dte[0:-1] = dt ; dte[-1] = dt[-1]
        # renormalize dte
        dte = dte / dte.sum()
        # this is the mean
        mean = (array * dte).sum()
        return mean

def temp_mean_multiscales_2d(array, time):
        # create a dt array from time
        dt = (time[1:]-time[:-1])
        # create a extended array for the last value
        dte = npy.zeros((dt.shape[0]+1))
        dte[0:-1] = dt ; dte[-1] = dt[-1]
        # renormalize dte
        dte = dte / dte.sum()
        dte2d = npy.zeros((array.shape))
        for kk in range(array.shape[1]):
                dte2d[:,kk] = dte
        # this is the mean
        mean = (array * dte2d).sum(0)
        return mean

### spatial operations

# compute the spatial average of a 2d field
# returns single value
def spatial_mean(array,mask,e1t=1,e2t=1):
        out = (array * mask * e1t * e2t ).sum() / ( mask * e1t * e2t ).sum()
        return out

# compute the spatial sum of a 2d field
# returns single value
def spatial_sum(array,mask,e1t=1,e2t=1):
        out = (array * mask * e1t * e2t ).sum()
        return out

# compute the zonal mean of a 2d field
# return 1d array
def zonal_mean(array,mask,e1t=1,e2t=1):
        out = (array * mask * e1t * e2t ).sum(1) / ( mask * e1t * e2t ).sum(1)
        return out

#######################################################
### Grid stuff
#######################################################

def metrics_reg(resolution,lon,lat):
        R=6371229. # from OPA phycst
        dlon = resolution * npy.pi / 180
        dlat = resolution * npy.pi / 180
        lon2d, lat2d = npy.meshgrid(lon,lat)
        e1 = R * dlon * npy.cos(npy.pi * lat2d / 180.)
        e2 = R * dlat * npy.ones((lat2d.shape))
        return e1, e2

#######################################################
### PLOTS
#######################################################

## plots are in PyRaf_plots.py
## this is keeped here (at present time) for backward compatibility
def nemo_global_plot(lon,lat,tab,contours,limits,myticks=None,name=None,filename='test.png'):
	import matplotlib.pylab as plt
	import matplotlib.cm as cm
	from matplotlib import rc, mpl
	import numpy as numpy
	from mpl_toolkits.basemap import Basemap
	#
	from matplotlib import rcParams
	rcParams['text.usetex']=True
	rcParams['text.latex.unicode']=True
	rcParams['font.family']='serif'
#	rcParams['font.serif']=['Arial']
	#rc('font',**{'family':'serif','serif':['Arial']})
	#rc('text', usetex=True)
	#rc('text', unicode=True)

	plt.rcParams['contour.negative_linestyle'] = 'solid'
	#
	tab=limit_range(tab,limits[0],limits[1])
	fig = plt.figure(figsize=[12.8 , 6.])
	ax  = fig.add_subplot(111)
	# background
	m = Basemap(projection='cyl',llcrnrlat=-80,urcrnrlat=80,\
	            llcrnrlon=0,urcrnrlon=420,resolution='c')
	m.drawcoastlines()
	m.fillcontinents(color='grey',lake_color='white')
	#m.bluemarble(scale=0.3)
	# contour filled
	norm = mpl.colors.Normalize(vmin=limits[0], vmax=limits[1])
	if len(limits) > 2:
		norm = mpl.colors.Normalize(vmin=limits[2], vmax=limits[3])
	import PyRaf_colormaps
	#pal = PyRaf_colormaps.gen_pal_Testu()
	pal = PyRaf_colormaps.gen_pal_diff()
	C = m.contourf(lon,lat,tab,25,cmap=pal,norm=norm)
#	C = m.contourf(lon,lat,tab,25,cmap=cm.gist_ncar,norm=norm)
	
	if myticks is None:
		cbar = plt.colorbar(C,format='%.2f',orientation='horizontal',shrink=0.8)
	else:
		cbar = plt.colorbar(C,format='%.2f',orientation='horizontal',shrink=0.8,ticks=myticks)
        plt.clim(limits[2],limits[3])
	# contour
	if len(contours) > 0:
		CS2 = m.contour(lon, lat, tab, contours, colors='k')
		plt.clabel(CS2, CS2.levels, inline=True, fmt='%.1f', fontsize=8)
	#
	# x axis
	locs, labels = plt.xticks()
	newlocs   = numpy.array([0,30,60,90,120,150,180,210,240,270,300,330,360,390,420],'f')
	newlabels = numpy.array([0,30,60,90,120,150,180,210,240,270,300,330,0,30,60],'i')
	plt.xticks(newlocs,newlabels)
	plt.xlabel('Longitude',fontsize=16)
	#
	# y axis
	locsY,labelsy = plt.yticks()
	newlocsy   = numpy.array([-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80],'f')
	newlabelsy = numpy.array([-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80],'i')
	plt.yticks(newlocsy,newlabelsy)
	plt.ylabel('Latitude',fontsize=16)	
	# title
	plt.title(name,fontsize=18)
	plt.savefig(filename)
	return fig

#######################################################
### Writing netCDF
#######################################################

def write_time(ncfile,time):
	# write time counter
	from netCDF4 import Dataset
	fid = Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')
	fid.description = 'file created with PyRaf (raphael.dussin@gmail.com)'
	# dimensions
	fid.createDimension('time', None)
	# variables
	times = fid.createVariable('time', 'f8', ('time',))
	# data
	times[:] = time
	fid.close()
	return None

## Write 3d field in a regular lon/lat file
## fill with lon,lat,depth and field
def write_3d_reg_file(ncfile,lon_array,lat_array,depth,time,var,varname):
	fid = Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')
	fid.description = 'file created with PyRaf (raphael.dussin@gmail.com)'
	# dimensions
	fid.createDimension('lat', lat_array.shape[0])
	fid.createDimension('lon', lon_array.shape[0])
	fid.createDimension('depth', depth.shape[0])
	fid.createDimension('time', None)
	# variables
	latitudes  = fid.createVariable('lat', 'f4', ('lat',))
	longitudes = fid.createVariable('lon', 'f4', ('lon',))
	zprof      = fid.createVariable('depth', 'f4', ('depth',))
	times      = fid.createVariable('time', 'f4', ('time',))
	variable   = fid.createVariable(varname, 'f8', ('time','lat','lon',))
	# data
	latitudes[:]    = lat_array
	longitudes[:]   = lon_array
	zprof[:]        = depth
	times[:]        = time
	variable[0,:,:,:] = var
	# close
	fid.close()
	return None

## Write 2d field in a regular lon/lat file
## fill with lon,lat and field
def write_2d_reg_file(ncfile,lon_array,lat_array,time,var,varname):
	fid = Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')
	fid.description = 'file created with PyRaf (raphael.dussin@gmail.com)'
	# dimensions
	fid.createDimension('lat', lat_array.shape[0])
	fid.createDimension('lon', lon_array.shape[0])
	fid.createDimension('time', None)
	# variables
	latitudes  = fid.createVariable('lat', 'f4', ('lat',))
	longitudes = fid.createVariable('lon', 'f4', ('lon',))
	times      = fid.createVariable('time', 'f4', ('time',))
	variable   = fid.createVariable(varname, 'f8', ('time','lat','lon',))
	# data
	latitudes[:]    = lat_array
	longitudes[:]   = lon_array
	times[:]        = time
	variable[0,:,:] = var
	# close
	fid.close()
	return None

## Write 3d field in a NEMO-like file
## fill with lon,lat,depth and field
def write_3d_nemo_file(ncfile,lon_array,lat_array,z_profile,time,var,varname):
	fid = Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')
	fid.description = 'file created with PyRaf (raphael.dussin@gmail.com)'
	# dimensions
	fid.createDimension('z', var.shape[0])
	fid.createDimension('y', var.shape[1])
	fid.createDimension('x', var.shape[2])
	fid.createDimension('time_counter', None)
	# variables
	latitudes  = fid.createVariable('nav_lat', 'f4', ('y','x',))
	longitudes = fid.createVariable('nav_lon', 'f4', ('y','x',))
	zprof      = fid.createVariable('z', 'f4', ('z',))
	times      = fid.createVariable('time_counter', 'f4', ('time_counter',))
	variable   = fid.createVariable(varname, 'f8', ('time_counter','z','y','x',))
	# data
	latitudes[:,:]    = lat_array
	longitudes[:,:]   = lon_array
        zprof[:]          = z_profile
	times[:]          = time
	variable[0,:,:,:] = var
	# close
	fid.close()
	return None

## Write 2d field in a NEMO-like file
## fill with lon,lat and field
def write_2d_nemo_file(ncfile,lon_array,lat_array,time,var,varname):
	fid = Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')
	fid.description = 'file created with PyRaf (raphael.dussin@gmail.com)'
	# dimensions
	fid.createDimension('y', var.shape[0])
	fid.createDimension('x', var.shape[1])
	fid.createDimension('time_counter', None)
	# variables
	latitudes  = fid.createVariable('nav_lat', 'f4', ('y','x',))
	longitudes = fid.createVariable('nav_lon', 'f4', ('y','x',))
	times      = fid.createVariable('time_counter', 'f4', ('time_counter',))
	variable   = fid.createVariable(varname, 'f4', ('time_counter','y','x',))
	# data
	latitudes[:,:]    = lat_array
	longitudes[:,:]   = lon_array
	times[:]          = time
	variable[0,:,:] = var
	# close
	fid.close()
	return None

## Write 2d field in a NEMO-like file
## fill with lon,lat and field
def write_2dMV_nemo_file(ncfile,lon_array,lat_array,time,var,varname,missing):
	fid = Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')
	fid.description = 'file created with PyRaf (raphael.dussin@gmail.com)'
	# dimensions
	fid.createDimension('y', var.shape[0])
	fid.createDimension('x', var.shape[1])
	fid.createDimension('time_counter', None)
	# variables
	latitudes  = fid.createVariable('nav_lat', 'f4', ('y','x',))
	longitudes = fid.createVariable('nav_lon', 'f4', ('y','x',))
	times      = fid.createVariable('time_counter', 'f4', ('time_counter',))
	variable   = fid.createVariable(varname, 'f4', ('time_counter','y','x',))
	variable.missing_value=missing
	#atts={'missing_value':missing}
        #fid.setncatts(atts)
        #fid.setncattr('missing_value',missing)
	# data
	latitudes[:,:]    = lat_array
	longitudes[:,:]   = lon_array
	times[:]          = time
	variable[0,:,:] = var
	# close
	fid.close()
	return None

## Write 1d z-profile in a NEMO-like file
def write_1d_profile(ncfile,lon_array,lat_array,depth,time,var,varname):
	fid = Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')
	fid.description = 'file created with PyRaf (raphael.dussin@gmail.com)'
	# dimensions
	fid.createDimension('y', 1)
	fid.createDimension('x', 1)
	fid.createDimension('z', var.shape[0])
	fid.createDimension('time_counter', None)
	# variables
	latitudes  = fid.createVariable('nav_lat', 'f4', ('y','x',))
	longitudes = fid.createVariable('nav_lon', 'f4', ('y','x',))
	zprof      = fid.createVariable('z', 'f4', ('z',))
	times      = fid.createVariable('time_counter', 'f4', ('time_counter',))
	variable   = fid.createVariable(varname, 'f4', ('time_counter','y','x',))
	# data
	latitudes[:,:]    = lat_array
	longitudes[:,:]   = lon_array
	times[:]          = time
        zprof[:]          = depth
	variable[:,:,:]   = var
	# close
	fid.close()
	return None

## Write 1d timeserie in a NEMO-like file
def write_1d_timeserie(ncfile,lon_array,lat_array,time,var,varname):
	fid = Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')
	fid.description = 'file created with PyRaf (raphael.dussin@gmail.com)'
	# dimensions
	fid.createDimension('y', 1)
	fid.createDimension('x', 1)
	fid.createDimension('time_counter', time.shape[0])
	# variables
	latitudes  = fid.createVariable('nav_lat', 'f4', ('y','x',))
	longitudes = fid.createVariable('nav_lon', 'f4', ('y','x',))
	times      = fid.createVariable('time_counter', 'f4', ('time_counter',))
	variable   = fid.createVariable(varname, 'f4', ('time_counter','y','x',))
	# data
	latitudes[:,:]    = lat_array
	longitudes[:,:]   = lon_array
	times[:]          = time
	variable[:,:,:]   = var
	# close
	fid.close()
	return None

#######################################################
### My Bluemarble
#######################################################

## My version of bluemarble :
## Allows the use of any background image
## optional definition of image coordinates (for splitted images)
## Use with a standard basemap instance bmap
def bluemarble_hr(bmap,image,imlonb=-180.,imlone=180.,imlatb=-90.,imlate=90.):
        " bmap is a basemap instance "
        " image is full path to image "
        " imlonb (image lon begin) is optional and default is -180. "
        " imlone (image lon end)   is optional and default is  180. "
        " imlatb (image lat begin) is optional and default is  -90. "
        " imlate (image lat end)   is optional and default is   90. "
        from PIL import Image
        from matplotlib.image import pil_to_array
        #
        # read full image
        pilImage = Image.open(image)
        width, height = pilImage.size
        nlons = width ; nlats = height
        #
        # define lon and lat arrays
        delta = (imlone - imlonb) / float(nlons)
        lons = np.arange(imlonb + 0.5*delta,imlone,delta)
        lats = np.arange(imlatb + 0.5*delta,imlate,delta)

        ilonw = npy.argmin(npy.abs(lons - bmap.llcrnrlon)) - 1
        ilone = npy.argmin(npy.abs(lons - bmap.urcrnrlon)) + 1
        ilatn = npy.argmin(npy.abs(lats - bmap.llcrnrlat)) - 1
        ilats = npy.argmin(npy.abs(lats - bmap.urcrnrlat)) + 1
        #
        # reduce image size
        image_domain = pilImage.crop((ilonw,height-ilats,ilone,height-ilatn))
        background   = pil_to_array(image_domain)
        im = bmap.imshow(background)
        return im

#######################################################
### The end...
### If you are happy with this lib, 
### have bugs to report,...
### write to raphael.dussin@gmail.com 
#######################################################
