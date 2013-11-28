# -*- coding: utf-8 -*-
# <nbformat>3</nbformat>

# <codecell>

import numpy as np
import netCDF4
import scipy.interpolate
import datetime



def shrink(a,b):
    """Return array shrunk to fit a specified shape by triming or averaging.
    
    a = shrink(array, shape)
    
    array is an numpy ndarray, and shape is a tuple (e.g., from
    array.shape). a is the input array shrunk such that its maximum
    dimensions are given by shape. If shape has more dimensions than
    array, the last dimensions of shape are fit.
    
    as, bs = shrink(a, b)
    
    If the second argument is also an array, both a and b are shrunk to
    the dimensions of each other. The input arrays must have the same
    number of dimensions, and the resulting arrays will have the same
    shape.
    Example
    -------
    
    >>> shrink(rand(10, 10), (5, 9, 18)).shape
    (9, 10)
    >>> map(shape, shrink(rand(10, 10, 10), rand(5, 9, 18)))        
    [(5, 9, 10), (5, 9, 10)]   
       
    """

    if isinstance(b, np.ndarray):
        if not len(a.shape) == len(b.shape):
            raise Exception, \
                  'input arrays must have the same number of dimensions'
        a = shrink(a,b.shape)
        b = shrink(b,a.shape)
        return (a, b)

    if isinstance(b, int):
        b = (b,)

    if len(a.shape) == 1:                # 1D array is a special case
        dim = b[-1]
        while a.shape[0] > dim:          # only shrink a
            if (dim - a.shape[0]) >= 2:  # trim off edges evenly
                a = a[1:-1]
            else:                        # or average adjacent cells
                a = 0.5*(a[1:] + a[:-1])
    else:
        for dim_idx in range(-(len(a.shape)),0):
            dim = b[dim_idx]
            a = a.swapaxes(0,dim_idx)        # put working dim first
            while a.shape[0] > dim:          # only shrink a
                if (a.shape[0] - dim) >= 2:  # trim off edges evenly
                    a = a[1:-1,:]
                if (a.shape[0] - dim) == 1:  # or average adjacent cells
                    a = 0.5*(a[1:,:] + a[:-1,:])
            a = a.swapaxes(0,dim_idx)        # swap working dim back

    return a

def rot2d(x, y, ang):
    '''rotate vectors by geometric angle'''
    xr = x*np.cos(ang) - y*np.sin(ang)
    yr = x*np.sin(ang) + y*np.cos(ang)
    return xr, yr

# <codecell>

def surf_vel_roms(x,y,url,date_mid=datetime.datetime.utcnow,hours_ave=24,tvar='ocean_time',lonlat_sub=1,time_sub=6):
    #url = 'http://testbedapps-dev.sura.org/thredds/dodsC/alldata/Shelf_Hypoxia/tamu/roms/tamu_roms.nc'

    #url='http://tds.ve.ismar.cnr.it:8080/thredds/dodsC/field2_test/run1/his'
    #####################################################################################

    nc = netCDF4.Dataset(url)
    mask = nc.variables['mask_rho'][:]
    lon_rho = nc.variables['lon_rho'][:]
    lat_rho = nc.variables['lat_rho'][:]
    anglev = nc.variables['angle'][:]

    desired_stop_date = date_mid+datetime.timedelta(0,3600.*hours_ave/2.)  # specific time (UTC)
    istop = netCDF4.date2index(desired_stop_date,nc.variables[tvar],select='nearest')   
    actual_stop_date=netCDF4.num2date(nc.variables[tvar][istop],nc.variables[tvar].units)   
    start_date=actual_stop_date-datetime.timedelta(0,3600.*hours_ave)
    istart = netCDF4.date2index(start_date,nc.variables[tvar],select='nearest') 

    uvar='u'
    vvar='v'
    isurf_layer = -1
    print('reading u...')
    u=np.mean(nc.variables[uvar][istart:istop:time_sub,isurf_layer,:,:],axis=0)
    print('reading v...')
    v=np.mean(nc.variables[vvar][istart:istop:time_sub,isurf_layer,:,:],axis=0)
    print('done reading data...')
    u = shrink(u, mask[1:-1, 1:-1].shape)
    v = shrink(v, mask[1:-1, 1:-1].shape)

    u, v = rot2d(u, v, anglev[1:-1, 1:-1])


    # <codecell>

    lon=lon_rho[1:-1,1:-1]
    lat=lat_rho[1:-1,1:-1]


    # <codecell>


    # <codecell>

    xx2,yy2=np.meshgrid(x,y)
    print('interpolating u to uniform grid...')
    ui=scipy.interpolate.griddata((lon.flatten(),lat.flatten()),u.flatten(),(xx2,yy2),method='linear',fill_value=0.0)
    print('interpolating v to uniform grid...')
    vi=scipy.interpolate.griddata((lon.flatten(),lat.flatten()),v.flatten(),(xx2,yy2),method='linear',fill_value=0.0)
    ui[np.isnan(ui)]=0.0
    vi[np.isnan(vi)]=0.0

    
    return ui,vi


