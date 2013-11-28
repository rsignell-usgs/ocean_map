# -*- coding: utf-8 -*-
# <nbformat>3</nbformat>

# <codecell>

import numpy as np
import netCDF4
import scipy.interpolate
import datetime
import roms_utils




def surf_vel_scb(x,y,url,lonlat_sub=1,time_sub=6):
    #url = 'http://testbedapps-dev.sura.org/thredds/dodsC/alldata/Shelf_Hypoxia/tamu/roms/tamu_roms.nc'

    #url='http://tds.ve.ismar.cnr.it:8080/thredds/dodsC/field2_test/run1/his'
    #####################################################################################

    nc = netCDF4.Dataset(url)
    mask = nc.variables['mask_rho'][:]
    lon_rho = nc.variables['lon_rho'][:]
    lat_rho = nc.variables['lat_rho'][:]
    anglev = nc.variables['angle'][:]


    uvar='u'
    vvar='v'
    isurf_layer = -1
    print('reading u...')
    u=np.mean(nc.variables[uvar][0:24:3,isurf_layer,:,:],axis=0)
    print('reading v...')
    v=np.mean(nc.variables[vvar][0:24:3,isurf_layer,:,:],axis=0)
    print('done reading data...')
    u = roms_utils.shrink(u, mask[1:-1, 1:-1].shape)
    v = roms_utils.shrink(v, mask[1:-1, 1:-1].shape)

    u, v = roms_utils.rot2d(u, v, anglev[1:-1, 1:-1])


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


