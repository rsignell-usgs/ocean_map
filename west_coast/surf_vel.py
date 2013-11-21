"""
Created on Wed Apr 18 16:02:24 2012

@author: rsignell
"""
import netCDF4
import numpy as np
import datetime
import scipy.interpolate

def surf_vel(x,y,url,date_mid=datetime.datetime.utcnow(),uvar='u',vvar='v',isurf_layer=0,lonvar='lon',latvar='lat',
    tvar='time',hours_ave=24,lon360=False,ugrid=False,lonlat_sub=1,time_sub=1):
            
    nc=netCDF4.Dataset(url)
    lon = np.remainder((nc.variables[lonvar][:]+180),360)-180
    lat = nc.variables[latvar][:]
#rps
    print lon.ndim==1
    if ugrid:
        lon2d=lon
        lat2d=lat
    elif lon.ndim==1:
        # ai and aj are logical arrays, True in subset region
        igood = np.where((lon>=x.min()) & (lon<=x.max()))
        jgood = np.where((lat>=y.min()) & (lat<=y.max()))
        bi=np.arange(igood[0].min(),igood[0].max()+1,lonlat_sub)
        print 'i am here' 
        bj=np.arange(jgood[0].min(),jgood[0].max()+1,lonlat_sub)
        print bi
        print bj
        [lon2d,lat2d]=np.meshgrid(lon[bi],lat[bj]) 
    elif lon.ndim==2:
        igood=np.where(((lon>=x.min())&(lon<=x.max())) & ((lat>=y.min())&(lat<=y.max())))
        bj=np.arange(igood[0].min(),igood[0].max()+1,lonlat_sub)
        bi=np.arange(igood[1].min(),igood[1].max()+1,lonlat_sub)
        lon2d=nc.variables[lonvar][bj,bi]
        lat2d=nc.variables[latvar][bj,bi]
    else:
        print 'uh oh'
        
    #desired_stop_date=datetime.datetime(2011,9,9,17,00)  # specific time (UTC)
    desired_stop_date=date_mid+datetime.timedelta(0,3600.*hours_ave/2.)  
    istop = netCDF4.date2index(desired_stop_date,nc.variables[tvar],select='nearest')   
    actual_stop_date=netCDF4.num2date(nc.variables[tvar][istop],nc.variables[tvar].units)   
    start_date=actual_stop_date-datetime.timedelta(0,3600.*hours_ave)
    istart = netCDF4.date2index(start_date,nc.variables[tvar],select='nearest') 
    
    if ugrid:
        u1=np.mean(nc.variables[uvar][istart:istop:time_sub,isurf_layer,:],axis=0)
        v1=np.mean(nc.variables[vvar][istart:istop:time_sub,isurf_layer,:],axis=0)
    else:
        print('reading u...')
        u1=np.mean(nc.variables[uvar][istart:istop:time_sub,isurf_layer,bj,bi],axis=0)
        print('reading v...')
        v1=np.mean(nc.variables[vvar][istart:istop:time_sub,isurf_layer,bj,bi],axis=0)

    xx2,yy2=np.meshgrid(x,y)
#rps
    print 'lon2d: ',np.shape(lon2d)
    print 'lat2d: ',np.shape(lat2d)
    print 'u1: ',np.shape(u1)
    print 'v1: ',np.shape(v1)
    print 'xx2: ',np.shape(xx2)
    print 'yy2: ',np.shape(yy2)

    ui=scipy.interpolate.griddata((lon2d.flatten(),lat2d.flatten()),u1.flatten(),(xx2,yy2),method='linear',fill_value=0.0)
    vi=scipy.interpolate.griddata((lon2d.flatten(),lat2d.flatten()),v1.flatten(),(xx2,yy2),method='linear',fill_value=0.0)
    ui[np.isnan(ui)]=0.0
    vi[np.isnan(vi)]=0.0

    
    return ui,vi
