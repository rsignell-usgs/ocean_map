#!/usr/bin/env/python
"""
merge_vel: script to merge US IOOS Ocean Model surface currents onto a 
common grid for visualization by <http://testbedwww.sura.org/ocean>

This script attempts to extract 24 hours of data centered on the current time, and 
then averages to create a surface current "snapshot" from which the tidal 
currents have been mostly removed.  If this is not possible, the latest available
24 hour period is used.  The data is obtained from OPeNDAP servers 
at OceanNOMADS (http://www.northerngulfinstitute.org/edac/ocean_nomads.php) as
well as at various IOOS Regional Associations. 

@author: rsignell@usgs.gov
"""
import netCDF4
import numpy as np
import datetime
import scipy.interpolate
import surf_vel_roms

def surf_vel(x,y,url,uvar='u',vvar='v',isurf_layer=0,lonvar='lon',latvar='lat',
    tvar='time',hours_ave=24,lon360=False,ugrid=False,lonlat_sub=1,time_sub=1):
            
    nc=netCDF4.Dataset(url)
    lon = nc.variables[lonvar][:]-360.*lon360
    lat = nc.variables[latvar][:]
    
    if ugrid:
        lon2d=lon
        lat2d=lat
    elif lon.ndim==1:
        # ai and aj are logical arrays, True in subset region
        igood = np.where((lon>=x.min()) & (lon<=x.max()))
        jgood = np.where((lat>=y.min()) & (lat<=y.max()))
        bi=np.arange(igood[0].min(),igood[0].max(),lonlat_sub)
        bj=np.arange(jgood[0].min(),jgood[0].max(),lonlat_sub)
        [lon2d,lat2d]=np.meshgrid(lon[bi],lat[bj]) 
    elif lon.ndim==2:
        igood=np.where(((lon>=x.min())&(lon<=x.max())) & ((lat>=y.min())&(lat<=y.max())))
        bj=np.arange(igood[0].min(),igood[0].max(),lonlat_sub)
        bi=np.arange(igood[1].min(),igood[1].max(),lonlat_sub)
        lon2d=nc.variables[lonvar][bj,bi]
        lat2d=nc.variables[latvar][bj,bi]
    else:
        print 'uh oh'
        
    #desired_stop_date=datetime.datetime(2011,9,9,17,00)  # specific time (UTC)
    desired_stop_date=datetime.datetime.utcnow()+datetime.timedelta(0,3600.*hours_ave/2.)  
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
    ui=scipy.interpolate.griddata((lon2d.flatten(),lat2d.flatten()),u1.flatten(),(xx2,yy2),method='linear',fill_value=0.0)
    vi=scipy.interpolate.griddata((lon2d.flatten(),lat2d.flatten()),v1.flatten(),(xx2,yy2),method='linear',fill_value=0.0)
    ui[np.isnan(ui)]=0.0
    vi[np.isnan(vi)]=0.0

    
    return ui,vi


# define lon/lat range and resolution of interpolation grid
x0=-130.103438
y0= 20.191999
x1= -60.885558
y1= 52.807669
gridWidth=501.
gridHeight=237.

#date_mid = datetime.datetime(2011,3,1,12,0)
date_now = datetime.datetime.utcnow()

if 0:
    x0=-98.0
    y0= 5.0
    x1= -55.0
    y1= 32.0
    dx=0.10
    dy=0.10
    
    nx=((x1-x0)/dx)+1
    ny=((y1-y0)/dy)+1
    
    gridWidth= nx
    gridHeight= ny
    
x=np.linspace(x0,x1,gridWidth)
y=np.linspace(y0,y1,gridHeight)

#  Rutgers ROMS ESPRESSO 
url='http://tds.marine.rutgers.edu:8080/thredds/dodsC/roms/espresso/2009_da/his'
print url
ut,vt=surf_vel_roms.surf_vel_roms(x,y,url,date_mid=date_now,hours_ave=24,time_sub=1)
ui=ut
vi=vt

url ='http://edac-dap3.northerngulfinstitute.org/thredds/dodsC/ncom_amseas_agg/AmSeas_Aggregation_best.ncd'
print url
ut,vt = surf_vel(x,y,url,uvar='water_u',vvar='water_v',isurf_layer=0,lon360=True,lonlat_sub=2)
ind = (ui==0)
ui[ind] = ut[ind]
vi[ind] = vt[ind]

url ='http://edac-dap3.northerngulfinstitute.org/thredds/dodsC/ncom_useast_agg/US_East_Aggregation_best.ncd'
print url
ut,vt = surf_vel(x,y,url,uvar='water_u',vvar='water_v',isurf_layer=0,lon360=True,lonlat_sub=2)
ind = (ui==0)
ui[ind] = ut[ind]
vi[ind] = vt[ind]

url='http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'
print url
ut,vt = surf_vel(x,y,url,lonvar='lonc',latvar='latc',isurf_layer=0,ugrid=True,time_sub=3)
ind = (ui==0)
ui[ind] = ut[ind]
vi[ind] = vt[ind]



for nam in ['michigan','huron','erie','ontario','superior']:
    url=('http://michigan.glin.net:8080/thredds/dodsC/glos/glcfs/%s/ncas_his3d' % nam)
    print url
    ut,vt = surf_vel(x,y,url,uvar='u',vvar='v',isurf_layer=0)
    ind = (ui==0)
    ui[ind] = ut[ind]
    vi[ind] = vt[ind]
    

url='http://ecowatch.ncddc.noaa.gov/thredds/dodsC/ncom/ncom_reg7_agg/NCOM_Region_7_Aggregation_best.ncd'
print url
ut,vt = surf_vel(x,y,url,uvar='water_u',vvar='water_v',isurf_layer=0,lon360=False,lonlat_sub=1)
ind = (ui==0)
ui[ind] = ut[ind]
vi[ind] = vt[ind]

js_header = '''var windData = {
timestamp: "12:00 pm on April 18, 2012",
x0: -130.103438,
y0: 20.191999,
x1: -60.885558,
y1: 52.807669,
gridWidth: 501.0,
gridHeight: 237.0,
field: [
'''

ui=ui.T   # transpose to convention for javascript
vi=vi.T   # transpose
ui=ui.flatten()
vi=vi.flatten()
nvals=len(ui)

ui[np.isnan(ui)]=0.0
vi[np.isnan(vi)]=0.0


nvals=len(ui)

f=open('ocean-data.js', 'w')
f.write('var windData = {\n')
#f.write('timestamp: "%s",\n' % '12:00 pm on April 17, 2012')
timestamp=(datetime.datetime.now()).strftime('%I:00 %p on %b %d, %Y')
f.write('timestamp: "%s",\n' % timestamp )
f.write('x0: %12.6f,\n' % x0) 
f.write('y0: %12.6f,\n' % y0)
f.write('x1: %12.6f,\n' % x1)
f.write('y1: %12.6f,\n' % y1)
f.write('gridWidth: %6.1f,\n' % gridWidth)
f.write('gridHeight: %6.1f,\n' % gridHeight)
f.write('field: [\n')
Lines = ['%4.3f,%4.3f,\n' % (ui[i],vi[i]) for i in range(nvals-1)]
f.writelines(Lines)
f.write('%4.3f,%4.3f\n' % (ui[-1],vi[-1]))
f.write(']\n}\n')
f.close()

