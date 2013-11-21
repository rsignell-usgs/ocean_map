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
import surf_vel
import surf_vel_roms
import surf_vel_scb

# define lon/lat range and resolution of interpolation grid
x0= -140.0
y0= 25.0
x1= -115.0
y1= 52.0
dx=0.07
dy=0.07

nx=((x1-x0)/dx)+1
ny=((y1-y0)/dy)+1

gridWidth = int(nx)
gridHeight = int(ny)
    
x=np.linspace(x0,x1,gridWidth)
y=np.linspace(y0,y1,gridHeight)

#date_mid = datetime.datetime(2011,3,1,12,0)

date_now = datetime.datetime.utcnow()

#  Rutgers ROMS ESPRESSO
#url='http://oceanmodeling.pmc.ucsc.edu:8080/thredds/dodsC/ccsnrt/fmrc/CCSNRT_Aggregation_best.ncd'

#ut,vt=surf_vel_roms.surf_vel_roms(x,y,url,date_mid=date_now,hours_ave=24,time_sub=1)

# CCROMS
url ='http://thredds.axiomalaska.com/thredds/dodsC/CA_FCST.nc'
ut,vt = surf_vel.surf_vel(x,y,url,lonvar='lon',latvar='lat',isurf_layer=0,ugrid=False,time_sub=3)

print url
#ut,vt=surf_vel_roms.surf_vel_roms(x,y,url,date_mid=date_now,hours_ave=24,time_sub=1)
ui=ut
vi=vt

if 0:
    #  JPL ROMS SCB 
    url='http://ourocean.jpl.nasa.gov:8080/thredds/dodsC/SCBfcst/scb_latest_fcst_roms.nc'
    print url
    ut,vt=surf_vel_scb.surf_vel_scb(x,y,url,lonlat_sub=2,time_sub=1)
    ui=ut
    vi=vt

if 0:
    url ='http://edac-dap3.northerngulfinstitute.org/thredds/dodsC/ncom_useast_agg/US_East_Aggregation_best.ncd'
    print url
    ut,vt = surf_vel.surf_vel(x,y,url,uvar='water_u',vvar='water_v',isurf_layer=0,lon360=True,lonlat_sub=2)

    ind = (ui==0)
    ui[ind] = ut[ind]
    vi[ind] = vt[ind]

if 0:
    url='http://www.smast.umassd.edu:8080/thredds/dodsCVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'
    print url
    ut,vt = surf_vel.surf_vel(x,y,url,lonvar='lonc',latvar='latc',isurf_layer=0,ugrid=True,time_sub=3)
    ind = (ui==0)
    ui[ind] = ut[ind]
    vi[ind] = vt[ind]

if 1:
    url='http://ecowatch.ncddc.noaa.gov/thredds/dodsC/ncom/ncom_reg7_agg/NCOM_Region_7_Aggregation_best.ncd'
    print url
    ut,vt = surf_vel.surf_vel(x,y,url,uvar='water_u',vvar='water_v',isurf_layer=0,lon360=False,lonlat_sub=1)
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

