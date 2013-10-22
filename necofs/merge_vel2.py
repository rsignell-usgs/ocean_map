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


# define lon/lat range and resolution of interpolation grid
x0= -75.9
y0= 35.1
x1= -56.6
y1= 46.0
dx=0.05
dy=0.05

nx=((x1-x0)/dx)+1
ny=((y1-y0)/dy)+1

gridWidth = int(nx)
gridHeight = int(ny)
    
x=np.linspace(x0,x1,gridWidth)
y=np.linspace(y0,y1,gridHeight)

#date_mid = datetime.datetime(2011,3,1,12,0)

date_now = datetime.datetime.utcnow()

url='http://www.smast.umassd.edu:8080/thredds/dodsC/FVCOM/NECOFS/Forecasts/NECOFS_GOM3_FORECAST.nc'
print url
ut,vt,actual_date_mid_est = surf_vel.surf_vel(x,y,url,date_mid=datetime.datetime.utcnow(),lonvar='lonc',latvar='latc',isurf_layer=0,ugrid=True,time_sub=3)

ui=ut
vi=vt


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
timestamp=actual_date_mid_est.strftime('%I:00 %p on %B %d, %Y')
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

