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
import scipy.io
import surf_vel

# define lon/lat range and resolution of interpolation grid
x0= -92.0
y0= 40.0
x1= -76.0
y1= 49.0
dx=0.05
dy=0.05

#date_mid = datetime.datetime(2011,3,1,12,0)
date_now = datetime.datetime.utcnow()

nx=((x1-x0)/dx)+1
ny=((y1-y0)/dy)+1

gridWidth= int(nx)
gridHeight= int(ny)
    
x=np.linspace(x0,x1,gridWidth)
y=np.linspace(y0,y1,gridHeight)


for nam in ['superior']:
    url=('http://michigan.glin.net:8080/thredds/dodsC/glos/glcfs/%s/ncas_his3d' % nam)
    print url
    ut,vt,u1,v1,lon2d,lat2d = surf_vel.surf_vel(x,y,url,date_mid=date_now,uvar='u',vvar='v',isurf_layer=0)
    ui=ut
    vi=vt
 
if 1: 
    for nam in ['michigan','huron','erie','ontario']:
        url=('http://michigan.glin.net:8080/thredds/dodsC/glos/glcfs/%s/ncas_his3d' % nam)
        print url
        ut,vt,u1,v1,lon2d,lat2d = surf_vel.surf_vel(x,y,url,date_mid=date_now,uvar='u',vvar='v',isurf_layer=0)
        ind = (ui==0)
        ui[ind] = ut[ind]
        vi[ind] = vt[ind]
    
ui2d=ui
vi2d=vi

ui=ui.T   # transpose to convention for javascript
vi=vi.T   # transpose
ui=ui.flatten()
vi=vi.flatten()

ui[np.isnan(ui)]=0.0
vi[np.isnan(vi)]=0.0

nvals=len(ui)

f=open('ocean-data.js', 'w')
f.write('var windData = {\n')

timestamp=(datetime.datetime.now()).strftime('%I:00 %p on %B %d, %Y')
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
