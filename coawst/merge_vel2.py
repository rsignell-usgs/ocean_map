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
import numpy as np
import datetime
import surf_vel_roms

# define lon/lat range and resolution of interpolation grid
x0= -100.0
y0= 15.0
x1= -55.0
y1= 48.0
dx=0.1
dy=0.1

nx=((x1-x0)/dx)+1
ny=((y1-y0)/dy)+1

gridWidth = int(nx)
gridHeight = int(ny)
    
x=np.linspace(x0,x1,gridWidth)
y=np.linspace(y0,y1,gridHeight)

#date_mid = datetime.datetime(2011,3,1,12,0)
date_mid = datetime.datetime.utcnow()

# load the COAWST surface currents
url='http://geoport.whoi.edu/thredds/dodsC/coawst_2_2/fmrc/coawst_2_2_best.ncd'
ut,vt,actual_date_mid_est=surf_vel_roms.surf_vel_roms(x,y,url,date_mid=date_mid,hours_ave=24,tvar='time1',lonlat_sub=1,time_sub=3)
print url
ui=ut
vi=vt

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

#scale factor
fac=1.0

ui=ui.T   # transpose to convention for javascript
vi=vi.T   # transpose
ui=ui.flatten()*fac
vi=vi.flatten()*fac
nvals=len(ui)

ui[np.isnan(ui)]=0.0
vi[np.isnan(vi)]=0.0


nvals=len(ui)

f=open('ocean-data.js', 'w')
f.write('var windData = {\n')
timestamp=actual_date_mid_est.strftime('%I:00 %p on %B %d, %Y')
f.write('timestamp: "%s",\n' % (timestamp))

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

#scipy.io.savemat('all.mat',mdict={'ui':ui,'vi':vi})