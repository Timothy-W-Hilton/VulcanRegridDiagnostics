import numpy as np
import netCDF4
import sys
import os

files = sys.argv[1:]
folder = os.getcwdu()
for f in files:
	nc=netCDF4.Dataset(os.path.join(folder, f))
	keys = nc.variables.keys()
	lat = nc.variables['LAT'][:]
	lon = nc.variables['LON'][:]
        nc.close()
	print 'lat shape: {}'.format(lat.shape)
	print 'lon shape: {}'.format(lon.shape)
	t = 0  # time step 0
        v = 0  # vertical cell 0
        nrows = lat.shape[2]
        ncols = lat.shape[3]
	print lat[t, v, 0,0],lon[t, v, 0,0]
	print lat[t, v, nrows-1,0],lon[t, v, nrows-1,0]
	print lat[t, v, 0,ncols-1],lon[t, v, 0,ncols-1]
	print lat[t, v, nrows-1,ncols-1],lon[t, v, nrows-1,ncols-1]
