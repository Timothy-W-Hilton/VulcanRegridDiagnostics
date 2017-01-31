import numpy as np
import netCDF4
import sys

files = sys.argv[1:]
folder = '/global/homes/g/gara/WRFV3/run/'
files = ['wrfout_d01_2015-05-24_00:00:00','wrfout_d02_2015-05-24_00:00:00','wrfout_d03_2015-05-24_00:00:00']
for i in files:
	nc=netCDF4.Dataset(folder+i)
	keys = nc.variables.keys()
	lat = nc.variables['XLAT']
	lon = nc.variables['XLONG']
	print lat.shape
	print lon.shape
	t = 0	
	print lat[t,0,0],lon[t,0,0]
	print lat[t,lat.shape[1]-1,0],lon[t,lon.shape[1]-1,0]
	print lat[t,0,lat.shape[2]-1],lon[t,0,lon.shape[2]-1]
	print lat[t,lat.shape[1]-1,lat.shape[2]-1],lon[t,lon.shape[1]-1,lon.shape[2]-1]
	print lat.shape
