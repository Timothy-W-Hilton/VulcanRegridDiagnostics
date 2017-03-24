import netCDF4
import sys
import os
from mpl_toolkits.basemap import Basemap
import pyproj

files = sys.argv[1:]
folder = os.getcwdu()
for f in files:
	nc = netCDF4.Dataset(os.path.join(folder, f))
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
        crnr_0_0 = (lon[t, v, 0,0], lat[t, v, 0,0])
        crnr_n_0 = (lon[t, v, nrows-1,0], lat[t, v, nrows-1,0])
        crnr_0_n = (lon[t, v, 0,ncols-1], lat[t, v, 0,ncols-1])
        crnr_n_n = (lon[t, v, nrows-1,ncols-1], lat[t, v, nrows-1,ncols-1])

print "using Basemap:\n"
# bmap = Basemap(projection='lcc', lat_1=30.0, lat_2=60.0,
#                lon_0=-121.76, lat_0=37.66001,
#                llcrnrlon=-124.0, llcrnrlat=32.0,
#                urcrnrlon=-116.0, urcrnrlat=42.0,
#                rsphere=(6378137.00,6356752.3142))
bmap = Basemap(projection='lcc', lat_1=30.0, lat_2=60.0,
               lon_0=-98, lat_0=37.66001,
               llcrnrlon=-123.40027, llcrnrlat=35.113796,
               urcrnrlon=-119.98279, urcrnrlat=40.192265,
               rsphere=(6378137.00,6356752.3142))

print "(0, 0)", crnr_0_0, bmap(*crnr_0_0)
print "(0, {})".format(nrows), crnr_0_n, bmap(*crnr_0_n)
print "({}, 0)".format(ncols), crnr_n_0, bmap(*crnr_n_0)
print "({}, {})".format(ncols, nrows), crnr_n_n, bmap(*crnr_n_n)

print "using pyproj:\n"
# from http://www.pkrc.net/wrf-lambert.html
# stemproj = pyproj.Proj(("+proj=lcc +lat_1=TRUELAT1 +lat_2=TRUELAT2"
#                         " +lat_0=MOAD_CEN_LAT +lon_0=STAND_LON"
#                         " +a=6370 +b=6370 +towgs84=0,0,0 +no_defs"))
stemproj = pyproj.Proj(("+proj=lcc +lat_1=30.0 +lat_2=60.0"
                        " +lat_0=37.66001 +lon_0=-98.0"
                        " +a=6378137.00 +b=6356752.3142 +towgs84=0,0,0 +no_defs"))
print "(0, 0)", crnr_0_0, stemproj(*crnr_0_0)
print "(0, {})".format(nrows), crnr_0_n, stemproj(*crnr_0_n)
print "({}, 0)".format(ncols), crnr_n_0, stemproj(*crnr_n_0)
print "({}, {})".format(ncols, nrows), crnr_n_n, stemproj(*crnr_n_n)
