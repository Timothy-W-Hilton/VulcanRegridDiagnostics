"""LCC projection: Lambert Conformal Conic
"""
import os
from mpl_toolkits.basemap import Basemap
import mpl_toolkits  # to get needed parameters for LCC projection
import pyproj
import pandas as pd

from map_vulcan_csv import vulcan_grid_mapper


def get_vulcan_domain_corners(fname='./vulcangrid.10.2012.csv'):
    """find and print the (longitude, latitude) coordinates of the four
    corners of the Vulcan grid.

    """
    df = pd.read_csv(fname)
    print df.loc[(df.i == df.i.min()) & (df.j == df.j.min())]
    print df.loc[(df.i == df.i.min()) & (df.j == df.j.max())]
    print df.loc[(df.i == df.i.max()) & (df.j == df.j.min())]
    print df.loc[(df.i == df.i.max()) & (df.j == df.j.max())]


fname_v_csv = os.path.join('/', 'project', 'projectdirs',
                           'm2319', 'Data', 'VULCAN',
                           'vulcangrid.10.2012.csv')
fname_ioapi_latlon = os.path.join('/', 'global', 'homes', 't',
                                  'twhilton', 'Code', 'Regrid',
                                  'VulcanRegrid', 'vulcan_latlon.nc')
mapper = vulcan_grid_mapper(fname_ioapi_latlon, fname_v_csv)
mapper.parse_vulcan_csv()

print "\n\n==================================================\n\n"
print "USING BASEMAP"
print "parameters needed for Vulcan Lambert conformal conic projection:\n"
print mpl_toolkits.basemap.projection_params['lcc']
map = Basemap(projection='lcc', lat_1=45.0, lat_2=33.0,
              lon_0=-97.0, lat_0=40.0,
              llcrnrlon=-122.9833065, llcrnrlat=22.26005742,
              urcrnrlon=-62.18676839, urcrnrlat=53.4190285)
print '(-122.9833065 W, 22.2605742 N) (m): ', map(-122.9833065, 22.26005742)
print '(-97.0 W, 39.0 N) (m): ', map(-97.0, 40.0)

print "\n\n==================================================\n\n"
print "USING PYPROJ"

vulcanproj = pyproj.Proj(("+proj=lcc +lat_1=33.0 +lat_2=45.0"
                          " +lat_0=40.0 +lon_0=-97.0"
                          " +a=6378137.00 +b=6356752.3142 "
                          "+towgs84=0,0,0 +no_defs"))
crnrs = mapper.get_csv_corners()
for k, this_corner in crnrs.items():
    print '{}: ({} E, {} N) (m): {}, {}'.format(k,
        *(this_corner + vulcanproj(*this_corner)))
print '(-97.0 W, 39.0 N) (m): ', vulcanproj(-97.0, 40.0)


g = pyproj.Geod(ellps='WGS84')
# Vulcan grid cells are 10 km.  Vulcan CSV provides cell centers, but
# I/O API wants the cell corner in xorig, yorig.  So we need to go
# 5000 m west and then 5000 m south from the SW corner cell in Vulcan
# CSV.  Output above shows that the southwest corner of Vulcan grid is
# in CSV cell (0, 507)
(lon0, lat0, az) = g.fwd(crnrs['0_n'][0], crnrs['0_n'][1], 270, 5000)
(lon_orig, lat_orig, az) = g.fwd(lon0, lat0, 180, 5000)
print '(lon orig, lat orig): {}, {}'.format(lon_orig, lat_orig)
print '(xorig, yorig): {}, {}'.format(*vulcanproj(lon_orig, lat_orig))
