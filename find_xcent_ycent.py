from mpl_toolkits.basemap import Basemap
import mpl_toolkits
import pyproj

print "parameters needed for Vulcan Lambert conformal conic projection:\n"
print mpl_toolkits.basemap.projection_params['lcc']
map = Basemap(projection='lcc', lat_1=45.0, lat_2=33.0,
              lon_0=-97.0, lat_0=39.0,
              llcrnrlon=-122.9833065, llcrnrlat=22.26005742,
              urcrnrlon=-62.18676839, urcrnrlat=53.4190285)
print '(-122.9833065 W, 22.26005742 N) (m): ', map(-122.9833065, 22.26005742)
print '(-97.0 W, 39.0 N) (m): ', map(-97.0, 39.0)

print "\n\n==================================================\n\n"
g = pyproj.Geod(ellps='WGS84')
(lon0, lat0, az) = g.fwd(-121.76, 37.76, 270, 225000)
# (-124.31276282453157, 37.73234744779659, 88.4371269105583)
(lon1, lat1, az) = g.fwd(lon0, lat0, 180, 225000)
# (-124.31276282453157, 35.704816457104535, 0.0)

print "parameters needed for STEM 9 km Lambert conformal conic projection:\n"
map = Basemap(projection='lcc', lat_1=60.0, lat_2=30.0,
              lon_0=-98.0, lat_0=37.66,
              llcrnrlon=-124.0, llcrnrlat=32.0,
              urcrnrlon=-116.0, urcrnrlat=42.0)
print "bottom left corner of domain (lon, lat):", lon1, lat1
print 'bottom left corner of domain (x, y) (m): ', map(lon1, lat1)
