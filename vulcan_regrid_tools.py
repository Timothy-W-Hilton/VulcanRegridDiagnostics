"""calculations, diagnostics to regrid Vulcan data to STEM 9 km grid

class vulcan_grid_mapper implements diagnostics to determine if the
I/O API GRIDDESC file
(https://www.cmascenter.org/ioapi/documentation/all_versions/html/GRIDDESC.html)
is correct. It parses Vulcan native grid longitude and latitude
coordinates from the Vulcan csv file
(http://vulcan.project.asu.edu/research.php).  It also parses latitude
and longitude grids from the I/O API latlon utility
(https://www.cmascenter.org/ioapi/documentation/all_versions/html/LATLON.html).
It then draws both grids to a global map to verify that they are the
same.

get_vulcan_domain_corners() parses the Vulcan CSV file and returns the
(longitude, latitude) coordinates of the four corners.

abbreviations:
LCC projection: Lambert Conformal Conic

"""
import os
from mpl_toolkits.basemap import Basemap
import pyproj
import pandas as pd
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
from IOAPIpytools import ioapi_pytools

class vulcan_grid_mapper(object):
    """parse data from vulcan csv file, I/O API latlon output, draw to map
    """

    def __init__(self, fname_ioapi_latlon_output, fname_vulcan_csv):
        """class constructor

        ARGS:
        fname_ioapi_latlon_output (str): full path to Vulcan CSV cile
        fname_vulcan_csv (str): full path to I/O API latlon output for
            grid VULCANGRID in Vulcan GRIDDESC file
        """
        self.fname_ioapi_latlon_output = fname_ioapi_latlon_output
        self.fname_vulcan_csv = fname_vulcan_csv

    def parse_ioapi_latlon(self):
        """ read longitude, latitude data from I/O API latlon output
        """
        nc = netCDF4.Dataset(self.fname_ioapi_latlon_output)
        self.i_lat = nc.variables['LAT'][0, 0, ...]
        self.i_lon = nc.variables['LON'][0, 0, ...]
        nc.close()

    def parse_vulcan_csv(self):
        """place longitude, latitude data from Vulcan CSV into 2-D arrays

        Data frame df must have columns i, j,
        ddx, ddy, with i and j containing the array indices, and ddx and
        ddy containing longitude and latitude, respectively.
        """
        df = pd.read_csv(self.fname_vulcan_csv)
        self.imax = df['i'].max()
        self.jmax = df['j'].max()
        self.v_lon = np.ndarray((self.imax, self.jmax))
        self.v_lat = np.ndarray((self.imax, self.jmax))
        for this_row in df.itertuples():
            self.v_lon[this_row.i - 1, this_row.j - 1] = this_row.ddx
            self.v_lat[this_row.i - 1, this_row.j - 1] = this_row.ddy

    def get_csv_corners(self):
        """get (lon, lat) coordinates of Vulcan grid corners
        """
        corners = {'0_0': (self.v_lon[0, 0], self.v_lat[0, 0]),
                   '0_{}'.format(self.jmax): (self.v_lon[0, self.jmax - 1],
                                              self.v_lat[0, self.jmax - 1]),
                   '{}_0'.format(self.imax): (self.v_lon[self.imax - 1, 0],
                                              self.v_lat[self.imax - 1, 0]),
                   '{}_{}'.format(self.imax, self.jmax): (
                       self.v_lon[self.imax - 1, self.jmax - 1],
                       self.v_lat[self.imax - 1, self.jmax - 1])}
        return corners

    def draw_map(self):
        """Draw Vulcan CSV and I/O API latlon grids to a map

        ARGS:
        v_lons (array-like): Vulcan CSV longitudes
        v_lats (array-like): Vulcan CSV latitudes
        i_lons (array-like): I/O API latlon longitudes
        i_lats (array-like): I/O API latlon latitudes
        """
        m_csv = Basemap(projection='robin',
                        lon_0=0,
                        resolution='l')  # "low" resolution
        m_csv.drawcoastlines()
        m_csv.drawmapboundary()
        m_csv.scatter(self.v_lon, self.v_lat,
                      latlon=True,
                      marker='o', color='black',
                      label='Vulcan CSV')
        m_csv.scatter(self.i_lon, self.i_lat,
                      latlon=True,
                      marker='*', color='blue',
                      label='I/O API latlon')
        plt.legend()
        return m_csv


def get_vulcan_domain_corners(fname='./vulcangrid.10.2012.csv'):
    """find and print the (longitude, latitude) coordinates of the four
    corners of the Vulcan grid.

    """
    df = pd.read_csv(fname)
    print df.loc[(df.i == df.i.min()) & (df.j == df.j.min())]
    print df.loc[(df.i == df.i.min()) & (df.j == df.j.max())]
    print df.loc[(df.i == df.i.max()) & (df.j == df.j.min())]
    print df.loc[(df.i == df.i.max()) & (df.j == df.j.max())]


def get_vulcan_griddesc_parameters(mapper):
    """calculate xorig, yorig for Vulcan GRIDDESC entry

    note: I tried using a spherical Earth
    (" +a=6370000.00 +b=6370000 "), and it made the native Vulcan--I/O
    API grid match worse.  -TWH
    """
    print "\n\n==================================================\n\n"
    print "USING PYPROJ"

    vulcanproj = pyproj.Proj(("+proj=lcc +lat_1=33.0 +lat_2=45.0"
                              " +lat_0=40.0 +lon_0=-97.0"
                              " +a=6378137.00 +b=6356752.3142 "
                              "+towgs84=0,0,0 +no_defs"))
    crnrs = mapper.get_csv_corners()
    for k, this_corner in crnrs.items():
        print '{}: ({} E, {} N) (m): {}, {}'.format(
            *((k.replace('_', ', '), ) +
              this_corner +
              vulcanproj(*this_corner)))
    print '(-97.0 W, 39.0 N) (m): ', vulcanproj(-97.0, 40.0)

    g = pyproj.Geod(ellps='WGS84')
    # Vulcan grid cells are 10 km.  Vulcan CSV provides cell centers, but
    # I/O API wants the cell corner in xorig, yorig.  So we need to go
    # 5000 m west and then 5000 m south from the SW corner cell in Vulcan
    # CSV.  Output above shows that the southwest corner of Vulcan grid is
    # in CSV cell (0, 507)
    LL = crnrs['0_{}'.format(mapper.jmax)]
    (lon0, lat0, az) = g.fwd(LL[0], LL[1], 270, 5000)
    (lon_orig, lat_orig, az) = g.fwd(lon0, lat0, 180, 5000)
    print '(lon orig, lat orig): {}, {}'.format(lon_orig, lat_orig)
    print '(xorig, yorig): {}, {}'.format(*vulcanproj(lon_orig, lat_orig))


if __name__ == "__main__":
    # parse Vulcan CSV grid coordinates, I/O API coordinates
    fname_v_csv = os.path.join('/', 'project', 'projectdirs',
                               'm2319', 'Data', 'VULCAN',
                               'vulcangrid.10.2012.csv')
    fname_ioapi_latlon = os.path.join('/', 'global', 'homes', 't',
                                      'twhilton', 'Code', 'Regrid',
                                      'VulcanRegrid', 'vulcan_latlon.nc')
    mapper = vulcan_grid_mapper(fname_ioapi_latlon, fname_v_csv)
    mapper.parse_vulcan_csv()
    mapper.parse_ioapi_latlon()
    # calculate xorig, yorig for Vulcan grid GRIDDESC entry
    get_vulcan_griddesc_parameters(mapper)
    ioapi_pytools.run_latlon(fname_griddesc='GRIDDESC_GARA',
                             fname_gridfile='vulcan_latlon.nc',
                             gridname='VULCANGRID')
    # draw map comparing Vulcan CSV grid to I/O API grid
    mapper.draw_map()
    plt.gcf().savefig('vulcan_csv_ioapi_latlon.png')

    ioapi_pytools.calculate_regrid_matrix(fname_griddesc='GRIDDESC_GARA',
                                          fname_matrix='vulcan_mat',
                                          fname_mattxt='vulcan_mat.txt',
                                          in_grid='VULCANGRID',
                                          out_grid='STEM_9KM_GRD',
                                          col_refinement=5,
                                          row_refinement=5)
    fname_vulcan_raw = os.path.join('/', 'project', 'projectdirs', 'm2319', 'transfer', 'reversed_vulcan_fossilCO2_stem9km_ioapi.nc')
    ioapi_pytools.run_regrid(fname_raw=fname_vulcan_raw, fname_regridded='vulcan_test.nc', fname_matrix='vulcan_mat', fname_mattxt='vulcan_mat.txt')
