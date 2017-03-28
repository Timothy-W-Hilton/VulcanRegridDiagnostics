import os
import numpy as np
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


class vulcan_grid_mapper(object):
    """parse data from vulcan csv file, I/O API latlon output
    """

    def __init__(self, fname_ioapi_latlon_output, fname_vulcan_csv):
        self.fname_ioapi_latlon_output = fname_ioapi_latlon_output
        self.fname_vulcan_csv = fname_vulcan_csv

    def parse_ioapi_latlon(self, fname_ioapi_latlon_output):
        """ read longitude, latitude data from I/O API latlon
        """
        nc = netCDF4.Dataset(self.fname_ioapi_latlon_output)
        self.i_lat = nc.variables['LAT'][0, 0, ...]
        self.i_lon = nc.variables['LON'][0, 0, ...]
        nc.close()

    def parse_vulcan_csv(self, fname_v_csv):
        """place longitude, latitude data from Vulcan CSV into 2-D arrays

        Data frame df must have columns i, j,
        ddx, ddy, with i and j containing the array indices, and ddx and
        ddy containing longitude and latitude, respectively.
        """
        df = pd.read_csv(fname_v_csv)
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
                   '0_n': (self.v_lon[0, self.jmax - 1],
                           self.v_lat[0, self.jmax - 1]),
                   'n_0': (self.v_lon[self.imax - 1, 0],
                           self.v_lat[self.imax - 1, 0]),
                   'n_n': (self.v_lon[self.imax - 1, self.jmax - 1],
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


if __name__ == "__main__":
    fname_v_csv = os.path.join('/', 'project', 'projectdirs',
                               'm2319', 'Data', 'VULCAN',
                               'vulcangrid.10.2012.csv')
    fname_ioapi_latlon = os.path.join('/', 'global', 'homes', 't',
                                      'twhilton', 'Code', 'Regrid',
                                      'VulcanRegrid', 'vulcan_latlon.nc')
    mapper = vulcan_grid_mapper(fname_v_csv, fname_ioapi_latlon)
    plt.gcf().savefig('vulcan_csv_ioapi_latlon.png')
