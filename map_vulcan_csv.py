import os
import numpy as np
import pandas as pd
import netCDF4
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


def parse_ioapi_latlon(fname_ioapi_latlon_output):
    """ read longitude, latitude data from I/O API latlon
    """
    nc = netCDF4.Dataset(fname_ioapi_latlon_output)
    lat = nc.variables['LAT'][0, 0, ...]
    lon = nc.variables['LON'][0, 0, ...]
    nc.close()
    return lon, lat


def vulcan_csv_2_latslons(fname_v_csv):
    """place longitude, latitude data from Vulcan CSV into 2-D arrays

    Data frame df must have columns i, j,
    ddx, ddy, with i and j containing the array indices, and ddx and
    ddy containing longitude and latitude, respectively.
    """
    df = pd.read_csv(fname_v_csv)
    imax = df['i'].max()
    jmax = df['j'].max()
    lons = np.ndarray((imax, jmax))
    lats = np.ndarray((imax, jmax))
    for this_row in df.itertuples():
        lons[this_row.i - 1, this_row.j - 1] = this_row.ddx
        lats[this_row.i - 1, this_row.j - 1] = this_row.ddy
    return (lons, lats)


def draw_map(v_lons, v_lats, i_lons, i_lats):
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
    m_csv.scatter(v_lons, v_lats,
                  latlon=True,
                  marker='o', color='black',
                  label='Vulcan CSV')
    m_csv.scatter(i_lons, i_lats,
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
    v_lons, v_lats = vulcan_csv_2_latslons(fname_v_csv)
    i_lons, i_lats = parse_ioapi_latlon(fname_ioapi_latlon)
    m_csv = draw_map(v_lons, v_lats, i_lons, i_lats)
    plt.gcf().savefig('vulcan_csv_ioapi_latlon.png')
