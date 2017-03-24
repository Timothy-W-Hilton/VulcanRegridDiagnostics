"""locate CASA grid cell nearest (lon, lat) point and plot its GEE

locates the nearest cell to an arbitrary (lon, lat) point.  Reads the
GEE time series for both the raw CASA data (from NACP site) and the
non-regridded I/O API CASA data file.  Plots a two-panel plot with
time series on identical axes.
"""

import netCDF4
import os
import matplotlib.pyplot as plt
import numpy as np
from stem_pytools.domain import find_nearest_stem_xy


def get_casa_raw(pt_lat=37, pt_lon=-121.25):
    """read CASA raw data, find grid center nearest to arguments

    ARGS:
    pt_lat: float; latitude coordinate of the point to locate
    pt_lon: float; longitude coordinate of the point to locate

    RETURNS:
    x (integer): x index of CASA grid cell nearest (pt_lat, pt_lon)
    y (integer): y index of CASA grid cell nearest (pt_lat, pt_lon)
    lon[x, y] (float): the longitude of [x, y] in the CASA grid
    lat[x, y] (float): the latitude of [x, y] in the CASA grid
    pt_gee (np.ndarray): CASA GEE time series for [x, y]
    """
    nc = netCDF4.Dataset(os.path.join("/", "project", "projectdirs",
                                      "m2319", "Data", "CASA",
                                      "GEE.3hrly.1x1.25.2015.nc"))
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    lon, lat = np.meshgrid(lon, lat)
    x, y = find_nearest_stem_xy(pt_lon, pt_lat, lon, lat)
    pt_gee = nc.variables['GEE'][:, x, y]
    nc.close()
    return x, y, lon[x, y], lat[x, y], pt_gee


def get_casa_ioapi(fname_latlon, fname_gee, pt_lat=37, pt_lon=-121.25):
    """read CASA I/O API data, find grid center nearest to arguments

    This expects the non-regridded CASA data in I/O API format.  The
    longitude and latitude coordinates for the I/O API CASA grid are
    read from the output of latlon
    (https://www.cmascenter.org/ioapi/documentation/all_versions/html/LATLON.html).

    ARGS:
    fname_latlon (string): full path to latlon output
    fname_gee (string): full path to the I/O API file containing GEE.
        The GEE variable must be CO2_FLUX
    pt_lat: float; latitude coordinate of the point to locate
    pt_lon: float; longitude coordinate of the point to locate

    RETURNS:
    x (integer): x index of CASA I/O API grid cell nearest (pt_lat, pt_lon)
    y (integer): y index of CASA I/O API grid cell nearest (pt_lat, pt_lon)
    lon[x, y] (float): the longitude of [x, y] in the CASA I/O API grid
    lat[x, y] (float): the latitude of [x, y] in the CASA I/O API grid
    pt_gee (np.ndarray): CASA I/O API GEE time series for [x, y]
    """
    nc = netCDF4.Dataset(fname_latlon)
    lon = nc.variables['LON'][0, 0, ...]
    lat = nc.variables['LAT'][0, 0, ...]
    nc.close()
    x, y = find_nearest_stem_xy(pt_lon, pt_lat, lon, lat)
    # x, y = (24, 28)
    nc = netCDF4.Dataset(fname_gee)
    pt_gee = nc.variables['CO2_FLUX'][:, 0, x, y]
    return x, y, lon[x, y], lat[x, y], pt_gee


if __name__ == "__main__":
    x_raw, y_raw, lon_raw, lat_raw, gee_raw = get_casa_raw()
    latlon_casa_raw = os.path.join('/', 'global', 'homes', 't',
                                   'twhilton', 'Code', 'Regrid',
                                   'VulcanRegrid', 'CASA',
                                   'casa_raw_grid.nc')
    gee_casa_raw_ioapi = os.path.join('/global', 'homes', 't',
                                      'twhilton', 'Code', 'Regrid',
                                      'VulcanRegrid', 'CASA',
                                      'tim_test3.nc')
    latlon_casa_9km = os.path.join('/', 'global', 'homes', 't',
                                   'twhilton', 'Code', 'Regrid',
                                   'VulcanRegrid', 'CASA',
                                   'casa_grid.nc')
    gee_casa_9km = os.path.join("/", "project", "projectdirs",
                                "m2319", "transfer",
                                "CASA_GEE_perm2_stem9km_ioapi.nc")
    (x_ioapi, y_ioapi, lon_ioapi,
     lat_ioapi, gee_ioapi) = get_casa_ioapi(latlon_casa_raw,
                                            gee_casa_raw_ioapi)
    (x_9km, y_9km, lon_9km,
     lat_9km, gee_9km) = get_casa_ioapi(latlon_casa_9km,
                                            gee_casa_9km)
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(8, 16))
    ax[0].plot(gee_raw)
    ax[0].set_title('CASA raw ({}, {}) ({}, {})'.format(
        x_raw, y_raw, lon_raw, lat_raw))
    ax[1].plot(gee_ioapi)
    ax[1].set_title('CASA native grid I/O API ({}, {}) ({}, {})'.format(
        x_ioapi, y_ioapi, lon_ioapi, lat_ioapi))
    ax[2].plot(gee_ioapi)
    ax[2].set_title('CASA 9 KM I/O API ({}, {}) ({}, {})'.format(
        x_9km, y_9km, lon_9km, lat_9km))
    for this_ax in ax:
        this_ax.set_ylabel('GEE')
    ax[0].set_ylim((0, -3.0e-7))
    ax[1].set_ylim((0, -3.0e11))
    ax[2].set_ylim((0, -3.0e11))
    ax[2].set_xlabel('tstep')
    fname = '{}_{}.pdf'.format(x_ioapi, y_ioapi)
    print 'saving {}'.format(fname)
    fig.savefig(fname)
