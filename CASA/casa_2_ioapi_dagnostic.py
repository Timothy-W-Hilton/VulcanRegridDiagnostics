import netCDF4
import os
import matplotlib.pyplot as plt
import numpy as np
from stem_pytools.domain import find_nearest_stem_xy

def get_casa_raw(pt_lat=37, pt_lon=-121.25):
    nc = netCDF4.Dataset("/project/projectdirs/m2319/Data/CASA/GEE.3hrly.1x1.25.2015.nc")
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    lon, lat = np.meshgrid(lon, lat)
    x, y = find_nearest_stem_xy(pt_lon, pt_lat, lon, lat)
    pt_gee = nc.variables['GEE'][:, x, y]
    nc.close()
    return x, y, lon[x, y], lat[x, y], pt_gee

def get_casa_raw_ioapi(pt_lat=37, pt_lon=-121.25):
    nc = netCDF4.Dataset('/global/homes/t/twhilton/Code/Regrid/VulcanRegrid/CASA/casa_raw_grid.nc')
    lon = nc.variables['LON'][0, 0, ...]
    lat = nc.variables['LAT'][0, 0, ...]
    nc.close()
    x, y = find_nearest_stem_xy(pt_lon, pt_lat, lon, lat)
    nc = netCDF4.Dataset('/global/homes/t/twhilton/Code/Regrid/VulcanRegrid/CASA/tim_test3.nc')
    pt_gee = nc.variables['CO2_FLUX'][:, 0, x, y]
    return x, y, lon[x, y], lat[x, y], pt_gee

if __name__ == "__main__":
    x_raw, y_raw, lon_raw, lat_raw, gee_raw = get_casa_raw()
    x_ioapi, y_ioapi, lon_ioapi, lat_ioapi, gee_ioapi = get_casa_raw_ioapi()

    fig, ax = plt.subplots(nrows=2, ncols=1)
    ax[0].plot(gee_raw)
    ax[0].set_title('CASA raw ({}, {}) ({}, {})'.format(
        x_raw, y_raw, lon_raw, lat_raw))
    ax[1].plot(gee_ioapi)
    ax[1].set_title('CASA native grid I/O API ({}, {}) ({}, {})'.format(
        x_ioapi, y_ioapi, lon_ioapi, lat_ioapi))
    for this_ax in ax:
        this_ax.set_ylabel('GEE')
    ax[1].set_xlabel('tstep')
