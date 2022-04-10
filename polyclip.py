#! /usr/bin/env python

"""
Find the slush limit in MODIS scenes. Based on the approach by Greuell and Knap (2000),
originally developed and applied to AVHRR scenes.

Testing possibilities to clip from a raster using a series of polygons

"""

import numpy as np
import pandas as pd
import os
import glob
import rasterio
import rasterio.mask
import fiona
import xarray as xr

import geopandas as gpd
# from shapely.geometry import MultiLineString
# from shapely.geometry import LineString

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

####################################################################################################################
inraster = r'N:/MODIS/mask_greenland_icesheet/dem_test.tif'
infolder = r'N:/MODIS/gis/test_poly/'
outfolder = r'N:/MODIS/gis/test_poly_out/'

infiletype = '.shp'

crs = 'EPSG:3413' # coordination system
# -------------------------------------------------------------------------------------------------
# check if output folder exists, if no create
isdir = os.path.isdir(outfolder)
if not isdir:
    os.mkdir(outfolder)

files_o = glob.glob(infolder + '*' + infiletype)  # create list of all relevant filenames in folder
print('found ' + str(len(files_o)) + ' shapefiles in ' + infolder)

dem = rasterio.open(inraster)  # read the ice sheet dem and create mask
CRS = dem.crs.data['init']

for ni, i in enumerate(files_o):
    #ipoly = gpd.read_file(i)
    with fiona.open(i, "r") as shapefile:
        ipoly = [feature["geometry"] for feature in shapefile]
    clip_dem, out_transform = rasterio.mask.mask(dem, ipoly, crop=True)
    clip_meta = dem.meta
    clip_meta.update({"driver": "GTiff",
        "height": clip_dem.shape[1],
        "width": clip_dem.shape[2],
        "transform": out_transform})
    clip_dem = np.where(clip_dem == dem.nodatavals, np.NaN, clip_dem)  # set all nodata points to NaN
    # write output raster (only development)
    
    # extract coordinates
    x_coords = np.arange(out_transform[2], out_transform[0] * clip_dem.shape[2] +
                         out_transform[2], out_transform[0]) - out_transform[4]/2
    y_coords = np.arange(out_transform[5], out_transform[4] * clip_dem.shape[1] +
                         out_transform[5], out_transform[4]) + out_transform[4]/2

    coord_mm = [np.min(x_coords), np.max(x_coords), np.min(y_coords), np.max(y_coords)]

    # remove one dimension from DEM
    clip_dem = clip_dem.squeeze()
    cdem = xr.Dataset({'dem': (('y', 'x'), clip_dem)}, coords={'y': y_coords, 'x': x_coords})

    with rasterio.open(outfolder + 'test'+str(ni)+'.tif', "w", **clip_meta) as dest:
        dest.write(cdem.to_array())
print('')
