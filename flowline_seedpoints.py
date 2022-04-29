#! /usr/bin/env python

"""
calculate flowlines on an ice sheet, based on maps of the flowfield in x and y direction
This piece of code reads a line shapefile and creates seedpoints from the lines. Along
all lines points are created at given intervals.
"""

import os
import numpy as np
import pandas as pd
import rasterio
import platform
import geopandas as gpd
from shapely.ops import unary_union

import flowline_funcs as flf

osys = platform.system()
print('operating system is %s' %osys)

if osys != 'Windows':
    seedfile = '/home/horstm/erc/vel_greenland_crop/something.shp'
    outfolder = '/home/horstm/erc/vel_greenland_crop_processed/'
else:
    #seedfile = r'N:/MODIS/gis/Greenland2400mContours.shp'  # needs to be a line shapefile, can contain many lines
    seedfile = r'N:/MODIS/polygons/seedline_v2.shp'
    outfolder = r'N:/MODIS/polygons/'

distance_delta = 15000  # distance of seedpoints located along seedfile polylines

# -------------------------------------------------------------------------------------------------
# check if output folder exists, if no create
isdir = os.path.isdir(outfolder)
if not isdir:
    os.mkdir(outfolder)

# ------------------------------------ read seedfile data ------------------------------------------
cont = gpd.read_file(seedfile)
#cont_sel = cont.loc[cont['ELEV'] == 2400]  # select contour linesyb elevation. Given contourline needs to exist.
cont['length'] = cont['geometry'].length  # determine length of the contour lines
seedlines = cont.loc[cont['length'] > 21000]  # remove very short contour lines
seedlines['geometry'] = seedlines['geometry'].simplify(1000)  # simplyfy the lines

# create the seedpoints
gdf_out = gpd.GeoDataFrame()
gdf_out['geometry'] = None
# Set the GeoDataFrame's coordinate system according to input
gdf_out.crs = cont.crs

store_ids = []
store_points = []
id_counter = 0
for index, row in seedlines.iterrows():
    distances = np.arange(0, row['geometry'].length, distance_delta)
    points = [row['geometry'].interpolate(distance) for distance in distances] + [row['geometry'].boundary[1]]
    points_id = np.arange(id_counter, len(points))
    store_ids.extend(points_id)
    store_points.extend(points)
    id_counter += len(points)

gdf_new = gpd.GeoDataFrame({'id': store_ids}, geometry=store_points, crs=cont.crs)
gdf_new.to_file(outfolder + 'seedpoints_v3.shp')
