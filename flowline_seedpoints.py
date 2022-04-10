#! /usr/bin/env python

"""
calculate flowlines on an ice sheet, based on maps of the flowfield in x and y direction
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
    infolder = '/home/horstm/erc/vel_greenland_crop/'
    outfolder = '/home/horstm/erc/vel_greenland_crop_processed/'
else:
    demmask = r'N:/MODIS/mask_greenland_icesheet/dem_test.tif'
    seedfile = r'N:/MODIS/gis/Greenland2400mContours.shp'  # needs to be a line shapefile, can contain many lines
    outfolder = r'N:/MODIS/gis/'

distance_delta = 5000  # distance of seedpoints located along seedfile polylines

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
seedlines['geometry'] = seedlines['geometry'].simplify(10000)  # vers strongly simplyfy the lines

# create the seedpoints
gdf_out = gpd.GeoDataFrame()
gdf_out['geometry'] = None
# Set the GeoDataFrame's coordinate system according to input
gdf_out.crs = cont.crs

for index, row in seedlines.iterrows():
    distances = np.arange(0, row['geometry'].length, distance_delta)
    points = [row['geometry'].interpolate(distance) for distance in distances] + [row['geometry'].boundary[1]]
    multipoint = unary_union(points)  # or new_line = LineString(points)
    gdf_out.loc[index, 'geometry'] = multipoint

# write output
gdf_out.to_file(outfolder + 'seedpoints5000.shp')



