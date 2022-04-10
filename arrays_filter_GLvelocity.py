#! /usr/bin/env python

"""
stack a set of arrays and filter them
"""

import os
import numpy as np
import pandas as pd
import glob
import rasterio
import xarray as xr
import platform

osys = platform.system()
print('operating system is %s' %osys)

if osys != 'Windows':
    infolder = '/home/horstm/erc/vel_greenland_crop/'
    outfolder = '/home/horstm/erc/vel_greenland_crop_processed/'
else:
    #infolder = r'N:/MODIS/vel_greenland_crop2/'
    infolder = r'N:/MODIS/vel_greenland_raw/'
    #outfolder = r'N:/MODIS/vel_greenland_crop_processed_test2/'
    outfolder = r'N:/MODIS/vel_greenland_processed/'

filetype = '.tif'

# dictionary which defines which output layer (value) gets filtered how (key)
vdict = {"ex": "vx",
         "ey": "vy"}

# original no_data value and new no_data value
# this code uses argmin() and an extreme negative no_data will be selected by argmin() over valid values.
nd = [-2000000000, 2000000]

# -------------------------------------------------------------------------------------------------
# check if output folder exists, if no create
isdir = os.path.isdir(outfolder)
if not isdir:
    os.mkdir(outfolder)

# -------------------------------------------------------------------------------------------------
files_o = glob.glob(infolder + '*' + filetype) # create list of all relevant file names in folder

print('found ' + str(len(files_o)) + ' files in ' + infolder)

files_i = []
for i in files_o:  # split file names by underscores to extract information on content and dates
    ci = i.split('_')
    files_i.append(ci)

files_np = np.asarray(files_i)
years = np.unique(files_np[:, -3])
vt = np.unique(files_np[:, -2])  # get all different types of arrays

sorted = pd.DataFrame(columns=vt, index=years)  # create dataframe

for vi, v in enumerate(vt):  # fill dataframe in ascending order with file names
    vw = np.where(files_np[:, -2] == v)
    for yi, y in enumerate(years):
        yw = np.where(files_np[vw[0], -3] == y)
        sorted[v][y] = files_o[int(vw[0][yw[0]])]

# check size of arrays (assuming all arrays are of identical size)
raster = rasterio.open(sorted['ex']['2015'])

stack = xr.DataArray(np.zeros((raster.height, raster.width, len(years), len(vt))),
                     dims=['x', 'y', 'z', 'v'], coords={"v": vt})

for vi, v in enumerate(vt):
    for rowi, row in enumerate(sorted[v]):  # open the relevant files and stack them
        raster = rasterio.open(row)
        rasterdata = raster.read(1)
        nd_ix = np.where(rasterdata == nd[0])  # find no_data
        rasterdata[nd_ix] = nd[1]              # assign new no_data value
        stack[:, :, rowi, vi] = rasterdata

ix_min_of_stack = stack.argmin('z')  # find the index of the minimum value along axis 2 (z-axis)
ix_min_of_stack = ix_min_of_stack.astype(np.int32)  # to 32 bit integer, otherwise issues when writing int64 to geoTIFF

# prepare output array
v_filtered = xr.DataArray(np.zeros((raster.height, raster.width, len(vdict))),
                          dims=['x', 'y', 'v_filt'], coords={"v_filt": list(vdict.values())})

# get only the optimal data (lowest error) for 'vx' and 'vy'
for ii, i in enumerate(vdict):
    v_filtered[:, :, ii] = stack.sel(z=ix_min_of_stack.sel(v=i), v=vdict[i])
    v_filtered = v_filtered.astype(np.float32)

# write output as geoTIFF
for i in vdict:
    # write both velocity composites to output
    filtered_dataset = rasterio.open(outfolder +'greenland_vel_mosaic200_2015-2018_'+vdict[i]+'_v02-composite-crop.tif',
            'w', driver='GTiff', height=v_filtered.shape[0], width=v_filtered.shape[1],  count=1,
            dtype=v_filtered.dtype, crs=raster.crs,  nodata=nd[1], transform=raster.transform)

    print('writing: ' + outfolder +'greenland_vel_mosaic200_2015-2018_'+vdict[i]+'_v02-composite-crop.tif')
    filtered_dataset.write(v_filtered.sel(v_filt=vdict[i]), 1)
    filtered_dataset.close()

    # write the arrays to output which indicate the year that was used (year of optimal quality) for all pixels
    filtered_dataset = rasterio.open(outfolder +'greenland_vel_mosaic200_2015-2018_'+i+'_v02-composite-crop.tif',
            'w', driver='GTiff', height=ix_min_of_stack.shape[0], width=ix_min_of_stack.shape[1],  count=1,
            dtype=ix_min_of_stack.dtype, crs=raster.crs,  nodata=nd[1], transform=raster.transform)

    print('writing: ' + outfolder +'greenland_vel_mosaic200_2015-2018_'+i+'_v02-composite-crop.tif')
    filtered_dataset.write(ix_min_of_stack.sel(v=i), 1)
    filtered_dataset.close()

#sorted.to_excel('/home/horstm/erc/vel_greenland_raw/sorted.xlsx')
