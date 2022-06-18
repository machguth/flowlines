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

import flowline_funcs as flf

osys = platform.system()
print('operating system is %s' %osys)

if osys != 'Windows':
    infolder = '/home/horstm/erc/vel_greenland_crop/'
    outfolder = '/home/horstm/erc/vel_greenland_crop_processed/'
else:
    # infilex = r'N:/MODIS/vel_greenland_crop_processed_test2/greenland_vel_mosaic200_2015-2018_vx_v02-composite-crop.tif'
    # infiley = r'N:/MODIS/vel_greenland_crop_processed_test2/greenland_vel_mosaic200_2015-2018_vy_v02-composite-crop.tif'
    # demmask = r'N:/MODIS/mask_greenland_icesheet/dem_test.tif'
    infilex = r'N:/MODIS/vel_greenland_500m/greenland_vel_mosaic500_2015-2018_vx_v02-composite-crop.tif'
    infiley = r'N:/MODIS/vel_greenland_500m/greenland_vel_mosaic500_2015-2018_vy_v02-composite-crop.tif'
    demmask = r'N:/MODIS/ArcticDEM_500m/arcticdem_mosaic_500m_v30_greenland_icesheet_geoidCorr.tif'
    seedfile = r'N:/MODIS/polygons/seedpoints_v3.2.shp'  # needs to be a point shapefile

    outfolder = r'N:/MODIS/polygons/'

version = 'v3.2test'  # Simple identifier, has no other function than being appended to the end of the filename

vmin = 3  # [m yr-1] minimum flow speed. Flowlines are ended when they reach areas of v < vmin
buff = 10000  # [m] buffer distance by which the flowlines get buffered (to create polygons from the flow lines)
flminlength = 30 # minimum required points in a flowline

# whenever a seedfile is specified, seedpoints are read from the seedfile
# if seedfile is not specified, then seed points get created along a vertical line as defined below
seedXcoord = -5050  # in coordinates of CRS
seedspacing = 15000

# -------------------------------------------------------------------------------------------------
# check if output folder exists, if no create
isdir = os.path.isdir(outfolder)
if not isdir:
    os.mkdir(outfolder)

# ------------------------------------ analyse/preparing the DEM data -----------------------------------------
demDR = rasterio.open(demmask)  # read the ice sheet dem and create mask
demDRnd = demDR.nodatavals  # get nodata value
CRS = demDR.crs.data['init']

# extract coordinates
x_coords = np.arange((demDR.get_transform())[0], (demDR.get_transform())[1] * demDR.shape[1] +
                     (demDR.get_transform())[0], (demDR.get_transform())[1]) - (demDR.get_transform())[5] / 2
y_coords = np.arange((demDR.get_transform())[3], (demDR.get_transform())[5] * demDR.shape[0] +
                     (demDR.get_transform())[3], (demDR.get_transform())[5]) + (demDR.get_transform())[5] / 2

coord_mm = [np.min(x_coords), np.max(x_coords), np.min(y_coords), np.max(y_coords)]

dem = demDR.read(1)  # read values
dem = np.where(dem == demDR.nodatavals, np.NaN, dem)

cs = (demDR.get_transform())[1]

# ------------------------ seed points: read from seed file or create ------------------------------------------
if 'seedfile' in locals():
    seedpoints = gpd.read_file(seedfile)

else:  # do not create seedpoints from seedfile but create points along a vertical line
    ycoords = np.arange(coord_mm[2]+buff+1, coord_mm[3]-buff-1, seedspacing)
    seedpoints = gpd.GeoDataFrame({'id': np.arange(0, len(ycoords))},
                                  geometry=gpd.points_from_xy([seedXcoord] * len(ycoords),ycoords))

# ------------------------------ read x and y velocity values -----------------------------------
xvaDR = rasterio.open(infilex)
yvaDR = rasterio.open(infiley)
xva = xvaDR.read(1)
yva = yvaDR.read(1)
xva = np.where(xva == xvaDR.nodatavals, np.NaN, xva)  # replace noData value with np.nan
yva = np.where(yva == yvaDR.nodatavals, np.NaN, yva)  # replace noData value with np.nan

if (xvaDR.get_transform())[0] != (demDR.get_transform())[0] or \
        (xvaDR.get_transform())[3] != (demDR.get_transform())[3]:
    print('')
    print('! Coordinate references of DEM mask and velocity field not identical. !')
    print('! Either x or y reference differ (upper/left vs. lower/right margin). !')
    print('xvaDR.get_transform():', xvaDR.get_transform())
    print('demDR.get_transform():', demDR.get_transform())
    quit()

# mask veloctiy fields to ice sheet extent
xva = np.where(np.isnan(dem), np.NaN, xva)
yva = np.where(np.isnan(dem), np.NaN, yva)

# Initiate DataFrame to store data on all flowlines
df_fl = pd.DataFrame(columns=['ID', 'n', 'X', 'Y', 'Z', 'v', 'd'])

# ***************************************** calculate the flowlines *******************************************

#for p in seedpoints:
for index, row in seedpoints.iterrows():

    print('starting with seedpoint: ', index)

    p = (row.geometry.x, row.geometry.y)  # create x/y tuple. a bit complicated but for compatibility to code below

    fl_id = index  # get ID of the flowline

    # ----------------------------- preparations --------------------------
    # initialize flowline and starting point
    flowline = []

    # initialize conflict tracking
    conflict = [False, 0]  # 1st element indicates if conflict (True = conflict), 2nd count of consecutive conflicts
    max_conflict = 3  # maximum count of consecutive conflicts. If exceeded then flowline is ended

    # initialize arrays of array-coordinates
    # Calculation of the flowline from one point of the flowline (not identical to grid cell coordinates) to the next
    # basically requires only information on one grid cell. However, there are cases of conflict, when the flowline
    # leaves grid cell 'n' (old) into direction fo grid cell 'n + 1' (new) but flowdirection of grid cell 'n + 1'
    # immediately tries to send back the flowline into grid cell 'n'. In this case the average flowdirection between
    # two grid cells needs to be calculated.
    x, y = [np.nan, np.nan], [np.nan, np.nan]

    # initialize array of the margin of the grid cells. As above, basically only the margins of one grid cell are needed
    # but in case of conflict the margins of the two conflicting grid cess are needed.
    m = np.empty((2, 4))
    m[:] = np.nan

    # Boolean variable: as long as True the flowline continues
    # Gets set to True if seedpoint passes initial check
    c = False

    # the algorithm keeps track of two grid cells, the one one that is currently being crossed (= new cell)
    # and the one that has just been crossed (=old cell). Their indices in a two element array are 0 (=new) and 1 (=old)
    # Typically the old grid cell is of no relevance, unless there is a conflict between flow directions between old
    # and the new cell. In that case it is possible that the line never enters the new cell but gets immediately
    # redirected into the old cell. Only in that case the old cell is used twice in a row.
    i = 0  # Start with 0.

    # ----------------------------- computation of flowlines ------------------------------
    # Check if there actually is a cell, if there are data and if velocity high enough

    x[0], y[0] = int(np.floor(p[0] - coord_mm[0]) / cs), \
                 int(np.floor(coord_mm[3] - p[1]) / cs)  # array coords. initial cell

    if 0 <= x[0] < dem.shape[1] and 0 <= y[0] < dem.shape[0]:
        if np.isfinite(dem[y[0], x[0]]):
            v = np.sqrt(xva[y[0], x[0]]**2 + yva[y[0], x[0]]**2)
            if v > vmin:
                c = True

    flowline.append([fl_id, 0, p[0], p[1], dem[y[0], x[0]], v, 0])  # information for initial point

    # ****************** main loop ***********************
    while c:

        # define margins of current cells
        m[0, :] = [x[0] * cs + coord_mm[0], (x[0] + 1) * cs + coord_mm[0],
                   coord_mm[3] - y[0] * cs, coord_mm[3] - (y[0] + 1) * cs]
        m[1, :] = [x[1] * cs + coord_mm[0], (x[1] + 1) * cs + coord_mm[0],
                   coord_mm[3] - y[1] * cs, coord_mm[3] - (y[1] + 1) * cs]

        # *** conflict solving ***
        # If conflict, then analyse whether the old or new cell have stronger flow (in x or y direction,
        # depending on which border the flowline is about to cross). If, for example, the new cell has stronger flow
        # (np.abs(yva[y[0], x[0]]) > np.abs(yva[y[1], x[1]]) or np.abs(xva[y[0], x[0]]) > np.abs(xva[y[1], x[1]]))
        # then the flowline gets pushed back into the old cell. Flowdirection and strength of the old cell, however,
        # cannot be used in this case because this would again send the flowline into the new cell. The conflict
        # is solved by using the average x and y flow directions of old and new cell.
        if conflict[0]:
            if x[0] == x[1]:  # means there is a conflict between two grid cells with differing y coords
                # which of the two conflicting cells has stronger flow in y direction?
                if np.abs(yva[y[0], x[0]]) > np.abs(yva[y[1], x[1]]):
                    i = 1  # flowline will move through second (=old) cell in array of current cells
                else:
                    i = 0  # flowline will move through first (=new) cell in array of current cells
            else:  # means there is a conflict between two grid cells with differing x coords
                # which of the two conflicting cells has stronger flow in x direction?
                if np.abs(xva[y[0], x[0]]) > np.abs(xva[y[1], x[1]]):
                    i = 1  # flowline will move through second (=old) cell in array of current cells
                else:
                    i = 0  # flowline will move through first (=new) cell in array of current cells
            # If there is a conflict, regardless of in x or y direction, average x and y flow directions are used
            xv, yv = np.mean([xva[y[0], x[0]], xva[y[1], x[1]]]), np.mean([yva[y[0], x[0]], yva[y[1], x[1]]])
        else:
            i = 0
            xv, yv = xva[y[0], x[0]], yva[y[0], x[0]]

        # calculate where the flowline leaves the current cell
        s = yv / xv  # slope of flow direction

        # distances to cell margins: d = x or y distance, l = absolute length, h = horizontal, v = vertical,
        # i = intercept
        if xv < 0:  # m[i, 0] relevant, heading towards the more negative (western) margin
            hd = m[i, 0] - p[0]
            yih = hd * s
            lh = np.sqrt(yih**2 + hd**2)
        else:  # m[i, 1] relevant, heading towards the less negative (eastern) margin
            hd = m[i, 1] - p[0]
            yih = hd * s
            lh = np.sqrt(yih**2 + hd**2)
        if yv < 0:  # m[i, 3] relevant, heading towards the more negative (southern) margin
            vd = m[i, 3] - p[1]
            xiv = vd / s
            lv = np.sqrt(xiv**2 + vd**2)
        else:  # m[i, 2] relevant, heading towards the less negative (northern) margin
            vd = m[i, 2] - p[1]
            xiv = vd / s
            lv = np.sqrt(xiv**2 + vd**2)

        if hd == 0:
            print('!!! hd = 0')
            break
        if vd == 0:
            print('!!! vd = 0')
            break

        # move information of current grid cell to second rows of x, y and m arrays
        x[1], y[1], m[1, :] = x[i], y[i], m[i, :]

        # define the next (=new) grid cell based on whether lv or lh is smaller
        # information on current (=old) grid cell gets overwritten in x[0], y[0] but is still stored in x[1], y[1]
        if lh < lv:
            p = [p[0] + hd, p[1] + yih]
            x[0], y[0] = x[i] + int(np.round(np.abs(hd)/hd)), y[i]
        else:
            p = [p[0] + xiv, p[1] + vd]
            x[0], y[0] = x[i], y[i] - int(np.round(np.abs(vd)/vd))

        # calculate current length of flowline and velocity
        if lh < lv:
            l = lh
        else:
            l = lv
        v = np.sqrt(xv**2 + yv**2)

        # append information on current (=old) point to flowline
        flowline.append([fl_id, len(flowline), p[0], p[1], dem[y[1], x[1]], v, l])

        # *** Conflict solving ****
        # check whether there is a conflict that needs to be solved
        # in the next step (computation of next flowline segment)
        if np.isfinite(x[0]) and np.isfinite(x[1]):  # skip checking during initial iteration
            if lh < lv:
                if np.sign(xva[y[0], x[0]]) == np.sign(xva[y[1], x[1]]):
                    conflict[0] = False
                    conflict[1] = 0
                else:
                    conflict[0] = True
                    conflict[1] += 1
            else:
                if np.sign(yva[y[0], x[0]]) == np.sign(yva[y[1], x[1]]):
                    conflict[0] = False
                    conflict[1] = 0
                else:
                    conflict[0] = True
                    conflict[1] += 1

        # check whether conditions are given to proceed to the calculation of the next segment
        # or if flowline needs to be ended
        if 0 <= x[0] < dem.shape[1] and 0 <= y[0] < dem.shape[0]:
            if np.isfinite(dem[y[0], x[0]]):
                # calculate local flowspeed
                v = np.sqrt(xva[y[0], x[0]]**2 + yva[y[0], x[0]]**2)
                if v > vmin:
                    c = True
                else:
                    c = False
            else:
                c = False
        else:
            c = False

        if conflict[1] > max_conflict:
            c = False

    df_fl = pd.concat([pd.DataFrame(np.asarray(flowline), columns=df_fl.columns), df_fl], ignore_index=True)

# create the polygons but also write the flowlines to output
flf.create_polygons(df_fl, CRS, buff, outfolder, version, flminlength)

# create map with the flowlines
flf.flowlines_plot(df_fl, dem, x_coords, y_coords, CRS, outfolder, version)
