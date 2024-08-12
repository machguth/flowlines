#! /usr/bin/env python

"""
crop a geoTIFF using gdalwarp
"""

import os
import subprocess
import glob
import re
import platform

osys = platform.system()
print('operating system is %s' %osys)

#resampling = 'near'
resampling = 'cubic'

# xmax = -104000  # 5000
# xmin = -456000  # -278000
# ymin = -2050000  # -3125000  # ymin and ymax were exchanged during the last processing
# ymax = -1425000  # -2050000
# xmax = 562000  # 5000
# xmin = 168000  # -278000
# ymin = -1193000  # -3125000  # ymin and ymax were exchanged during the last processing
# ymax = -1058000  # -2050000
# xmax = 618000  # NE
# xmin = 145000  #
# ymin = -1341000  #
# ymax = -1017000  #
# xmax = 3000  #
# xmin = -318000  #
# ymin = -2741000  #
# ymax = -2320000  #
xmax = 665000  #
xmin = -601500  #
ymin = -3256000  #
ymax = -826000  #

s_srs = 'EPSG:3413'  # source projection
t_srs = 'EPSG:3413'  # target projection

tr_x = 500
tr_y = 500

# input can be a folder or a single file. Single files recognized by the '.' in file name
if osys != 'Windows':
    # ingrid = '/home/horstm/erc/vel_greenland_raw/'
    ingrid = '/home/horstm/ownCloud/GIS/dem_data/dem_arctic/arcticdem_mosaic_100m_v30_greenland_icesheet_geoidCorr_GapFilled.tif'
else:
    #ingrid = r'N:/MODIS/vel_greenland_raw/'
    #ingrid = r'N:/MODIS/vel_greenland_processed/'
    ingrid = r'N:/MODIS/ArcticDEM_100m/arcticdem_mosaic_100m_v30_greenland_icesheet_GeoidCorr_GapFilled.tif'
    # ingrid = r'C:/Users/machguth/switchdrive/GIS/dem_data/dem_arctic/arcticdem_mosaic_100m_v30_greenland_icesheet_geoidCorr.tif'
    #ingrid = r'N:/MODIS/ArcticDEM_100m/arcticdem_mosaic_100m_v30_greenland_icesheet_geoidCorr.tif'

filetype = '.tif'  # not relevant if single file

# output can be a folder or a single file
if osys != 'Windows':
    # outgrid = '/home/horstm/erc/vel_greenland_crop/'
    outgrid = '/home/horstm/erc/mask_greenland_icesheet/dem_southwest_test.tif'
else:
    #outgrid = r'N:/MODIS/vel_greenland_500m/'
    outgrid = r'N:/MODIS/ArcticDEM_500m/arcticdem_mosaic_500m_v30_greenland_icesheet_geoidCorr_GapFilled.tif'
    # outgrid = r'E:/MODIS/mask_greenland_icesheet/dem_northwest.tif'
    #outgrid = r'N:/MODIS/mask_greenland_icesheet/dem_NE.tif'
    # outgrid = r'E:/MODIS/mask_greenland_icesheet/test.tif'

# -------------------------------------------------------------------------------------------------
test = ingrid.split('.')  # check whether there is a dot in ingrid. If not, it is a folder

if len(test) == 1:
    folder = True
    isdir = os.path.isdir(outgrid)
    if not isdir:
        os.mkdir(outgrid)
    files = glob.glob(ingrid + '*' + filetype)  # create list of all relevant filenames in folder
else:
    folder = False
    files = [ingrid]
# -------------------------------------------------------------------------------------------------

grid = '-s_srs "' + s_srs + '" -t_srs "' + t_srs + '" -te ' \
       + str(xmin) + ' ' + str(ymin) + ' ' + str(xmax) + ' ' + str(ymax) \
       + ' -tr ' + str(tr_x) + ' ' + str(tr_y)

for f in files:

    f_el = re.split(r'[./\\]', f)

    if folder:
        outfile = outgrid + f_el[-2] + '-crop.' + f_el[-1]
    else:
        outfile = outgrid

    # Crop using gdalwarp
    cmd = 'gdalwarp -r ' + resampling + ' -of GTiff ' + grid + ' ' + f + ' ' + outfile

    print('** CMD ** ' + cmd)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, shell=True)
    (stdout, stderr) = p.communicate()
    p.wait()
    if p.returncode != 0:
        print(f + ': gdalwarp failed: ' + str(stderr))
