#! /usr/bin/env python

"""
calculate flowlines on an ice sheet, based on maps of the flowfield in x and y direction
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon


# #####################################################################################################################
def create_polygons(df_fl, CRS, buff, outfolder, version, flminlength):  # create the polygons

    # check minimum length is at least set to 2, otherwise LineString will issue errors
    if flminlength == 1:
        flminlength = 2

    # convert to LineString
    flid = np.unique(df_fl['ID'])
    # prepare GeoDataFrame to store final polygons
    gdf_poly = gpd.GeoDataFrame()
    gdf_fl = gpd.GeoDataFrame()
    gdf_poly['geometry'] = None
    gdf_fl['geometry'] = None
    # Set the GeoDataFrames' coordinate system according to the DEM
    gdf_poly.crs = CRS
    gdf_fl.crs = CRS

    for f in flid:
        flt = df_fl.loc[df_fl['ID'] == f]
        geo = [xy for xy in zip(flt.X, flt.Y)]
        if len(geo) > flminlength:  # minimum lentgh of a flowline
            s = LineString(geo)
            final = s.buffer(buff, join_style=2)
            gdf_poly.loc[f, 'geometry'] = final
            gdf_fl.loc[f, 'geometry'] = s

    # write output
    gdf_poly.to_file(outfolder + 'Ys_polygons_' + version + '.shp')
    gdf_fl.to_file(outfolder + 'flowlines_' + version + '.shp')


# #####################################################################################################################
def flowlines_plot(df_fl, dem, x_coords, y_coords, CRS, outfolder, version):  # create map with the flowlines

    # text formatting
    if plt.rcParams["text.usetex"]:
        fmt = r'%.0f'
    else:
        fmt = '%.0f'

    # calculate spacing of contour lines
    rel = np.isfinite(dem)
    interval = (np.floor((np.max(dem[rel]) - np.min(dem[rel])) / 800) + 1) * 100
    if interval > 200:
        interval = 200
    levels = np.arange(np.ceil(np.min(dem[rel]) / interval) * interval,
                       np.ceil(np.max(dem[rel]) / interval) * interval, interval)

    # calculate the bounding box
    bbox = ((np.min(x_coords), np.max(x_coords), np.min(y_coords), np.max(y_coords)))

    # set the range of velocity values
    vmax = df_fl['v'].max()
    if vmax > 400:
        vmax = 400
    norm = plt.Normalize(df_fl['v'].min(), vmax)

    # create the map
    fig, ax = plt.subplots(figsize=(18, 24))   #
    ax.set_xlabel('easting (' + CRS + ')')
    ax.set_ylabel('northing (' + CRS + ')')
    plt.title('Flowlines Test')
    ax.set_xlim(bbox[0], bbox[1])
    ax.set_ylim(bbox[2], bbox[3])

    flid = np.unique(df_fl['ID'])

    for f in flid:
        flt = df_fl.loc[df_fl['ID'] == f]
        points = np.array([flt['X'], flt['Y']]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        # create the line collection
        lc = LineCollection(segments, cmap='rainbow', norm=norm)
        # Set the values used for color mapping
        lc.set_array(flt['v'].astype(float))
        lc.set_linewidth(1.2)
        line = ax.add_collection(lc)

    cs1 = ax.contour(dem[:, :], levels=levels, extent=bbox, origin='upper',
                     zorder=2, colors='black', alpha=0.6)
    ax.clabel(cs1, cs1.levels, inline=True, fmt=fmt, fontsize=8)

    plt.gca().set_aspect('equal', adjustable='box')

    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    fig.colorbar(line, cax=cax, label='velocity (m yr$^{-1}$)')

    plt.tight_layout()
    plt.savefig(outfolder + 'flowlines_' + version + '.png')
    plt.close()











