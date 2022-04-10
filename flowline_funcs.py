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
def create_polygons(df_fl, CRS, buff, outfolder):  # create the polygons

    # convert to LineString
    flid = np.unique(df_fl['ID'])
    # prepare GeoDataFrame to store final polygons
    final_gdf = gpd.GeoDataFrame()
    final_gdf['geometry'] = None
    # Set the GeoDataFrame's coordinate system according to the DEM
    final_gdf.crs = CRS

    for f in flid:
        flt = df_fl.loc[df_fl['ID'] == f]
        geo = [xy for xy in zip(flt.X, flt.Y)]
        s = LineString(geo)
        final = s.buffer(buff, join_style=2)
        final_gdf.loc[f, 'geometry'] = final

    # write output
    final_gdf.to_file(outfolder + 'testpoly3.shp')


# #####################################################################################################################
def flowlines_plot(df_fl, dem, x_coords, y_coords, CRS, outfolder):  # create map with the flowlines

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
    norm = plt.Normalize(df_fl['v'].min(), df_fl['v'].max())

    # create the map
    fig, ax = plt.subplots(figsize=(9, 12))   #
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
        lc.set_array(flt['v'])
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
    plt.savefig(outfolder + 'flowlines_testplot.png')
    plt.close()











