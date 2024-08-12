# flowlines
Determine seedpoints on an ice sheet or glaciers, calculate flowlines based on velocity fields, convert flowlines 
into polygons that can then be used to determine slush limits.

1. `Crop_gdalwarp.py` can be used to crop a grid and bring it to the same grid as the MODIS data that will later be used. This piece of code is not only linked to the flowline calculation.
2. `flowline_seepoints` is used to compute seedpoints at a defined distance along a line polygon.
3. `flowline.py` calculates the flowlines and the related flowline-polygons, based on velocity fields as well as the seepoints.

`arrays_filter_GLvelocity.py` is used beforehand, to compute average flow fields in the case that multiple data
sets for different points in time exist. 

## Output
The final output are the flowlines and the flowline polygons which are provided as shapefiles (.shp)

## Remarks
If there are no specially high requirements to the accuracy of the flowines, then also relatively old flow fields
work fine. So far, only data from publications from the years 2015 and 2016 have been used.