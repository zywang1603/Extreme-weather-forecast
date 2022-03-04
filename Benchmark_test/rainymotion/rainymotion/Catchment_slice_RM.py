# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 10:42:13 2019

@author: imhof002

Script to find pixels corresponding to catchment extent in complete radar image
of KNMI.
"""
from osgeo import gdal
from osgeo import gdal_array
from osgeo import ogr, osr

import subprocess
import shapefile
import numpy as np
import os

def catchment_slice(X, catchment_filenames, leadtime):
    """Slice the output array in x- and y-direction in order to only save the 
    extent of one or multiple catchments. Currently, this tool only works for
    the KNMI radar datasets.
    
    Parameters
    -------------
    X : 3D-array, which is the output of the steps.py ensemble nowcast per time
        step.The array consists of:
            1. n ensemble members
            2. n rows
            3. n cols
    
    catchment_filenames: A list with the location (directory + filename) of 
        catchment shapefiles. The shapefiles will be used and rasterized in a
        loop.
        
    leadtime: The number of lead times in the nowcast.
        
    Returns
    -------------
    array_out: 4D-array, though similar to the input, but with sliced rows and cols. 
        Hence, this array consists of:
            1. n catchments
            2. n ensemble members
            3. n rows (sliced)
            4. n cols (sliced)
            
    chunks: The chunk sizes for saving, this will be the same as maxshape.
    
    maxshape: The number of cols, rows and leadtimes per sliced catchment.
    
    """
    

    ##########
    # Initial settings
    ##########
    
    # Set the info from the NL_RAD data
    xmin_nlrad = 0
    xmax_nlrad = 700
    ymin_nlrad = -4415
    ymax_nlrad = -3650
    
    ########
    # Reproject the shapefile
    ########
    
    for i in range(0, len(catchment_filenames)):
        #########################################################
        # Time to read the shapefile
        #########################################################
        filename = catchment_filenames[i]
    
        #########################################################
        # Get the projection of the shapefiles
        ##########################################################
    
        driver = ogr.GetDriverByName('ESRI Shapefile')
        dataset = driver.Open(filename)
    
        # from Layer
        layer = dataset.GetLayer()
        spatialRef = layer.GetSpatialRef()
        # from Geometry
        feature = layer.GetNextFeature()
        geom = feature.GetGeometryRef()
        spatialRef = geom.GetSpatialReference()
    	
    	
        #########################################################
        # Reproject the geometry of the shapefile
        #########################################################
    
        # set file names
        infile = filename
        dirname = os.path.dirname(filename)
        basename = os.path.basename(filename)
        basenametxt = os.path.splitext(basename)[0]
        outfile = os.path.join(dirname, basenametxt+'_Reprojected.shp')			
    	
        # The projection will be based on the projection of the given data set	
    	
        ShapefileProj = '+proj=stere +lat_0=90 +lon_0=0.0 +lat_ts=60.0 +a=6378.137 +b=6356.752 +x_0=0 +y_0=0'
    	
        subprocess.call(['ogr2ogr', '-t_srs', ShapefileProj, outfile, infile])
    
    ###########
    # Loop to rasterize the shapefile and then slice the radar array with the 
    # extent of the rasterized shapefile.
    ###########
    
    # This is going to be the output
    array_out = []
    chunks = []
    maxshape = []
    
    for i in range(0, len(catchment_filenames)):
        #########################################################
        # Time to read the shapefile
        #########################################################
        filename = catchment_filenames[i]
    
        # set file names in order to obtain the reprojected shapefile, which 
        # was made with the catchment_medata functionality.
        dirname = os.path.dirname(filename)
        basename = os.path.basename(filename)
        basenametxt = os.path.splitext(basename)[0]
        outfile = os.path.join(dirname, basenametxt+'_Reprojected.shp')			
    	
        #########################################################
        # Define the new, reprojected shapefile as sfnew
        # We are going to rasterize this shapefile, since the resulting radar raster (after
        # rasterizing the hdf5-dataset) should be exactly the same as this raster of this shapefile
        #########################################################
    	
#        sfnew = shapefile.Reader(outfile)
    
        # 1. Define pixel_size and NoData value of new raster
        NoData_value = -9999
    
        # 2. Open Shapefile
        driver = ogr.GetDriverByName('ESRI Shapefile')
        sf = driver.Open(outfile)
        source_layer = sf.GetLayer()
        # Obtain xmin, xmax, ymin and ymax as floats
        x_min, x_max, y_min, y_max = source_layer.GetExtent()
        
        # 3. Round the values for xmin, xmax, ymin and ymax (make a rough rounding)
        # in order to not get a too small extent. Finally, make integer values of it.
        xmin = int(round(x_min, 0) - 1)
        xmax = int(round(x_max, 0) + 1)
        ymin = int(round(y_min, 0) - 1)
        ymax = int(round(y_max, 0) + 1)
        
        x_res = 1
        y_res = 1
        
        cols = int( (xmax - xmin) / x_res )
        rows = int( (ymax - ymin) / y_res )
        
        # 4. Slice the array with the known coordinates - which equal the col and row number!
        xslice_min = xmin
        xslice_max = xmax
        yslice_max = (ymin - ymax_nlrad) * -1 # + 3650 in order to compensate for the 3650 
        yslice_min = (ymax - ymax_nlrad) * -1 # ymin and ymax are flipped, because the indexation is flipped as well
        
        array_new = X[:,yslice_min:yslice_max,xslice_min:xslice_max]
        array_out.append(array_new)
        
        chunks.append([(leadtime,rows,cols)])
        maxshape.append([leadtime,rows,cols])
        
        array_new = None
        xmin = None
        xmax = None
        ymin = None
        ymax = None
        xslice_min = None
        xslice_max = None
        yslice_min = None
        yslice_max = None

    return array_out, chunks, maxshape
    
