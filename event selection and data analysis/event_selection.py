import h5py
import numpy as np
import os
import random
from datetime import datetime, timedelta

from osgeo import gdal
from osgeo import gdal_array
from osgeo import ogr, osr

import subprocess
#import shapefile

prefix = "E://Extreme weather data_KNMI//RAD_NL25_RAC_MFBS_EM_5min//"
# use gauge adjusted data for event selection
# Method 1
def eventSelection_maxN(N, t):
    # Ruben's way for event
    # Input: whole dataset; Output: picked sample's starting date
    # N: total number of sample we want
    # N/2 sample of highest average, N/2 sample of highest peak intensity
    # Loop through each possible starting time point
    # Step 1: calculate average precipitation
    # Step 2: update the ordered dict
    result = []
    result1 = {}
    result2 = {}
    year = "2008"
    for m in range(7, 9):
        for d in range(1, 32):
            for h in range(24):
                for mi in [0]:  # every 30 mins
                    month = str(m) if m >= 10 else "0" + str(m)
                    date = str(d) if d >= 10 else "0" + str(d)
                    hour = str(h) if h >= 10 else "0" + str(h)
                    minute = str(mi) if mi >= 10 else "0" + str(mi)
                    start_time = year + month + date + hour + minute
                    print(start_time)
                    lead_times, obs_times = eventGeneration(start_time, obs_time=5, lead_time=t)
                    sample_times = obs_times + lead_times
                    if os.path.exists(prefix + year + "//" + month + "//" + \
                                   "RAD_NL25_RAC_MFBS_EM_5min_" + sample_times[0] + "_NL.h5") and \
                        os.path.exists(prefix + year + "//" + month + "//" + \
                                   "RAD_NL25_RAC_MFBS_EM_5min_" + sample_times[-1] + "_NL.h5"):

                        prep_average = 0
                        peak_average = 0
                        for time in sample_times:
                            filename = prefix + year + "//" + month + "//" + \
                                       "RAD_NL25_RAC_MFBS_EM_5min_" + time + "_NL.h5"
                            f = h5py.File(filename)['image1']['image_data']
                            f = np.array(f)
                            # f = np.where(f == 65535, -1, f)
                            # f = np.ma.masked_where(f <= 0, f)
                            f[f == 65535] = 0
                            f = f[100:500, 100:500] * 12 / 100
                            [rows, cols] = f.shape
                            prep_average += np.sum(f)/(rows*cols)
                            peak_average += np.max(f)
                        prep_average = prep_average/len(sample_times)
                        peak_average = peak_average/len(sample_times)

                        if len(result1)>=N:
                            result1_min = min(result1, key=result1.get)
                            if result1[result1_min] <= prep_average:
                                result1.pop(result1_min)
                                result1[start_time] = prep_average
                        else:
                            result1[start_time] = prep_average

                        if len(result2) >= N:
                            result2_min = min(result2, key=result2.get)
                            if result2[result2_min] <= peak_average:
                                result2.pop(result2_min)
                                result2[start_time] = peak_average
                        else:
                            result2[start_time] = peak_average

    result1 = sorted(result1.items(), key = lambda item:item[1], reverse=True)
    result2 = sorted(result2.items(), key = lambda item:item[1], reverse=True)
    i1 = 0
    i2 = 0
    while len(result) < N:
        while result1[i1][0] in result:
            i1+=1
        result.append(result1[i1][0])
        while result2[i2][0] in result:
            i2+=1
        result.append(result2[i2][0])
    return result

# Method 2
def eventSelection_threshold(X, Y, Z, A, B, t):
    # E.C. thesis's event selection
    # Input: whole dataset; Output: picked sample's starting date
    # Loop through each frame and label the frame X% >= A, Y%>=B
    # Loop through each possible starting time point
    # Label the sequence, >Z heavy rainfall
    result_light = []
    result_heavy = []
    for year in ["2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021"]:
        for m in range(1, 13):
            print("year/month:", year, "/", m)
            for d in range(1, 32):
                for h in range(24):
                    for mi in range(0, 60, 60):  # every 30 mins
                        month = str(m) if m >= 10 else "0" + str(m)
                        date = str(d) if d >= 10 else "0" + str(d)
                        hour = str(h) if h >= 10 else "0" + str(h)
                        minute = str(mi) if mi >= 10 else "0" + str(mi)
                        start_time = year + month + date + hour + minute
                        if os.path.exists(prefix + year + "//" + month + "//" + \
                                          "RAD_NL25_RAC_MFBS_EM_5min_" + start_time + "_NL.h5"):
                            lead_times, obs_times = eventGeneration(start_time, obs_time=5, lead_time=t)
                            sample_times = obs_times + lead_times
                            sample_label = { i : 0 for i in sample_times }
                            if os.path.exists(prefix + year + "//" + month + "//" + \
                                              "RAD_NL25_RAC_MFBS_EM_5min_" + sample_times[0] + "_NL.h5") and \
                                    os.path.exists(prefix + year + "//" + month + "//" + \
                                                   "RAD_NL25_RAC_MFBS_EM_5min_" + sample_times[-1] + "_NL.h5"):
                                heavy_frame = 0
                                no_rain = False
                                for time in sample_times:
                                    filename = prefix + year + "//" + month + "//" + \
                                                   "RAD_NL25_RAC_MFBS_EM_5min_" + time + "_NL.h5"
                                    f = h5py.File(filename)['image1']['image_data']
                                    f = np.array(f)
                                    f = np.where(f == 65535, -1, f)
                                    f = f[100:500, 100:500] * 12 / 100
                                    f_valid = np.ma.masked_where(f < 0,f)
                                    f_light = np.ma.masked_where(f <= A, f)
                                    f_heavy = np.ma.masked_where(f <= B, f)
                                    nAll = f_valid.count()
                                    nLight = f_light.count()
                                    nHeavy = f_heavy.count()
                                    if nLight/nAll >= X/100:
                                        sample_label[time] = 1
                                        if nHeavy/nAll >= Y/100:
                                            sample_label[time] = 2
                                            heavy_frame += 1
                                    else:
                                        no_rain = True
                                if not no_rain:
                                    if heavy_frame >= Z:
                                        print(start_time, " heavy")
                                        result_heavy.append(start_time)
                                    else:
                                        print(start_time, " light")
                                        result_light.append(start_time)
    return result_light, result_heavy

# Method 3
def eventSelection_importance(p_s, p_m, t, q_min):
    # Importance sampling scheme, Nature
    # Input: whole dataset; Output: picked sample's starting date
    # Loop through each possible starting time point
    # Step 1: generate t*h*w sample
    # Step 2: calculate the probability of sample: q_n = min{1, q_min + (m/c)*sum(x_sat)}, x_sat = 1-exp(-x/s)
    # Step 3: Store the data
    result = []
    year = "2008"
    for m in range(7, 8):
        for d in range(3, 32):
            for h in range(24):
                for mi in range(0, 60, 30):  # every 30 mins
                    month = str(m) if m >= 10 else "0" + str(m)
                    date = str(d) if d >= 10 else "0" + str(d)
                    hour = str(h) if h >= 10 else "0" + str(h)
                    minute = str(mi) if mi >= 10 else "0" + str(mi)
                    start_time = year + month + date + hour + minute
                    lead_times, obs_times = eventGeneration(start_time, obs_time=5, lead_time=t)
                    sample_times = obs_times + lead_times
                    x_sat_sum = 0
                    if os.path.exists(prefix + year + "//" + month + "//" + \
                                   "RAD_NL25_RAC_MFBS_EM_5min_" + sample_times[0] + "_NL.h5") and \
                        os.path.exists(prefix + year + "//" + month + "//" + \
                                   "RAD_NL25_RAC_MFBS_EM_5min_" + sample_times[-1] + "_NL.h5"):

                        for time in sample_times:
                            filename = prefix + year + "//" + month + "//" + \
                                       "RAD_NL25_RAC_MFBS_EM_5min_" + time + "_NL.h5"
                            f = h5py.File(filename)['image1']['image_data']
                            f = np.array(f)
                            # f = np.where(f == 65535, -1, f)
                            f = np.ma.masked_where(f == 65535, f)
                            c = np.ma.count(f)
                            f = f * 12 / 100
                            x_sat = np.sum(1-np.exp((-1/p_s)*f))
                            x_sat_sum += x_sat
                        q = min(1, q_min + (p_m*x_sat_sum/(c*(t+5))))
                        picked = random.choices([True,False], weights=[q,1-q])[0]
                        print(start_time, " probability:", q, "Picked:", picked)
                        if picked: result.append(start_time)
    return result

# utils
def eventGeneration(start_time, obs_time = 4 ,lead_time = 72):
    # Generate event based on starting time point, return a list: [[t-4,...,t-1,t], [t+1,...,t+72]]
    # Get the start year, month, day, hour, minute
    year = int(start_time[0:4])
    month = int(start_time[4:6])
    day = int(start_time[6:8])
    hour = int(start_time[8:10])
    minute = int(start_time[10:12])
    #print(datetime(year=year, month=month, day=day, hour=hour, minute=minute))
    times = [(datetime(year, month, day, hour, minute) + timedelta(minutes=5 * (x+1))) for x in range(lead_time)]
    lead = [dt.strftime('%Y%m%d%H%M') for dt in times]
    times = [(datetime(year, month, day, hour, minute) - timedelta(minutes=5 * x)) for x in range(obs_time)]
    obs = [dt.strftime('%Y%m%d%H%M') for dt in times]
    obs.reverse()
    return lead, obs


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
        outfile = os.path.join(dirname, basenametxt + '_Reprojected.shp')

        # The projection will be based on the projection of the given data set

        ShapefileProj = '+proj=stere +lat_0=90 +lon_0=0.0 +lat_ts=60.0 +a=6378137 +b=6356752 +x_0=0 +y_0=0'

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
        outfile = os.path.join(dirname, basenametxt + '_Reprojected.shp')

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
        xmin = int(round(x_min, 0) - 1)//1000
        xmax = int(round(x_max, 0) + 1)//1000
        ymin = int(round(y_min, 0) - 1)//1000
        ymax = int(round(y_max, 0) + 1)//1000
        #print("xmin", xmin, xmax, ymin, ymax)
        x_res = 1
        y_res = 1

        cols = int((xmax - xmin) / x_res)
        rows = int((ymax - ymin) / y_res)

        # 4. Slice the array with the known coordinates - which equal the col and row number!
        xslice_min = xmin
        xslice_max = xmax
        yslice_max = (ymin - ymax_nlrad) * -1  # + 3650 in order to compensate for the 3650
        yslice_min = (ymax - ymax_nlrad) * -1  # ymin and ymax are flipped, because the indexation is flipped as well

        print(filename, xslice_min, xslice_max, yslice_min, yslice_max)

        array_new = X[:, yslice_min:yslice_max, xslice_min:xslice_max]
        array_out.append(array_new)

        chunks.append([(leadtime, rows, cols)])
        maxshape.append([leadtime, rows, cols])

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


# Catchment_slice test
prefix_catchment = "E:/Extreme weather data_KNMI/catchments/"
catchment_names = ["Hupsel.shp", "Regge.shp", "GroteWaterleiding.shp", "Aa.shp", "Reusel.shp", "Luntersebeek.shp", "Dwarsdiep.shp", "HHRijnland.shp", "Beemster.shp", "DeLinde.shp", "Delfland.shp"]

catchment_file = []
for name in catchment_names:
    catchment_file.append(prefix_catchment + name)
#print(catchment_file)
year = "2008"
month = "08"
date = "03"
hour = "22"
minute = "00"
start_time = year + month + date + hour + minute
lead_times, obs_times = eventGeneration(start_time, obs_time=5, lead_time=24)
sample_times = obs_times + lead_times
path = prefix + year + "//" + month + "//" + "RAD_NL25_RAC_MFBS_EM_5min_" + sample_times[0] + "_NL.h5"
X = h5py.File(path)['image1']['image_data']
X = np.array(X)
X = np.expand_dims(X, axis=0)
array_out, chunks, maxshape = catchment_slice(X, catchment_file, 24)
print(array_out[1])

#print(chunks)
#print(maxshape)

# Function test
import sys
if len(sys.argv) > 1:
    if sys.argv[1] == "max":
        if len(sys.argv)>2: t = sys.argv[2]
        else: t = 72
        if len(sys.argv)>3: N = sys.argv[3]
        else: N = 10
        result = eventSelection_maxN(N, t)
    elif sys.argv[1] == "importance":
        result = eventSelection_importance(1, 0.1, 10, 2e-4)
    elif sys.argv[1] == "threshold":
        result_light, result_heavy = eventSelection_threshold(X=20, Y=5, Z=10, A=1, B=20, t=20)
        result = result_light + result_heavy
    print("Total number of sample: ", len(result))

#result = eventSelection_maxN(20, 24)
#result = eventSelection_importance(1, 0.1, 10, 2e-4)
#result = eventSelection_threshold(X=10, Y=0.1, Z=2, A=1, B=30, t=24)
#with open("selectedTime_threshold_2010_all.txt", "w") as f:
    #for r in result:
        #if r: f.write(r + "\n")
