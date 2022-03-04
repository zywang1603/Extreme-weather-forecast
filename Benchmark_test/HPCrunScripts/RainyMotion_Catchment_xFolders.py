# -*- coding: utf-8 -*-
"""
Created on Fri Jan 4 10:05:58 2019

@author: imhof002

RainyMotion script to run all four modules and to use the sliced outputs per 
catchment.
In addition, this script manages to loop through a main directory and opens
all folders in order to make nowcasts per event for 'n' selected events.

For simplicity, normally only the initial data has to be changed.
"""

###########
# Import package
###########

# import rainymotion library
from rainymotion import models, metrics, utils, Catchment_slice_RM

# import accompanying libraries
from collections import OrderedDict
import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
from os import listdir
from os.path import isfile, join
import time

###############################################################################
############
# Initial data - Basically only this has to be changed
############

# Set the directory with the input radar files (Note that this is the top directory,
# which contains multiple directories with all the cases).
Dir_in = "c:\\Users\\imhof_rn\\Documents\\Nowcasting\\Input\\Processor_1"

# Set number of rows and number of cols in the full radar image.
num_rows = 765
num_cols = 700

# Catchment filenames and directories
catchments = True # Put on false when you don't want any slicing for catchments (i.e. you will use the full output)
# If catchments = 'False', uncomment the next two lines.
catchment_filenames = ["c:\\Users\\imhof_rn\\Documents\\GIS\\Catchments_for1H\\Hupsel.shp", "c:\\Users\\imhof_rn\\Documents\\GIS\\Catchments_for1H\\stroomgebied_Regge.shp", "c:\\Users\\imhof_rn\\Documents\\GIS\\Catchments_for1H\\GroteWaterleiding.shp", "c:\\Users\\imhof_rn\\Documents\\GIS\\Catchments_for1H\\Aa.shp", "c:\\Users\\imhof_rn\\Documents\\GIS\\Catchments_for1H\\Reusel.shp", "c:\\Users\\imhof_rn\\Documents\\GIS\\Catchments_for1H\\het_molentje.shp", "c:\\Users\\imhof_rn\\Documents\\GIS\\Catchments_for1H\\Luntersebeek.shp", "c:\\Users\\imhof_rn\\Documents\\GIS\\Catchments_for1H\\Dwarsdiep.shp", "c:\\Users\\imhof_rn\\Documents\\GIS\\Catchments_for1H\\AfwaterendgebiedBoezemsysteem.shp", "c:\\Users\\imhof_rn\\Documents\\GIS\\Catchments_for1H\\HHRijnland.shp", "c:\\Users\\imhof_rn\\Documents\\GIS\\Catchments_for1H\\Beemster.shp", "c:\\Users\\imhof_rn\\Documents\\GIS\\Catchments_for1H\\DeLinde.shp"] # Put here the locations of the shapefiles
catchment_names = ['Hupsel', 'Regge', 'GroteWaterleiding', 'Aa', 'Reusel', 'Molentje', 'Luntersebeek', 'Dwarsdiep', 'Delfland', 'Rijnland', 'Beemster', 'Linde'] # A list of catchment names.

# Set the output file directory. The names will later be added based on the catchment names. 
File_out_dir = "c:\\Users\\imhof_rn\\Documents\\Nowcasting\\Results\\1hours" 

# Set the number of lead times
leadtime = 72
###############################################################################

Start_loop = time.time()

##########
# We start with defining the functions itself
##########

# Model run setups

def persistence(data_instance, eval_instance, results_instance):
    
    for i in range(0,len(catchment_names)):
        results_instance[i].create_group("/Persistence/")

    for key in sorted(list(eval_instance.keys())):

        inputs = np.array([ data_instance[key] for key in eval_instance[key][0][-1:] ]) / 100.0

        model = models.Persistence()

        model.input_data = inputs

        nowcast = model.run()
        
        nowcast_slice, chunks_slice, maxshape_slice = Catchment_slice_RM.catchment_slice(nowcast, catchment_filenames, leadtime)
        
        inputs = None
        nowcast = None
        
        for i in range(0,len(catchment_names)):
            results_instance[i]["/Persistence/"].create_dataset(key,
                                                             data=nowcast_slice[i],
                                                             dtype="float16",
                                                             maxshape=maxshape_slice[i],
                                                             compression="gzip")
        nowcast_slice = None
        chunks_slice = None
        maxshape_slice = None


# create ground truth predictions
def ground_truth(data_instance, eval_instance, results_instance):

    for i in range(0,len(catchment_names)):
        results_instance[i].create_group("/GT/")

    for key in sorted(list(eval_instance.keys())):

        ground_truth = np.array([ data_instance[key] for key in eval_instance[key][1] ]) / 100.0

        GT_slice, chunks_slice, maxshape_slice = Catchment_slice_RM.catchment_slice(ground_truth, catchment_filenames, leadtime)
        
        ground_truth = None
        
        for i in range(0, len(catchment_names)):
            results_instance[i]["/GT/"].create_dataset(key,
                                                    data=GT_slice[i],
                                                    dtype="float16",
                                                    maxshape=maxshape_slice[i], compression="gzip")
            
        GT_slice = None
        chunks_slice = None
        maxshape_slice = None        

def optical_flow(data_instance, eval_instance, results_instance, model_name):

    if model_name == "Sparse":
        model = models.Sparse()

    elif model_name == "SparseSD":
        model = models.SparseSD()

    elif model_name == "Dense":
        model = models.Dense()

    elif model_name == "DenseRotation":
        model = models.DenseRotation()

    for i in range(0,len(catchment_names)):
        results_instance[i].create_group("/{}/".format(model_name))

    for key in sorted(list(eval_instance.keys())):

        inputs = np.array([ data_instance[key] for key in eval_instance[key][0] ]) / 100.0

        model.input_data = inputs

        nowcast = model.run()
        
        nowcast_slice, chunks_slice, maxshape_slice = Catchment_slice_RM.catchment_slice(nowcast, catchment_filenames, leadtime)

        inputs = None
        nowcast = None
        
        for i in range(0, len(catchment_names)):
            results_instance[i]["/{}/".format(model_name)].create_dataset(key,
                                                             data=nowcast_slice[i],
                                                             dtype="float16",
                                                             maxshape=maxshape_slice[i],
                                                             compression="gzip")

        nowcast_slice = None
        chunks_slice = None
        maxshape_slice = None

##################
# Also the definition of some metrics
##################

# load a mask which maps RY product coverage
mask = [ [1] * num_cols for _ in range(num_rows)]
mask = np.array([mask for i in range(12)])

# Verification block
def calculate_CSI(obs, sim, thresholds=[0.125, 0.250, 0.500, 1.000]):

    result = {}

    for threshold in thresholds:
        result[str(threshold)] = [metrics.CSI(obs[i], sim[i], threshold=threshold) for i in range(obs.shape[0])]

    return result

def calculate_MAE(obs, sim):

    return [metrics.MAE(obs[i], sim[i]) for i in range(obs.shape[0])]

def calculate_metrics_dict(eval_instance, results_instance,
                           model_names=["Persistence", "Sparse", "SparseSD", "Dense", "DenseRotation"]):

    metrics_dict = OrderedDict()

    for model_name in model_names:

        metrics_dict[model_name] = OrderedDict()

        for key in sorted(list(eval_instance.keys())):

            metrics_dict[model_name][key] = {model_name: {"CSI": None, "MAE": None}}

            # observed ground-truth
            o = results_instance["GT"][key][()]

            # results of nowcasting
            s = results_instance[model_name][key][()]

            # convert values from depth (mm) to intensity (mm/h)
            o = utils.depth2intensity(o)
            s = utils.depth2intensity(s)

            # mask arrays
            o = np.ma.array(o, mask=mask)
            s = np.ma.array(s, mask=mask)

            metrics_dict[model_name][key][model_name]["CSI"] = calculate_CSI(o, s)
            metrics_dict[model_name][key][model_name]["MAE"] = calculate_MAE(o, s)

    return metrics_dict

##########
# The loop and actual work
##########

# We loop through every folder and make the nowcast per folder
for n in range(0, len(os.listdir(Dir_in))):
    Start = time.time()
    
    # Get the directory name of the case of this loop
    Case_name = os.listdir(Dir_in)[n]
    Case_dir = os.path.join(Dir_in, Case_name)
    # Make the output directory
    try:
        if not os.path.exists(os.path.join(File_out_dir, Case_name)):
            os.makedirs(os.path.join(File_out_dir, Case_name)) # Make this directory in the output folder
    except OSError:
        print('Error: Creating directory', Case_name)
        
    print('Loop starts for case', Case_name)
    print('State so far:', float(n) / float(len(os.listdir(Dir_in))), '% of the total procedure completed.')
    
    ############
    # Set output filenames
    ############
    File_out = []
    for i in range(0,len(catchment_names)):
        File_out.append(os.path.join(File_out_dir, Case_name, catchment_names[i]+'.h5'))
    
    ############
    # Import data
    ############
    
    # First make a list of all h5-files
    # Then import them and combine them in one h5-file. 
    combined = []
    dates = []
    
    print('Reading all files in directory')
    
    for filename in os.listdir(Case_dir):
        f=h5py.File(os.path.join(Case_dir, filename), 'r')
        dset=f['image1']['image_data']
        f_copy=np.copy(dset) #copy the content
    #    f_copy = np.where(f_copy == 65535, 0, f_copy)
        combined.append(f_copy)
        dates.append(filename.split("_")[4].split(".")[0])
        f.close()
    
    print(len(dates), 'files in this directory.')
    
    # create dictionary with timestep indexes
    # dictionary structure: {"t": [ [t-24, t-23,..., t-1], [t+1,...,t+12] ]}
    #eval_idx = np.load("M:\My Documents\RainyMotion\data\eval_dict.npy").item()
    eval_idx = {}
    for i in range(23, (len(dates) - 71)):
        date = dates[i]
        old_date = []
        new_date = []
        for old_num in range(0,24):
            old_date.append(dates[i-(24-old_num)])
        for new_num in range(0,72):
            new_date.append(dates[i+new_num])
        eval_idx[date] = [old_date,new_date]
        
        old_date = None
        new_date = None
    
    data = {}
    for j in range(0,len(dates)):
        date = dates[j]
        data[date] = combined[j]
    
    ############
    # Storage of results
    ############
    
    # create placeholder (or load previously calculated) results
    results = []
    
    for i in range(0,len(catchment_names)):
        results.append(h5py.File(File_out[i]))
    
    #################
    # Run the model
    #################
    
    print('Start with the model - this will take some time')
    
    ground_truth(data, eval_idx, results)
    
    print("Ground truth is calculated")
    gt_time = time.time()
    print('Process took', (gt_time-Start)/3600.0, 'hours')
    
#    persistence(data, eval_idx, results)
#    
#    print("Persistence is calculated")
#    per_time = time.time()
#    print('Process took', (per_time - gt_time)/3600.0, 'hours')
    
    optical_flow(data, eval_idx, results, "Sparse")
    
    print("Sparse model is calculated")
    sparse_time = time.time()
    print('Process took', (sparse_time - gt_time)/3600.0, 'hours')
    
#    optical_flow(data, eval_idx, results, "SparseSD")
#    
#    print("Sparse SD model is calculated")
#    
#    sparseSD_time = time.time()
#    print('Process took', (sparseSD_time - sparse_time)/3600.0, 'hours')
#    
#    optical_flow(data, eval_idx, results, "Dense")
#    
#    print("Dense model is calculated")
#    dense_time = time.time()
#    print('Process took', (dense_time - sparseSD_time)/3600.0, 'hours')
    
    optical_flow(data, eval_idx, results, "DenseRotation")
    
    print("Dense Rotation model is calculated")
    denseRot_time = time.time()
    print('Process took', (denseRot_time - sparse_time)/3600.0, 'hours')
    
    for i in range(0, len(results)):
        results[i].close()
    
    End = time.time()
    print('This case took in total', (End - Start)/3600.0, 'hours')

    # Finally and for certainty, delete all variables
    Case_name = None
    Case_dir = None
    File_out = None
    combined = None
    dates = None
    eval_idx = None
    data = None
    results = None
    
# And we're done!
End_loop = time.time()
print('Total process took', (End_loop - Start_loop)/3600.0, 'hours')

print("Done!")


