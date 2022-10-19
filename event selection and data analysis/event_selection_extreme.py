from datetime import datetime, timedelta
import os
import h5py
import numpy as np
import csv
# open the file in the write mode
file = open('catchment_Aa.csv', 'w')
# create the csv writer
writer = csv.writer(file, delimiter=' ')
# catchment coordinate
catchment={
'Hupsel': (473, 476, 421, 424),
'Regge': (450, 494, 371, 410),
'GroteWaterleiding': (456, 467, 409, 421),
'Aa': (386, 431, 467, 515),
'Reusel': (369, 381, 487, 517),
'Luntersebeek': (391, 406, 422, 432),
'Dwarsdiep': (427, 442, 290, 304),
'HHRijnland': (332, 343, 429, 444),
'Beemster': (340, 351, 370, 382),
'DeLinde': (414, 439, 322, 340),
'Delfland': (296, 322, 432, 457)}
def daterange(start, end, step=timedelta(1)):
    curr = start
    while curr < end:
        yield curr
        curr += step
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
# data analysis
catchmentName = 'Aa'
prefix = "E://Extreme weather data_KNMI//RAD_NL25_RAC_MFBS_EM_5min//"
lead_time = 36
time_interval = 5
start = datetime(2015, 7, 28, 10, 35)
end = datetime(2015, 7, 28, 13, 5)
dates = daterange(start, end, step = timedelta(minutes = time_interval))
for time in dates:
    event_start_time = time.strftime('%Y%m%d%H%M')
    event_lead_list, event_obs_list = eventGeneration(event_start_time, 5 ,lead_time)
    event_exist = True
    year_start = event_lead_list[0][0:4]
    month_start = event_lead_list[0][4:6]
    year_end = event_lead_list[-1][0:4]
    month_end = event_lead_list[-1][4:6]
    if os.path.exists(prefix + year_start + "//" + month_start + "//" + \
                      "RAD_NL25_RAP_5min_" + event_lead_list[0] + ".h5") and \
        os.path.exists(prefix + year_end + "//" + month_end + "//" + \
                        "RAD_NL25_RAP_5min_" + event_lead_list[-1] + ".h5"):
        event_average = 0
        for t in event_lead_list:
            year = t[0:4]
            month = t[4:6]
            filename = prefix + year + "//" + month + "//" + \
                       "RAD_NL25_RAP_5min_" + t + ".h5"
            f = h5py.File(filename)['image1']['image_data']
            f = np.array(f)
            xmin, xmax, ymin, ymax = catchment[catchmentName]
            f = np.where(f == 65535, 0, f)
            f = f[ymin:ymax, xmin:xmax] / 100
            if t == event_lead_list[0]:
                f_accum = f
            else:
                f_accum += f
        event_average = np.mean(f_accum)
        #writer.writerow([event_start_time, event_average])
        #writer.writerow([event_average])
        print(event_start_time, event_average)
file.close()


