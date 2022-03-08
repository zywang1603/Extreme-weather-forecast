import h5py
import numpy as np
import os
import random
from datetime import datetime, timedelta

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
    for m in range(7, 8):
        for d in range(3, 15):
            for h in range(24):
                for mi in range(0, 60, 30):  # every 30 mins
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
    print(result1)
    print(result2)
    while len(result) <= N:
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
    year = "2008"
    for m in range(7, 8):
        for d in range(3, 11):
            for h in range(24):
                for mi in range(0, 60, 30):  # every 30 mins
                    month = str(m) if m >= 10 else "0" + str(m)
                    date = str(d) if d >= 10 else "0" + str(d)
                    hour = str(h) if h >= 10 else "0" + str(h)
                    minute = str(mi) if mi >= 10 else "0" + str(mi)
                    start_time = year + month + date + hour + minute
                    print(start_time)
                    lead_times, obs_times = eventGeneration(start_time, obs_time=5, lead_time=t)
                    sample_times = obs_times + lead_times
                    sample_label = { i : 0 for i in sample_times }
                    if os.path.exists(prefix + year + "//" + month + "//" + \
                                      "RAD_NL25_RAC_MFBS_EM_5min_" + sample_times[0] + "_NL.h5") and \
                            os.path.exists(prefix + year + "//" + month + "//" + \
                                           "RAD_NL25_RAC_MFBS_EM_5min_" + sample_times[-1] + "_NL.h5"):
                        heavy_frame = 0
                        for time in sample_times:
                            filename = prefix + year + "//" + month + "//" + \
                                           "RAD_NL25_RAC_MFBS_EM_5min_" + time + "_NL.h5"
                            f = h5py.File(filename)['image1']['image_data']
                            f = np.array(f)
                            f = np.where(f == 65535, -1, f)
                            f = f[100:500, 100:500] * 12 / 100
                            f_light = np.ma.masked_where(f <= A, f)
                            f_heavy = np.ma.masked_where(f <= B, f)
                            [rows, cols] = f.shape
                            rc = rows*cols
                            nLight = f_light.count()
                            nHeavy = f_heavy.count()
                            if nLight/rc >= X/100:
                                sample_label[time] = 1
                                if nHeavy/rc >= Y/100:
                                    sample_label[time] = 2
                                    heavy_frame += 1
                            else:
                                break
                        if heavy_frame >= Z:
                            result_heavy.append(start_time)
                        else:
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
                            # f = np.ma.masked_where(f <= 0, f)
                            f[f == 65535] = 0
                            f = f[100:500, 100:500] * 12 / 100
                            [rows, cols] = f.shape
                            x_sat = np.sum(1-np.exp((-1/p_s)*f))
                            x_sat_sum += x_sat
                        q = min(1, q_min + (p_m*x_sat_sum/(rows*cols*(t+5))))
                        print(start_time, " probability:", q)
                        picked = random.choices([True,False], weights=[q,1-q])
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

# Function test
lead, obs = eventGeneration('200910080600')
#print(lead)
result = eventSelection_importance(1, 0.1, 10, 2e-4)
print(result)
#result_light, result_heavy = eventSelection_threshold(X = 20, Y = 5, Z = 10, A = 1, B = 20, t=20)
#print(result_heavy)
#print(result_light)
#result_max = eventSelection_maxN(20, 20)
#print(result_max)