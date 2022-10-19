import csv
'''
count = {'<0.01':0, '0.01 - 0.1': 0, '0.1 - 1': 0, '1 - 2': 0, '2 - 3': 0, '3 - 4': 0, '4 - 5':0, '5 - 6': 0, '6 - 7': 0, '7 - 8': 0, '8 - 9':0, '9 - 10':0, '> 10':0}
with open('catchment_Delfland_RAP1.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if row:
            v = float(row.copy()[0][13:])
            if v > 10:
                count['> 10'] += 1
            elif v > 9:
                count['9 - 10'] += 1
            elif v > 8:
                count['8 - 9'] += 1
            elif v > 7:
                count['7 - 8'] += 1
            elif v > 6:
                count['6 - 7'] += 1
            elif v > 5:
                count['5 - 6'] += 1
            elif v > 4:
                count['4 - 5'] += 1
            elif v > 3:
                count['3 - 4'] += 1
            elif v > 2:
                count['2 - 3'] += 1
            elif v > 1:
                count['1 - 2'] += 1
            elif v > 0.1:
                count['0.1 - 1'] += 1
            elif v > 0.01:
                count['0.01 - 0.1'] += 1
            else:
                count['<0.01'] += 1
with open('catchment_Delfland_RAP2.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if row:
            v = float(row.copy()[0][13:])
            if v > 10:
                count['> 10'] += 1
            elif v > 9:
                count['9 - 10'] += 1
            elif v > 8:
                count['8 - 9'] += 1
            elif v > 7:
                count['7 - 8'] += 1
            elif v > 6:
                count['6 - 7'] += 1
            elif v > 5:
                count['5 - 6'] += 1
            elif v > 4:
                count['4 - 5'] += 1
            elif v > 3:
                count['3 - 4'] += 1
            elif v > 2:
                count['2 - 3'] += 1
            elif v > 1:
                count['1 - 2'] += 1
            elif v > 0.1:
                count['0.1 - 1'] += 1
            elif v > 0.01:
                count['0.01 - 0.1'] += 1
            else:
                count['<0.01'] += 1
# Delfland RAP

count_RAP = {'<0.01': 953707, '0.01 - 0.1': 154676, '0.1 - 1': 191206, '1 - 2': 42297, '2 - 3': 13936, '3 - 4': 6067, '4 - 5': 2539, '5 - 6': 1289, '6 - 7': 619, '7 - 8': 403, '8 - 9': 297, '9 - 10': 116, '> 10': 46}
total_number = sum(count_RAP.values())
percent_Delf_RAP = count_RAP.copy()
for key in percent_Delf_RAP:
    value = percent_Delf_RAP[key]
    percent_Delf_RAP[key] = (value, round(value/total_number*100,3))
print("Delfland RAP:", percent_Delf_RAP)
'''

# Aa MFBS
'''
count = {'<0.01':0, '0.01 - 0.1': 0, '0.1 - 1': 0, '1 - 2': 0, '2 - 3': 0, '3 - 4': 0, '4 - 5':0, '5 - 6': 0, '6 - 7': 0, '7 - 8': 0, '8 - 9':0, '9 - 10':0, '> 10':0}
with open('catchment_Delfland_MFBS.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if row:
            v = float(row.copy()[0][13:])
            if v > 10:
                count['> 10'] += 1
            elif v > 9:
                count['9 - 10'] += 1
            elif v > 8:
                count['8 - 9'] += 1
            elif v > 7:
                count['7 - 8'] += 1
            elif v > 6:
                count['6 - 7'] += 1
            elif v > 5:
                count['5 - 6'] += 1
            elif v > 4:
                count['4 - 5'] += 1
            elif v > 3:
                count['3 - 4'] += 1
            elif v > 2:
                count['2 - 3'] += 1
            elif v > 1:
                count['1 - 2'] += 1
            elif v > 0.1:
                count['0.1 - 1'] += 1
            elif v > 0.01:
                count['0.01 - 0.1'] += 1
            else:
                count['<0.01'] += 1'''
count_aa = {'<0.01': 964119, '0.01 - 0.1': 138860, '0.1 - 1': 162975, '1 - 2': 47460, '2 - 3': 22830, '3 - 4': 12309, '4 - 5': 6839, '5 - 6': 4301, '6 - 7': 2464, '7 - 8': 1604, '8 - 9': 1120, '9 - 10': 758, '> 10': 1784}
total_number = sum(count_aa.values())
percent_Aa_MFBS = count_aa.copy()
for key in percent_Aa_MFBS:
    value = percent_Aa_MFBS[key]
    percent_Aa_MFBS[key] = (value, round(value/total_number*100,3))
count_aa_rain = count_aa.copy()
count_aa_rain.pop('<0.01')
count_aa_rain.pop('0.01 - 0.1')
total_number_rain_aa = sum(count_aa_rain.values())
for key in count_aa_rain:
    value = count_aa_rain[key]
    count_aa_rain[key] = (value, round(value/total_number_rain_aa*100,3))
print("Aa MFBS:", percent_Aa_MFBS)
print("Aa total rain event:", total_number_rain_aa)
print('Aa rain:', count_aa_rain)
print(' ')
'''
count = {'<0.01':0, '0.01 - 0.1': 0, '0.1 - 1': 0, '1 - 2': 0, '2 - 3': 0, '3 - 4': 0, '4 - 5':0, '5 - 6': 0, '6 - 7': 0, '7 - 8': 0, '8 - 9':0, '9 - 10':0, '> 10':0}
with open('catchment_Delfland_MFBS_fixed.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if row:
            v = float(row.copy()[0][13:])
            if v > 10:
                count['> 10'] += 1
            elif v > 9:
                count['9 - 10'] += 1
            elif v > 8:
                count['8 - 9'] += 1
            elif v > 7:
                count['7 - 8'] += 1
            elif v > 6:
                count['6 - 7'] += 1
            elif v > 5:
                count['5 - 6'] += 1
            elif v > 4:
                count['4 - 5'] += 1
            elif v > 3:
                count['3 - 4'] += 1
            elif v > 2:
                count['2 - 3'] += 1
            elif v > 1:
                count['1 - 2'] += 1
            elif v > 0.1:
                count['0.1 - 1'] += 1
            elif v > 0.01:
                count['0.01 - 0.1'] += 1
            else:
                count['<0.01'] += 1

print(count)'''

count_delf = {'<0.01': 931117, '0.01 - 0.1': 137779, '0.1 - 1': 176938, '1 - 2': 52492, '2 - 3': 27415, '3 - 4': 15358, '4 - 5': 9126, '5 - 6': 5447, '6 - 7': 3745, '7 - 8': 2630, '8 - 9': 1511, '9 - 10': 1078, '> 10': 2787}
total_number = sum(count_delf.values())
percent_delf_MFBS = count_delf.copy()
for key in percent_delf_MFBS:
    value = percent_delf_MFBS[key]
    percent_delf_MFBS[key] = (value, round(value/total_number*100,3))

count_delf_rain = count_delf.copy()
count_delf_rain.pop('<0.01')
count_delf_rain.pop('0.01 - 0.1')
total_number_rain_delf = sum(count_delf_rain.values())
for key in count_delf_rain:
    value = count_delf_rain[key]
    count_delf_rain[key] = (value, round(value/total_number_rain_delf*100,3))
print("Delf MFBS:", percent_delf_MFBS)
print("Delf total rain event:", total_number_rain_delf)
print('Delf rain', count_delf_rain)
print('')

count_dwar = {'<0.01': 963855, '0.01 - 0.1': 120230, '0.1 - 1': 171416, '1 - 2': 50811, '2 - 3': 25301, '3 - 4': 14333, '4 - 5': 7920, '5 - 6': 4976, '6 - 7': 2951, '7 - 8': 1968, '8 - 9': 1164, '9 - 10': 768, '> 10': 1730}
percent_dwar_MFBS = count_dwar.copy()
for key in percent_dwar_MFBS:
    value = percent_dwar_MFBS[key]
    percent_dwar_MFBS[key] = (value, round(value/total_number*100,3))
count_dwar_rain = count_dwar.copy()
count_dwar_rain.pop('<0.01')
count_dwar_rain.pop('0.01 - 0.1')
total_number_rain_dwar = sum(count_dwar_rain.values())
for key in count_dwar_rain:
    value = count_dwar_rain[key]
    count_dwar_rain[key] = (value, round(value/total_number_rain_dwar*100,3))
print("Dwar MFBS:", percent_dwar_MFBS)
print("Dwar total rain event:", total_number_rain_dwar)
print('Dwar rain', count_dwar_rain)
print('')

# Regge
'''
count = {'<0.01':0, '0.01 - 0.1': 0, '0.1 - 1': 0, '1 - 2': 0, '2 - 3': 0, '3 - 4': 0, '4 - 5':0, '5 - 6': 0, '6 - 7': 0, '7 - 8': 0, '8 - 9':0, '9 - 10':0, '> 10':0}
with open('catchment_Regge_MFBS.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        if row:
            v = float(row.copy()[0][13:])
            if v > 10:
                count['> 10'] += 1
            elif v > 9:
                count['9 - 10'] += 1
            elif v > 8:
                count['8 - 9'] += 1
            elif v > 7:
                count['7 - 8'] += 1
            elif v > 6:
                count['6 - 7'] += 1
            elif v > 5:
                count['5 - 6'] += 1
            elif v > 4:
                count['4 - 5'] += 1
            elif v > 3:
                count['3 - 4'] += 1
            elif v > 2:
                count['2 - 3'] += 1
            elif v > 1:
                count['1 - 2'] += 1
            elif v > 0.1:
                count['0.1 - 1'] += 1
            elif v > 0.01:
                count['0.01 - 0.1'] += 1
            else:
                count['<0.01'] += 1

print(count)'''
count_regge = {'<0.01': 955728, '0.01 - 0.1': 135300, '0.1 - 1': 171608, '1 - 2': 48820, '2 - 3': 24289, '3 - 4': 12778, '4 - 5': 7454, '5 - 6': 4339, '6 - 7': 2553, '7 - 8': 1509, '8 - 9': 1041, '9 - 10': 649, '> 10': 1355}
percent_regge_MFBS = count_regge.copy()
for key in percent_regge_MFBS:
    value = percent_regge_MFBS[key]
    percent_regge_MFBS[key] = (value, round(value/total_number*100,3))
count_regge_rain = count_regge.copy()
count_regge_rain.pop('<0.01')
count_regge_rain.pop('0.01 - 0.1')
total_number_rain_regge = sum(count_regge_rain.values())
for key in count_regge_rain:
    value = count_regge_rain[key]
    count_regge_rain[key] = (value, round(value/total_number_rain_regge*100,3))
print("Regge MFBS:", percent_regge_MFBS)
print("Regge total rain event:", total_number_rain_regge)
print('Regge rain', count_regge_rain)
print('')