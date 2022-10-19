import pandas as pd
df = pd.read_csv('catchment_extreme.csv', sep = ' ', header=None, names=['time', 'prep'])
#df = df.sort_values(by='prep', ascending=False)
#print(df['time'])
record = []
for t in enumerate(df['time']):
    #print(t[1])
    if t[1]-5 not in record and t[1]-10 not in record and str(t[1])[-2:]!="55":
        record.append((t[1]))
for i in range(len(record)):
    record[i] = str(record[i])

print(record)
print(len(record))
record = []
df = pd.read_csv('catchment_extreme_2.csv', sep = ' ', header=None, names=['time', 'prep'])
#df = df.sort_values(by='prep', ascending=False)
#print(df['time'])
for t in enumerate(df['time']):
    #print(t[1])
    if t[1]-5 not in record and t[1]-10 not in record and t not in record and str(t[1])[-2:]!="55":
        record.append((t[1]))
for i in range(len(record)):
    record[i] = str(record[i])

print(record)
print(len(record))

df = pd.read_csv('catchment_extreme_3.csv', sep = ' ', header=None, names=['time', 'prep'])
#df = df.sort_values(by='prep', ascending=False)
#print(df['time'])
for t in enumerate(df['time']):
    #print(t[1])
    if t[1]-5 not in record and t[1]-10 not in record and t not in record and str(t[1])[-2:]!="55":
        record.append((t[1]))
for i in range(len(record)):
    record[i] = str(record[i])

print(record)
print(len(record))