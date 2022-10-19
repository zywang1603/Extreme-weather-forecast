import pandas as pd
total_rainy_event = 288338
extreme_event = int(total_rainy_event * 0.01)
rainy_event = int(total_rainy_event * 0.05)
df = pd.read_csv('catchment_Delfland_MFBS.csv', sep = ' ', header=None, names=['time', 'prep'])
df = df.sort_values(by='prep', ascending=False)

df_extreme = df.iloc[0:extreme_event].reset_index(drop=True)
df_rainy = df.iloc[extreme_event:rainy_event].reset_index(drop=True)

df_extreme_trai = df_extreme[df_extreme['time']//(1e8) < 2015]
df_extreme_vali = df_extreme[(df_extreme['time']//(1e8) >= 2015) & (df_extreme['time']//(1e8) < 2018)]
df_extreme_test = df_extreme[df_extreme['time']//(1e8) >= 2019]

df_rainy_trai = df_rainy[df_rainy['time']//(1e8) < 2015]
df_rainy_vali = df_rainy[(df_rainy['time']//(1e8) >= 2015) & (df_rainy['time']//(1e8) < 2018)]
df_rainy_test = df_rainy[df_rainy['time']//(1e8) >= 2018]

df_trai = pd.concat([df_extreme_trai, df_rainy_trai]).reset_index(drop=True)
df_vali = pd.concat([df_extreme_vali, df_rainy_vali]).reset_index(drop=True)
df_test = pd.concat([df_extreme_test, df_rainy_test]).reset_index(drop=True)

#df_trai.to_csv('training_Regge08-14.csv', index = False, header = False)
#df_vali.to_csv('validation_Delfland15-17_20.csv', index = False, header = False)
#df_test.to_csv('testing_Delfland18-20_20.csv', index = False, header = False)

#df_extreme_trai.to_csv('training_Delfland08-14_ext.csv', index = False, header = False)
#df_extreme_vali.to_csv('validation_Delfland15-17_ext.csv', index = False, header = False)
#df_extreme_test.to_csv('testing_Delfland18-20_ext.csv', index = False, header = False)

print('Training Set', len(df_trai))
#print(df_trai)
print('Validation Set', len(df_vali))
#print(df_vali)
print('Testing Set', len(df_extreme_test))
#print(df_test)

print(df_extreme_test)


