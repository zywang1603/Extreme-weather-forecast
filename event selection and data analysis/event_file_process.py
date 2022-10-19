import csv
root_dir = "E:/Extreme weather Ruben code/Events/"
for season in ['Spring', 'Autumn', 'Summer', 'Winter']:
    for hour in ['1hours', '3hours', '6hours', '24hours']:
        for catchment in ['Delfland']:
            filename = root_dir + season + '/' + hour + '/' + catchment + '_' + hour + '.csv'
            print("Season:", season, "Duration:", hour)
            # initializing the titles and rows list
            fields = []
            rows = []
            # reading csv file
            with open(filename, 'r') as csvfile:
                # creating a csv reader object
                csvreader = csv.reader(csvfile)
                # extracting field names through first row
                fields = next(csvreader)
                # extracting each data row one by one
                for row in csvreader:
                    if not row==[]:
                        rows.append(row)

            # printing the field names
            #print('Field names are:' + ', '.join(field for field in fields))

            #print(rows)
            for row in rows:
                print(row[0], row[1])