file_path = "exoenzyme.csv" #Altered for each functional group
results = {}

with open(file_path,'r') as data:
    for line in data.readlines():
        if 'Host' in line:
            continue
        else:
            tmpdat = line.split(',')
            key = line[line.index(',') + 1: len(line)].rstrip()
            if key in results:
                results[key].append(tmpdat[0])
            else:
                results[key] = [tmpdat[0]]

print()
for items in results:
    print("{},{},{}".format(items,len(results[items]),','.join(results[items])))
