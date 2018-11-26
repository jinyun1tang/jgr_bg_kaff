import csv
sand=[]
clay=[]
soil_type=[]
k1s=0
k2s=0
with open('data/franz_soil.csv', 'rb') as f:
    reader = csv.reader(f,delimiter=',')
    k=0
    for row in reader:
        print row[0],row[1],row[2]
        if k>=1:
            sand.append(float(row[1]))
            clay.append(float(row[2]))
        k=k+1
import numpy as np

pct_sand=np.array(sand)
pct_clay=np.array(clay)

kt=len(pct_sand)
k=0
while k < kt:

    #damm 1
    #su 1
    #supeca 1
    #dam 2
    #su 2
    #supeca 2
    k=k+1
