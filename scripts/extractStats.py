import numpy as np
import pandas as pd
import os
import sys

arguments = len(sys.argv) - 1 
print ("The script is called with %i arguments" % (arguments))

if (arguments < 5): 
	exit

filename = sys.argv[1]
delim = sys.argv[2]
parallel = int(sys.argv[3])
iso512 = int(sys.argv[4])
sortingSeperate = int(sys.argv[5])

print ("Filename %s"%filename)
print ("File has delimieter %s" %(delim))
print ("File has parallel run header %i"%parallel)
print ("File is isotropic statistic %i"%iso512)
print ("File contains sorting time seperately %i"%sortingSeperate)



##### Strong duct 180
df = pd.read_csv(filename,  delimiter=delim) 
if 'Unnamed: 23' in df.keys():
	df = df.drop(['Unnamed: 23'], axis=1)

if (iso512!=1): 
	df = df.drop(["timeStamp", "inputMode", "dataPath", "rmsFilename", "timeStep", "hMin", "hMax", "runType"], axis=1)
else:
	df = df.drop(["timeStamp", "inputMode", "dataPath", "rmsFilename", "timeStep", "hMin", "hMax", "runType"], axis=1)

# Remove empty column(s) if there are any
#df = df.dropna(axis='columns')
#df = df.fillna(0)
#df = df.fillna(0)
# Insert column for TotalTime without Loading
if (parallel==1):
	if (sortingSeperate==1):
		df['totalTimeWithoutLoad'] = df['sortTimeGlobal'] + df['watershedTimeGlobal'] + df['communicationTimeGlobal']
		df['localTimeWithoutLoad'] = df['sortTimeLocalAvg'] + df['watershedTimeLocalAvg'] + df['communicationTimeLocalAvg']
	else:
		df['totalTimeWithoutLoad'] = df['watershedTimeGlobal'] + df['communicationTimeGlobal']
		df['localTimeWithoutLoad'] = df['watershedTimeLocalAvg'] + df['communicationTimeLocalAvg']

groupby = ['totalSizeX', 'totalSizeY', 'totalSizeZ', 
		'blockSizeX', 'blockSizeY', 'blockSizeZ', 
        'numNodesX', 'numNodesY', 'numNodesZ',
        'numNodesTotal', 'hSamples']
if (sortingSeperate==1):
	grouby.append('usedBucketing')

notgroupby = []
#print(df.keys())
for key in df.keys():
    if key not in groupby:
        notgroupby.append(key)

print("Avering over:", notgroupby)

# Filter out hSamples = 100
group = df.groupby(groupby);
mean_df = group.mean().reset_index()
mean_df.rename(columns=lambda x: x if x in groupby else x +"_mean", inplace=True)
std_df = group.std().reset_index()
std_df.rename(columns=lambda x: x if x in groupby else x +"_std", inplace=True)
# Information is not needed twice
std_df = std_df.drop(groupby, axis=1)
new_df = pd.concat([mean_df,std_df],axis=1)
fileout = filename[:-4] + "_stats.csv"
new_df.to_csv(fileout, sep=",", index=False,encoding='utf-8')
print("Wrote to file: %s"%fileout)