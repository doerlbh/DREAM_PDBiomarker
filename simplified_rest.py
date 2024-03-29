# simplified_rest.py

# Rest feature extraction for training data: demographics: syn10146552, table: syn10146553
import sys
sys.path.insert(0,"/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages")

import synapseclient
# use syn = synapseclient.login() if you've already set up your config file
# syn = synapseclient.login(email="doerlbh@gmail.com", password="12345678", rememberMe=True)
syn = synapseclient.login()

import pandas as pd
import json
import numpy as np

# read in the healthCodes of interest from demographics training table
demo_syntable = syn.tableQuery("SELECT * FROM syn10146552")
demo = demo_syntable.asDataFrame()
healthCodeList = ", ".join( repr(i) for i in demo["healthCode"]) 

# Query 'walking training table' for rest data recordIDs and healthCodes. 
INPUT_REST_ACTIVITY_TABLE_SYNID = "syn10146553"
actv_rest_syntable = syn.tableQuery(('SELECT "recordId", "healthCode", "deviceMotion_walking_rest.json.items" FROM {0} WHERE healthCode IN ({1}) AND "deviceMotion_walking_rest.json.items" is not null LIMIT 500').format(INPUT_REST_ACTIVITY_TABLE_SYNID, healthCodeList))
actv_rest = actv_rest_syntable.asDataFrame()
actv_rest['idx'] = actv_rest.index

######################
# Download JSON Files
######################
# bulk download rest JSON files containing sensor data
rest_json_files = syn.downloadTableColumns(actv_rest_syntable, "deviceMotion_walking_rest.json.items")
items = rest_json_files.items()

# create pandas dataframe of JSON filepaths and filehandleIDs
rest_json_files_temp = pd.DataFrame({"deviceMotion_walking_rest.json.items": [i[0] for i in items], "rest_json_file": [i[1] for i in items]})

# convert ints to strings for merging
actv_rest["deviceMotion_walking_rest.json.items"] = actv_rest["deviceMotion_walking_rest.json.items"].astype(str)

# merge IDs/healthCodes with JSON data
actv_rest_temp = pd.merge(actv_rest, rest_json_files_temp, on="deviceMotion_walking_rest.json.items")

# move the files after downloading
import shutil
for file_handle_id, path in items():
    shutil.move(path, "/Users/DoerLBH/Dropbox/git/DREAM_PDBiomarker/data")

####################
# Feature Extraction
####################
## PLACE YOUR FEATURE EXTRACTION CODE HERE ##
## THIS IS AN EXAMPLE CODE PLACEHOLDER ##

# temp example of feature extraction
# get average(mean) x-coordinate UserAcceleration for each file
x_accel = [] # initialize empty list for storing x-acceleration values
avg_items = [] # initialize empty list for storing mean values

# loop through each row in dataframe to read in json file
# grab the userAcceleration x-values and calculate the means
for row in actv_rest_temp["rest_json_file"]: 
	with open(row) as json_data:
		data = json.load(json_data)
		for item in data:
			x = item.get("userAcceleration").get("x")
			x_accel.append(x)
		avg = np.mean(x_accel)
		avg_items.append(avg)

# create new column in dataframe from list of averaged values
actv_rest_temp["meanXaccel"] = avg_items

# Remove unnecessary columns
actv_rest1 = actv_rest_temp.drop(["deviceMotion_walking_rest.json.items", "idx", "rest_json_file"], axis=1)

#### END FEATURE EXTRACTION SNIPPET EXAMPLE ##
