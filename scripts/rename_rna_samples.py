import sys
import re
import os

infile = sys.argv[1]

filename_path = str(infile)
pathparts = filename_path.split('/')
path = pathparts[0] + '/' + pathparts[1] + '/' + pathparts[2] + '/' 
filename = filename_path.split('/')[3]
nameparts = filename.split('_')
if nameparts[1] == "1" or nameparts[1] == "2":
	sampleID = nameparts[0] + '_' + nameparts[1]
	barcode = nameparts[2]
	lane = nameparts[3]
	readnum = nameparts[4]
	filetype = nameparts[5].split('.')[1] + '.' + nameparts[5].split('.')[2]
else:
	sampleID = nameparts[0]
	barcode = nameparts[1]
	lane = nameparts[2]
	readnum = nameparts[3]
	filetype = nameparts[4].split('.')[1] + '.' + nameparts[4].split('.')[2]
newname = sampleID + '_' + lane + '_' + readnum + '.' + filetype
newname_path = path + newname
os.rename(filename_path, newname_path)
