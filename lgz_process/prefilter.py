# Pre-filter the RGZ output to select particular subject sets

from __future__ import print_function

import sys
from astropy.table import Table
import pandas as pd
import json
import os

def flatten_subject_info(subject_data):

    result = subject_data.values()[0]
    result.update({'id': subject_data.keys()[0]})
    return result

t=Table.read('/data/lofar/DR2/catalogues/LoTSS_DR2_v100.srl.fits')

try:
    classfile_in = sys.argv[1]
    filter_field = sys.argv[2]

except:
    print("\nUsage: "+sys.argv[0]+" classifications_infile filter_field")
    sys.exit(0)

lines=open(classfile_in).readlines()

t=t[t['Mosaic_ID']==filter_field]

sourcenames=list(t['Source_Name'])
print('Total sources in mosaic',len(sourcenames))

# Read in classification CSV and expand JSON fields
classifications = pd.read_csv(classfile_in)

# make outfile...
os.mkdir(filter_field)
os.chdir(filter_field)
outfile=open(filter_field+'.csv','w')
outfile.write(lines[0])

classifications['metadata_json'] = [json.loads(q) for q in classifications.metadata]
classifications['annotations_json'] = [json.loads(q) for q in classifications.annotations]
classifications['subject_data_json'] = [json.loads(q) for q in classifications.subject_data]

classifications['subject_data']=classifications["subject_data"].map(json.loads)

classifications["subject_data"]=classifications["subject_data"].map(flatten_subject_info)

# now we have the subject data, filter on that

for index,c in classifications.iterrows():
    subj=c.subject_data['source_name']
    if subj in sourcenames:
        outfile.write(lines[index+1])

outfile.close()
