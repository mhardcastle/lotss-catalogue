# -------------------------------------------------------------
# Panoptes Marking Export Script
#
# This is a quick script to export the contents of a LOFAR Galaxy Zoo workflow export.
# As input it requires a workflow export file. Note that care needs to be taken about
# workflow version numbers as the Zooniverse export includes everything. I have found it easiest
# to manually remove earlier workflow version entires from the csv before running this,
# rather than filtering using this script.
##
# The output from the scipt is four csv files: a list of association click data, a list of optical ID clicks,
# a list of subject information, and a list of problem clicks.
#
# Original version by J. Croston 17/11/17, modified M. Hardcastle 2020-2022
#
# It follows various examples from Zooniverse, including the Andromeda project example from
# https://github.com/zooniverse/Data-digging
# which was written by: Cliff Johnson (lcj@ucsd.edu), based on scripts by Brooke Simmons 
# -------------------------------------------------------------

from __future__ import print_function
import sys
from collections import defaultdict
from astropy.table import Table
import numpy as np

def flatten_subject_info(subject_data):

  result = subject_data.values()[0]
  result.update({'id': subject_data.keys()[0]})
  return result

try:
    classfile_in = sys.argv[1]

except:
    print("\nUsage: "+sys.argv[0]+" classifications_infile")
    sys.exit(0)

    
markfile_out_assoc = 'lofar-associations.csv'
markfile_out_ids = 'lofar-opticalids.csv'
output_probs = 'lofar-classification-problems.csv'
nsubjfile='lofar-allclassifications.csv'

import pandas as pd
import json

t=Table.read('source_lr.fits')

sources=set(t['Source_Name'])

# Read in classification CSV and expand JSON fields
classifications = pd.read_csv(classfile_in)

classifications['metadata_json'] = [json.loads(q) for q in classifications.metadata]
classifications['annotations_json'] = [json.loads(q) for q in classifications.annotations]
classifications['subject_data_json'] = [json.loads(q) for q in classifications.subject_data]

classifications['subject_data']=classifications["subject_data"].map(json.loads)

classifications["subject_data"]=classifications["subject_data"].map(flatten_subject_info)

# First deal with blends, artefacts, too-zoomed-in and host-galaxy-broken cases

classifications['n_problems']=[ len(q[2]['value']) for q in classifications.annotations_json ]

# Calculate number of markings per classification
# Note: index of annotations_json ("q" here) corresponds to task number (i.e., 0)
classifications['n_markings'] = [ len(q[0]['value']) for q in classifications.annotations_json ]

# Select only classifications from most recent workflow version
# iclass = classifications[classifications.workflow_version >23.62]

iclass = classifications

# Output markings from classifications in iclass to new list of dictionaries (prep for pandas dataframe)
# Applicable for workflows with marking task as first task, and outputs data for circular markers (x,y,r)
clist=[]
c2list=[]
plist=[]


# Note: index of annotations_json ("q" here) corresponds to task number (i.e., 0)
classifications['na_markings'] = [ len(q[0]['value']) for q in classifications.annotations_json ]
classifications['ni_markings'] = [ len(q[1]['value']) for q in classifications.annotations_json ]


iclass=classifications

# Output markings from classifications in iclass to new list of dictionaries (prep for pandas dataframe)
# Applicable for workflows with marking task as first task, and outputs data for circular markers (x,y,r)


# Now deal with all the complicated classifications
nlist=[]

subjdict=defaultdict(list)

outfile=open('rescue.txt','w')
rescued={}

for index,c in iclass.iterrows():

    subj=c.subject_data['source_name']
    user=c.user_name
    ra=float(c.subject_data['ra'])
    dec=float(c.subject_data['dec'])
    size=float(c.subject_data['#size'])
    if subj not in sources:
      if subj not in rescued:
        print('Attempting to rescue subject',subj,'which is missing from table')
        print('RA is',ra,'and DEC is',dec)
        dra=t['RA']-ra
        ddec=t['DEC']-dec
        dist=3600*np.sqrt((np.cos(ra*np.pi/180)*dra)**2.0+ddec**2.0)
        index=np.argmin(dist)
        oldsubj=subj
        subj=t[index]['Source_Name']
        print('Rescued, new source name is',subj,'distance was',dist[index])
        outfile.write('%s %s %f\n' % (subj,oldsubj,dist[index]))
        rescued[oldsubj]=subj
      else:
        subj=rescued[subj]
        #else:
      #  raise RuntimeError('Failed to find matching source, nearest counterpart is at %f arcsec',np.sqrt(np.min(dist))*3600.0)

    lusers=subjdict.get(subj,[])
    if user in lusers:
        continue
    else:
        nlist.append({'classification_id':c.classification_id, 'user_name':c.user_name,'subject_ids':c.subject_ids,'source_name':subj,'ra':ra, 'dec':dec,'size':size})
        lusers.append(user)
        subjdict[subj]=lusers

        if c['ni_markings'] > 0:
            for q in c.annotations_json[1]['value']:
                clist.append({'classification_id':c.classification_id, 'user_name':user, 'user_id':c.user_id,'created_at':c.created_at, 'subject_ids':c.subject_ids, 'source_name':subj,'x':q['x'], 'y':q['y'], 'frame':q['frame']})
            
        if c['na_markings'] > 0:
            for q in c.annotations_json[0]['value']:
                c2list.append({'classification_id':c.classification_id, 'user_name':c.user_name, 'user_id':c.user_id,'created_at':c.created_at, 'subject_ids':c.subject_ids, 'source_name':subj,'x':q['x'], 'y':q['y'], 'frame':q['frame']})
 
    
        for q in c.annotations_json[2]['value']:
            plist.append({'classification_id':c.classification_id, 'user_name':c.user_name,'subject_ids':c.subject_ids,'source_name':subj,'problem':q})

            
col_order=['classification_id','user_name','user_id','created_at','subject_ids','source_name','x','y','frame']
out=pd.DataFrame(clist)[col_order]
out.to_csv(markfile_out_ids,index_label='mark_id')
            
out=pd.DataFrame(c2list)[col_order]
out.to_csv(markfile_out_assoc,index_label='mark_id')
                   
pcol_order=['subject_ids','source_name','user_name','problem']
pout=pd.DataFrame(plist)[pcol_order]
pout.to_csv(output_probs,index_label='mark_id')

ncol_order=['classification_id','user_name','subject_ids','source_name','ra','dec','size']
nout=pd.DataFrame(nlist)[ncol_order]
nout.to_csv(nsubjfile,index_label='new_id')
