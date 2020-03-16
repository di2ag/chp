import os
import sys
import subprocess
import json

sys.path.append('/home/cyakaboski/src/python/projects/bkb-pathway-provider/core')
#sys.path.append('/home/ncats/live/core')

from query import Query

SERVICE_PATH = '/home/cyakaboski/src/python/projects/bkb-pathway-provider/core'
#SERVICE_PATH = '/home/ncats/live/core'

def saveLikeNcatsClient(query, json_file=None):
    #-- Make dictionary
    query_dict = {'name': query.name,
                         'genetic_evidence': query.evidence,
                         'demographic_evidence': query.meta_evidence,
                         'genetic_targets': query.targets,
                         'demographic_targets': query.meta_targets}

    #-- If no file provided save to working directory.
    if json_file is None:
        json_file = os.path.join(os.getcwd(), '{}.json'.format(query.name))

    with open(json_file, 'w') as out_file:
        json.dump(query_dict, out_file)

    return json_file

#-- Setup Query 
query0 = Query(name='query1',
               evidence={"var_5'FLANK=": "True",
                           "var_3'UTR=": "True"},
               targets=list(),
               meta_evidence=list(),
               meta_targets=[('Survival_Time', '>=', 365)])

#-- Save the query.
json_file = saveLikeNcatsClient(query0)

#-- Build system command
command = ['python3', os.path.join(SERVICE_PATH, 'bkb-service.py'),
           '--f', json_file]
subprocess.run(command)


#-- Setup Query 
query0 = Query(name='query2',
               evidence={"var_5'FLANK=": "True",
                           "var_3'UTR=": "True"},
               targets=list(),
               meta_evidence=[('Age_of_Diagnosis', '>=', 20000)],
               meta_targets=[('Survival_Time', '>=', 365)])

#-- Save the query.
json_file = saveLikeNcatsClient(query0)

#-- Build system command
command = ['python3', os.path.join(SERVICE_PATH, 'bkb-service.py'),
           '--f', json_file]
subprocess.run(command)
