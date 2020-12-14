import os
import sys
import subprocess

sys.path.append('/home/cyakaboski/src/python/projects/bkb-pathway-provider/core')

from query import Query

DRIVER_PATH = '/home/cyakaboski/src/python/projects/bkb-pathway-provider/core'

#-- Setup Query 
for i in range(2):
    query0 = Query(name='query{}'.format(i),
                   evidence=dict(),
                   targets=list(),
                   meta_evidence=[('Age_of_Diagnosis', '<=', 9706)],
                   meta_targets=[('Survival_Time', '>=', 943)])

    #-- Save the query.
    pickle_file, json_file = query0.save(os.getcwd())

    #-- Build system command
    command = ['python3', os.path.join(DRIVER_PATH, 'driver.py'),
               '--config_file', os.path.join(DRIVER_PATH, 'driver.config'),
               '--headless',
               '--query_file', pickle_file,
               '--save_dir', os.getcwd()]
    subprocess.run(command)

