import requests
import os
import pickle
import json
import time

query_path = '/home/cyakaboski/src/python/projects/bkb-pathway-provider/tests/reasonerStdTests/sample_query.pk'

if os.path.exists(query_path):
    with open(query_path, 'rb') as f_:
        reasoner_std = pickle.load(f_)


payload = {'query': reasoner_std}
start_time = time.time()
r = requests.post('http://129.170.69.138/submitQuery/', json=payload)
print('response time = {}'.format(time.time() - start_time))
