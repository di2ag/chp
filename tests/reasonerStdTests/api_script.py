import requests
import os
import pickle
import json

query_path = '/home/cyakaboski/src/python/projects/bkb-pathway-provider/tests/reasonerStdTests/sample_query.pk'

if os.path.exists(query_path):
    with open(query_path, 'rb') as f_:
        reasoner_std = pickle.load(f_)


payload = {'query': reasoner_std}
r = requests.post('http://127.0.0.1:8000/submitQuery/', json=payload)
