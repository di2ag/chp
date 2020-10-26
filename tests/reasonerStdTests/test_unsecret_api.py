import requests
import os
import pickle
import json

query_path = '/home/ghyde/bkb-pathway-provider/tests/reasonerStdTests/sample_query1.pk'

if os.path.exists(query_path):
    with open(query_path, 'rb') as f_:
        reasoner_std = pickle.load(f_)

reasoner_std['reasoner_id'] = 'unsecret'
payload = {'query': reasoner_std}
r = requests.post('http://10.3.4.0/submitQuery/', json=payload)
print(r.text)
