import csv
import tqdm
import requests
import threading
import queue


##########################################################
# curie_match
# Input: drug
# Output: string chembl id or 'NOT FOUND' if not found
#--------------------------------------------------------
# Description: Drug chembl ids are returned using MyChem.info
# and must fit criteria for chembl_id.

def curie_match(drug):
    url = 'http://mychem.info/v1/query?q=chembl.pref_name:{}&fields=chembl.molecule_chembl_id'.format(drug)
    r = requests.get(url)
    r_json = r.json()
    hits = r_json['hits']
    if len(hits) > 0:
        curie = 'ChEMBL:' + hits[0]['chembl']['molecule_chembl_id']
    else:
        curie = 'NOT FOUND'
    return [drug,curie]


if __name__ == "__main__":

    drug_names_file = '/home/public/data/ncats/data_drop_03-13-2020/drug_names.csv'

    # reading initial patient data
    with open(drug_names_file, 'r') as csv_file:
        reader = csv.reader(csv_file)
        drug_names = next(reader)
        drugs = [drug for drug in drug_names]

    curie_matches = []
    for drug in drugs:
        curie_matches.append(curie_match(drug))

    f = open('drug_curie_map.csv', 'w')
    f.write('Drug,Curie\n')

    for d in curie_matches:
        row = str(d[0])+','+str(d[1])+'\n'
        f.write(row)
    f.close()
