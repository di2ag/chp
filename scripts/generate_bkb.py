import time
import logging
import pickle

from pybkb.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB
from pybkb.common.bayesianKnowledgeBase import BKB_I_node, BKB_component, BKB_S_node
from pybkb.python_base.fusion import fuse
from pybkb.python_base.fusion_collapse import collapse_sources
from chp_data.bkb_handler import BkbDataHandler

from chp.reasoner import Reasoner
from chp.query import Query
from chp.util import process_operator
from chp.patientBKFProcessor import PatientProcessor

#-- Setup logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

t1 = time.time()

bkb_handler = BkbDataHandler()
with open(bkb_handler.patient_data_pk_path, 'rb') as f_:
    patient_data = pickle.load(f_)
    patient_processor = PatientProcessor()
    patient_processor.loadFromPatientData(bkb_handler.patient_data_pk_path, patient_fraction=.2)

hashes = []
bkfs = []
#-- Process BKFs
bad_genes = patient_processor.processPatientBKF_v2(interpolation_model='bigram',
                                                        interpolation_selection='frequency_based',
                                                        gene_subset_top_k=2)
#-- Extract patient bkfs based on hashes
for patient, bkf in zip(patient_processor.requested_patients, patient_processor.bkfs):
    #bkf.makeGraph()
    #input()
    bkfs.append(bkf)
    hashes.append(str(patient.patientHash))
#-- Add interpolator
bkfs.append(patient_processor.interpolator)
#patient_processor.interpolator.makeGraph()
hashes.append('interpolator')
fused_bkb = fuse(bkfs,
             [1 for _ in range(len(bkfs))],
             hashes)
#fused_bkb.makeGraph()

fused_bkb.save('fused_bkb_sample.bkb', use_pickle=True)
logger.info('Make BKB Ok.')
logger.info('Elapsed Time: {} sec'.format(time.time() - t1))
'''
drugs = set()
for patient, data in patient_data.items():
    drugs.add(data['Drug_Name(s)'])

print(drugs)

print(patient_data[list(patient_data.keys())[0]].keys())
'''

col_bkb = collapse_sources(fused_bkb)
#col_bkb.makeGraph()

reasoner = Reasoner(fused_bkb=fused_bkb, patient_data=patient_data)
query = Query(
              meta_evidence= [('Drug_Name(s)', '==', 'HERCEPTIN')],
              meta_targets = [('Survival_Time', '>=', 1000)])

res = reasoner.analyze_query(query, target_strategy='explicit', interpolation='standard')
res.bkb.makeGraph()
print('done')
