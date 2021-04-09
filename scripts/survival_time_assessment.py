import pickle
import logging
from logging.handlers import QueueHandler, QueueListener
import time
import os
import sys
import tqdm
from multiprocessing import Pool, Value, Queue

from chp_data.patient_bkb_builder import PatientBkbBuilder
from chp_data.bkb_handler import BkbDataHandler
from pybkb.common.bayesianKnowledgeBase import bayesianKnowledgeBase as BKB

from chp.reasoner import ChpDynamicReasoner
from chp.query import Query

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

start_time = time.time()
# Load patient data
bkb_handler = BkbDataHandler()
# Read patient data
with open(bkb_handler.patient_data_pk_path, 'rb') as f_:
    patient_data = pickle.load(f_)

# Filter out patients with MX staging.
_patient_data = {}
for pat, feature_dict in patient_data.items():
    if feature_dict['path_m'] == 'MX':
        print('True')
        continue
    _patient_data[pat] = feature_dict

patient_data = _patient_data

# Setup Discretization Scheme
#CRITERIA = 'linspace'
CRITERIA = 'equal_frequency'
discretization_scheme = {
    "EFO:0000714": {
        "nbins": 5,
        "criteria": CRITERIA
        }
}

# Setup BKB builder
builder = PatientBkbBuilder(
    patient_data,
    bkb_handler,
    discretization_scheme=discretization_scheme,
    num_genes=None,
    num_drugs=None,
)


# Set build option: gene or drug
BUILD = 'gene'
RESET = True

if os.path.exists('collapsed_{}.bkb'.format(BUILD)) and not RESET:
    logger.info('Using saved collapsed bkb found in current directory,')
    collapsed_bkb = BKB().load('collapsed_{}.bkb'.format(BUILD), use_pickle=True)
else:
    if BUILD == 'gene':
        # Build gene bkb
        collapsed_bkb = builder.build_gene_bkb(
            interpolate=True,
            save_fragments=True,
            prelink_features = ["EFO:0000714"]
        )
    elif BUILD == 'drug':
        # Build drug bkb
        collapsed_bkb = builder.build_drug_bkb(
            interpolate=True,
            #save_fragments=True,
        )

    # Save the bkbs
    logger.info('Trying to save')
    collapsed_bkb.save('collapsed_{}.bkb'.format(BUILD), use_pickle=True)
    logger.info('Total Build Time: {}'.format(time.time() - start_time))

# Set up inference(s)
'''
evidence = {'_ENSEMBL:ENSG00000155657': 'True',
            'CHEMBL:CHEMBL83': 'True',
           }
'''
targets = ["EFO:0000714"]

# Run top N frequently mutated genes.
top_genes = ['_{}'.format(gene_curie) for _, gene_curie in sorted([(count, gene_curie) for gene_curie, count in builder.all_gene_counts.items()], reverse=True)[:1000]]
top_drugs = ['{}'.format(drug_curie) for _, drug_curie in sorted([(count, drug_curie) for drug_curie, count in builder.all_drug_counts.items()], reverse=True)[:10]]

with open('top_gene_drugs.pk', 'wb') as f_:
    pickle.dump((top_genes, top_drugs), f_)

# Set up Reasoner
dynamic_reasoner = ChpDynamicReasoner(bkb_handler=bkb_handler,
                                      patient_bkb_builder=builder,
                                      gene_prelinked_bkb_override=collapsed_bkb)

# Query run function to be used with multiprocessing
def worker_fn(gene, drug, print_bkb=False, print_inf=False):
    evidence = {
        gene: 'True',
        drug: 'True'
    }
    _query = Query(
        evidence=evidence,
        targets=targets
    )
    # Run query(s)
    _query =  dynamic_reasoner.run_query(_query)

    if print_bkb:
        _query.bkb.makeGraph()
    if print_inf:
        updates = _query.result.process_updates()
        for target_comp, state_prob in updates.items():
            for target_state in state_prob:
                print(target_comp, target_state)
                for idx in range(_query.result.number_of_inferences(target_comp, target_state)):
                    inf = _query.result.get_inference(target_comp, target_state, idx)
                    _query.result.graph_inference(inf)

    # Print/Collect results
    return {(gene, drug): _query.result.process_updates()}

def worker_init(q):
    # all records from worker processes go to qh and then into q
    qh = QueueHandler(q)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(qh)

def logger_init():
    q = Queue()
    # this is the handler for all log records
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(levelname)s: %(asctime)s - %(process)s - %(message)s"))

    # ql gets records from the queue and sends them to the handler
    ql = QueueListener(q, handler)
    ql.start()

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    # add the handler to the logger so records from this process are handled
    logger.addHandler(handler)

    return ql, q

def main():
    q_listener, q = logger_init()

    logging.info('Starting multiprocess run.')
    pool = Pool(20, worker_init, [q])
    result = pool.starmap(
        worker_fn,
        tuple([(gene, drug) for gene in top_genes for drug in top_drugs]),
    )
    pool.close()
    pool.join()
    q_listener.stop()
    return result

def bkb_debug():
    for gene in top_genes:
        for drug in top_drugs:
            print(gene)
            print(drug)
            res = worker_fn(gene, drug, print_bkb=False, print_inf=True)
            input('Continue?')

start_time = time.time()
result_list = main()
#bkb_debug()

# Process result
res = {}
for _res in result_list:
    res.update(_res)
'''
# Setup query(s)
res = {}
for i, gene in enumerate(top_genes):
    logger.info('On GENE: {}/{}'.format(i+1, len(top_genes)))
    for j, drug in enumerate(top_drugs):
        logger.info('On DRUG: {}/{}'.format(j+1,len(top_drugs))) 
        evidence = {
            gene: 'True',
            drug: 'True'
        }
        _query = Query(
            evidence=evidence,
            targets=targets
        )
        # Run query(s)
        _query =  dynamic_reasoner.run_query(_query)

        # Print/Collect results
        #_query.result.summary(include_contributions=False)
        res[(gene, drug)] = _query.result.process_updates()
        #print(_query.result.process_updates())
'''
logger.info('Total Multiprocessing time: {}'.format(time.time() - start_time))
with open('results-{}-{}.pk'.format(CRITERIA, time.time()), 'wb') as f_:
    pickle.dump(res, f_)
