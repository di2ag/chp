from chp_client import get_client
from chp_client.query import build_query
import itertools
import tqdm
import logging
import pickle
import random
import json

from chp.trapi_interface import TrapiInterface

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Set number of queries to build
NUM_QUERIES = 10

# Set seed
random.seed(111)

# Get client
#client = get_client()

# Get curies
logger.info('Getting curies.')
curies = TrapiInterface().get_curies()
#curies = client.curies()
logger.info('Got curies.')

# Build all simple single gene, single drug, breast cancer, survival queries.
queries = []
for _ in range(NUM_QUERIES):
    drug = [random.choice(list(curies["biolink:Drug"].keys()))]
    q = build_query(
        therapeutic=drug,
        num_gene_wildcards=1,
        one_hop=True,
    )
    print(json.dumps(q, indent=2))
    queries.append(q)

# Pickle the queries
with open('random_gene_one_hop_queries.pk', 'wb') as f_:
    pickle.dump(queries, f_)
