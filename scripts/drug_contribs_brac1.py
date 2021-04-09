import json
import glob
import logging
import os
import pickle
from collections import defaultdict

from chp.trapi_interface import TrapiInterface

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# OPTIONS
RUN = False
RESULTS_FILENAME = 'drug_contrib_results.pk'

def get_year_from_filename(filename):
    split1 = filename.split('.')[0]
    split2 = split1.split('_')[-1]
    return split2


def parse_results(results):
    parsed = defaultdict(list)
    for _file, res in results.items():
        year = get_year_from_filename(_file)
        kg = res["message"]["knowledge_graph"]
        _results = res["message"]["results"]
        for _res in _results[1:]:
            eb = _res["edge_bindings"]
            # Get KG edge
            kg_edge = kg["edges"][eb["e1"][0]["id"]]
            contrib = kg_edge["attributes"][0]["value"]
            drug = kg_edge["subject"]
            parsed[year].append((contrib, drug))

    # Sort based on abs
    sorted_parsed = {}
    for year, res in parsed.items():
        sorted_parsed[year] = sorted(res, key=lambda x: abs(x[0]), reverse=True)
    return sorted_parsed

def main():
    if not os.path.exists(RESULTS_FILENAME) or RUN:
        results = {}
        for query_json in glob.glob('dw_*.json'):
            # Load in query
            with open(query_json, 'r') as f_:
                query = json.load(f_)
            query = query["message"]
            print(json.dumps(query, indent=2))
            interface = TrapiInterface(query=query)
            interface.build_chp_queries()
            interface.run_chp_queries()
            res = interface.construct_trapi_response()
            results[query_json] = res
        # Write File
        with open(RESULTS_FILENAME, 'wb') as f_:
            pickle.dump(results, f_)
    else:
        #Load results file
        with open(RESULTS_FILENAME, 'rb') as f_:
            results = pickle.load(f_)

    # Parse results
    parsed = parse_results(results)

    with open('parsed-{}'.format(RESULTS_FILENAME), 'wb') as f_:
        pickle.dump(parsed, f_)

if __name__ == '__main__':
    main()
