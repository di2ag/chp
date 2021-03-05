import json
import glob
import logging

from chp.trapi_interface import TrapiInterface

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def main():
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
    print(results)

if __name__ == '__main__':
    main()
