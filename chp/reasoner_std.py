'''
Source code developed by DI2AG.
Thayer School of Engineering at Dartmouth College
Authors:    Dr. Eugene Santos, Jr
            Mr. Chase Yakaboski,
            Mr. Gregory Hyde,
            Dr. Keum Joo Kim
'''


import json
import itertools
import tqdm
import numpy as np
import logging
import csv

#-- Integrators
from chp.integrator.default_handler import DefaultHandler
from chp.integrator.exploring_agent import ExploringHandler
from chp.integrator.unsecret_agent import UnsecretHandler
from chp.integrator.ranking_agent import RankingHandler
from chp.integrator.relay9_22 import Relay9_22
from chp.integrator.wildcard_handler import WildCardHandler

#-- Setup logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

class ReasonerStdHandler:
    def __init__(self,
                 source_ara,
                 json_query=None,
                 dict_query=None,
                 hosts_filename=None,
                 num_processes_per_host=0,
                 max_results=10):
        self.integrator = source_ara
        self.json_query = json_query
        self.hosts_filename = hosts_filename
        self.num_processes_per_host = num_processes_per_host
        self.max_results = max_results
        if json_query is not None or dict_query is not None:
            if dict_query is None:
                self.query = json.loads(json_query)
            else:
                self.query = dict_query
        else:
            self.query = None
        self.handler = self.getHandler()

    def getCuries(self):
        """ Returns the avaliable curies and there associated names.
        """
        with open(self.handler.bkb_data_handler.curies_path, 'r') as curies_file:
            return json.load(curies_file)
        '''
        curies = {
            "gene": [],
            "chemical_substance": [],
            "phenotypic_feature": [],
            "disease": [],
        }

        # Read in the gene curie map
        with open(self.handler.bkb_data_handler.gene_curie_path, 'r') as gene_file:
            reader = csv.reader(gene_file)
            next(reader)
            for row in reader:
                curies["gene"].append({"name": row[0],
                                       "curie": row[1]})
        # Read in the drug curie map
        with open(self.handler.bkb_data_handler.drug_curie_path, 'r') as drug_file:
            reader = csv.reader(drug_file)
            next(reader)
            for row in reader:
                curies["chemical_substance"].append({"name": row[0],
                                                     "curie": row[1]})

        # Add diseases
        curies["disease"].append({"name": 'breast_cancer',
                                  "curie": 'MONDO:0007254'})

        # Add phenotypic features, i.e. outcomes
        curies["phenotypic_feature"].append({"name": 'survival_time',
                                             "curie": 'EFO:0000714'})

        return curies
        '''

    def checkQuery(self):
        return True

    def checkWildCardQuery(self):
        """ Checks if the query graph has a wildcard (no curie specified) for a given node.

        Currently only supports single gene wildcard queries.
        """
        if self.query is not None:
            qg = self.query["query_graph"]
            for _, node in qg["nodes"].items():
                if node["type"] == 'gene' and 'curie' not in node:
                    return True
        return False

    def getHandler(self):
        if self.integrator == 'default':
            if self.checkWildCardQuery():
                handler = WildCardHandler(self.query, hosts_filename=self.hosts_filename, num_processes_per_host=self.num_processes_per_host,
                                          max_results=self.max_results)
            else:
                handler = DefaultHandler(self.query, hosts_filename=self.hosts_filename, num_processes_per_host=self.num_processes_per_host)
        elif self.integrator == 'exploring':
            handler = ExploringHandler(self.query, hosts_filename=self.hosts_filename, num_processes_per_host=self.num_processes_per_host)
        elif self.integrator == 'unsecret':
            handler = UnsecretHandler(self.query, hosts_filename=self.hosts_filename, num_processes_per_host=self.num_processes_per_host)
        elif self.integrator == 'ranking':
            handler = RankingHandler(self.query, hosts_filename=self.hosts_filename, num_processes_per_host=self.num_processes_per_host)
        elif self.integrator == 'relay9_22':
            handler = Relay9_22(self.query, hosts_filename=self.hosts_filename, num_processes_per_host=self.num_processes_per_host)
        else:
            raise NotImplementedError('Not integrated with {}.'.format(self.integrator))
        return handler

    def buildChpQueries(self):
        return self.handler.buildQueries()

    def runChpQueries(self):
        return self.handler.runQueries()

    def constructDecoratedKG(self):
        return self.handler.constructDecoratedKG()
