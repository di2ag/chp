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
    def __init__(self,source_ara, json_query=None, dict_query=None, hosts_filename=None, num_processes_per_host=0):

        self.integrator = source_ara
        self.json_query = json_query
        self.hosts_filename = hosts_filename
        self.num_processes_per_host = num_processes_per_host
        if dict_query is None:
            self.query = json.loads(json_query)
        else:
            self.query = dict_query
        self.handler = self.getHandler()

    def checkQuery(self):
        return True

    def checkWildCardQuery(self):
        """ Checks if the query graph has a wildcard (no curie specified) for a given node.

        Currently only supports single gene wildcard queries.
        """
        qg = self.query["query_graph"]
        for node in qg["nodes"]:
            if node["type"] == 'gene' and 'curie' not in node:
                return True
        return False

    def getHandler(self):
        if self.integrator == 'default':
            if self.checkWildCardQuery():
                handler = WildCardHandler(self.query, hosts_filename=self.hosts_filename, num_processes_per_host=self.num_processes_per_host)
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
