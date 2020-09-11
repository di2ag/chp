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

#-- Integrators
from chp.integrator.exploring_agent import ExploringHandler
from chp.integrator.unsecret_agent import UnsecretHandler
from chp.integrator.ranking_agent import RankingHandler

class ReasonerStdHandler:
    def __init__(self, source_ara, json_query=None, dict_query=None, hosts_filename=None, num_processes_per_host=0):
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

    def getHandler(self):
        if self.integrator == 'exploring':
            handler = ExploringHandler(self.query, hosts_filename=self.hosts_filename, num_processes_per_host=self.num_processes_per_host)
        elif self.integrator == 'unsecret':
            handler = UnsecretHandler(self.query, hosts_filename=self.hosts_filename, num_processes_per_host=self.num_processes_per_host)
        elif self.integrator == 'ranking':
            handler = RankingHandler(self.query, hosts_filename=self.hosts_filename, num_processes_per_host=self.num_processes_per_host)
        else:
            raise NotImplementedError('Not integrated with {}.'.format(self.integrator))
        return handler

    def buildChpQueries(self):
        return self.handler.buildQueries()

    def runChpQueries(self):
        return self.handler.runQueries()

    def constructDecoratedKG(self):
        return self.handler.constructDecoratedKG()
