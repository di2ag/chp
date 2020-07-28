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
from chp.integrator.explorer_agent import ExplorerHandler
from chp.integrator.unsecret_agent import UnsecretHandler
from chp.integrator.ranking_agent import RankingHandler

class ReasonerStdHandler:
    def __init__(self, source_ara, json_query=None, dict_query=None):
        self.integrator = source_ara
        self.json_query = json_query
        if dict_query is None:
            self.query = json.loads(json_query)
        else:
            self.query = dict_query
        self.handler = self.getHandler()

    def checkQuery(self):
        return True

    def getHandler(self):
        if self.integrator == 'explorer':
            handler = ExplorerHandler(self.query)
        elif self.integrator == 'unsecret':
            handler = UnsecretHandler(self.query)
        elif self.integrator == 'ranking':
            handler = RankingHandler(self.query)
        else:
            raise NotImplementedError('Not integrated with {}.'.format(self.integrator))
        return handler

    def buildChpQueries(self):
        return self.handler.buildQueries()

    def runChpQueries(self):
        return self.handler.runQueries()

    def constructDecoratedKG(self):
        return self.handler.constructDecoratedKG()
