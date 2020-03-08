import os
import sys

from pybkb import bayesianKnowledgeBase as BKB

#-- Change to your local data folder
NCATS_DIR = '/home/cyakaboski/src/python/pojects/bkb-pathway-provider'

sys.path.append(os.path.join(NCATS_DIR, 'core'))

from query import Query
from reasoner import Reasoner

if __name__ == '__main__':
    bkb = BKB()
    bkb.load('/home/public/data/ncats/fusedPatientBKBSmall/fusion.bkb')

    reasoner = Reasoner(bkb)

    query0 = Query(evidence={'_sourceFusion_14_02_2020_16:35:04__mut_CACNA1H=': '[0]_-7750956075882572850'},
                   targets=['mut_CACNA1H=', 'mu-STD>=CACNA1H<=mu+STD='],
                   type='updating'
                  )

    result = reasoner.analyze_query(query0)
