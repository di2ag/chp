import sys
import os

from pybkb import bayesianKnowledgeBase as BKB
from pybkb.core.python_base.reasoning import checkMutex

#!!!!!!!!!! Change to your local bkb-pathway-core directory.
#sys.path.append('/home/ghyde/bkb-pathway-provider/core')
sys.path.append('/home/cyakaboski/src/python/projects/bkb-pathway-provider/core')

from reasoner import Reasoner
from query import Query

#-- Initalize a BKB
fused_bkb = BKB()

#-- Load in the fused bkb from our datafiles
fused_bkb.load('/home/public/data/ncats/AxleBKBS/AxleDemoBKB3-14-20/fusion.bkb')

#-- Here are the associated patient data files
patient_data_file = '/home/public/data/ncats/AxleBKBS/AxleDemoBKB3-14-20/patient_data.pk'
withheld_patients_file = '/home/public/data/ncats/AxleBKBS/AxleDemoBKB3-14-20/withheldPatients.csv'
gene_var_direct_file = '/home/public/data/ncats/AxleBKBS/AxleDemoBKB3-14-20/geneFreqAnal.csv'

#-- Instaniate reasoner
reasoner = Reasoner(fused_bkb=fused_bkb, gene_var_direct=gene_var_direct_file, max_new_ev=2)

#-- Set the patient data file 
reasoner.set_src_metadata(patient_data_file)
#reasoner.collapsed_bkb.makeGraph()

#-- If you want to see what genetic or demographic evidence is avaliable, uncomment the line below
#print(reasoner.metadata_ranges)

#-- Make a query (evidence is for genetic info, and meta_ is for demographic info)
query0 = Query(name='sample_independ',
               evidence={"var_5'FLANK=": "True",
                         "var_3'UTR=": "True"},
               targets=None,
               meta_evidence=[('Age_of_Diagnosis', '>=', 20000)],
               meta_targets=[('Survival_Time', '>=', 943)])

#-- Run the query.
query = reasoner.analyze_query(query0, check_mutex=False, target_strategy='explicit', interpolation='independence')

#-- Return the report
query.getReport()
query.save(os.getcwd())
