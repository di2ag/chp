import sys

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
fused_bkb.load('/home/public/data/ncats/BabelBKBs/collapsed50/fusion.bkb')

#-- Here are the associated patient data files
patient_data_file = '/home/public/data/ncats/BabelBKBs/collapsed50/patient_data.pk'
withheld_patients_file = '/home/public/data/ncats/BabelBKBs/collapsed50/withheldPatients.csv'

#-- Instiante reasoner
reasoner = Reasoner(fused_bkb=fused_bkb)

#-- Set the patient data file 
reasoner.set_src_metadata(patient_data_file)
#reasoner.collapsed_bkb.makeGraph()

#-- If you want to see what genetic or demographic evidence is avaliable, uncomment the line below
#print(reasoner.metadata_ranges)

#-- Make a query (evidence is for genetic info, and meta_ is for demographic info)
query0 = Query(evidence={'_mut_RPTN': 'True'},
               targets=list(),
               meta_evidence=[('Age_of_Diagnosis', '>=',20000)],
               meta_targets=[('Survival_Time', '>=', 300)],
               interpolation='correlation')

#-- Run the query.
query = reasoner.analyze_query(query0, check_mutex=False)

#-- Return the report
query.getReport()
#print(query.result.print_contributions())
print(query.result.completed_inferences_report())
#query.bkb.makeGraph()

#-- Check for mutex if you want to.
#print(checkMutex(query.bkb))
