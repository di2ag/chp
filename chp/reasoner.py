import pickle
import logging

from pybkb.python_base.reasoning.reasoning import updating
from pybkb.python_base.reasoning.joint_reasoner import JointReasoner
from pybkb.python_base.learning.bkb_builder import LinkerBuilder

from chp_data.patient_bkb_builder import PatientBkbBuilder

from chp.mixins.reasoner.chp_joint_reasoner_mixin import ChpJointReasonerMixin
from chp.mixins.reasoner.chp_dynamic_reasoner_mixin import ChpDynamicReasonerMixin

logger = logging.getLogger(__name__)
#logger.setLevel(logging.INFO)

# Base Reasoner Class

class BaseReasoner:
    def __init__(self,
                 bkb_handler,
                 hosts_filename=None,
                 num_processes_per_host=None,
                 venv=None,
                 patient_bkb_builder=None,
                 gene_prelinked_bkb_override=None,
                 drug_prelinked_bkb_override=None,
                ):
        """ The base reasoner class for CHP.

            :param bkb_handler: The CHP Data handler that holds the paths to all important BKB files.
            :type bkb_handler: chp_data.bkb_handler.BkbHandler
            :param hosts_filename: If using distributed reasoning, this is the filename that contains
                all the IP addresses to run reasoning.
            :type hosts_filename: str
            :param num_processes_per_host: The number of processes that can be run on each host.
            :type num_processes_per_host: int
            :param venv: Virtual env path where to run distrubted reasoning.
            :type venv: str
            :param patient_bkb_builder: Chp Data's Patient BKB Builder that has already been initialized by some
            some other means. Used primarily in obtaining reasoning results for internal analysis.
            :type patient_bkb_builder: chp_data.PatientBkbBuilder
            :param gene_prelinked_bkb_override: A gene bkb that is used to override the bkb loaded from the
            bkb handler. Used primarily in obtaining reasoning results for internal analysis.
            :type gene_prelinked_bkb_override: pybkb.bayesianKnowledgeBase
            :param drug_prelinked_bkb_override: A drug bkb that is used to override the bkb loaded from the
            bkb handler. Used primarily in obtaining reasoning results for internal analysis.
            :type drug_prelinked_bkb_override: pybkb.bayesianKnowledgeBase
        """
        self.bkb_handler = bkb_handler
        self.hosts_filename = hosts_filename
        self.num_processes_per_host = num_processes_per_host
        self.venv = venv
        self.patient_bkb_builder = patient_bkb_builder
        self.gene_prelinked_bkb_override = gene_prelinked_bkb_override
        self.drug_prelinked_bkb_override = drug_prelinked_bkb_override

        # Run base reasoner setup
        self._setup_base_reasoner()

    def _setup_base_reasoner(self):
        """ Reads in patient BKB builder if one is not passed and loads patient data into
        the appropriate format. Lastly, it runs the specific reasoner setup method.
        """
        # Read in raw patient data
        if self.patient_bkb_builder is None:
            with open(self.bkb_handler.patient_data_pk_path, 'rb') as patient_file:
                self.raw_patient_data = pickle.load(patient_file)
            # Load in the CHP Data Patient data builder
            self.patient_bkb_builder = PatientBkbBuilder(
                                self.raw_patient_data,
                                self.bkb_handler,
                               )
            logger.info('Constructed Patient Bkb Builder.')
        # For readability
        self.patient_data = self.patient_bkb_builder.patient_data
        features = self.patient_data[list(self.patient_data.keys())[0]].keys()
        non_gene_features = [feature for feature in features if 'ENSEMBL' not in feature]

        # Setup reasoner mixin
        self._setup_reasoner()

    def _setup_reasoner(self):
        pass

    def run_query(self, query):
        pass


# Mixin classes

class ChpJointReasoner(ChpJointReasonerMixin, BaseReasoner):
    pass

class ChpDynamicReasoner(ChpDynamicReasonerMixin, BaseReasoner):
    pass
