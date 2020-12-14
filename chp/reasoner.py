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
        """
        self.bkb_handler = bkb_handler
        self.hosts_filename = hosts_filename
        self.num_processes_per_host = num_processes_per_host
        self.venv = venv

        # Run base reasoner setup
        self._setup_base_reasoner()

    def _setup_base_reasoner(self):
        # Read in raw patient data
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
