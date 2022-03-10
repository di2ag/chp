from chp.trapi_interface import TrapiInterface
from chp.apps import *
from collections import defaultdict
import time
from trapi_model.biolink.constants import *
import json

def get_app_config(query):
    return ChpApiConfig

def get_trapi_interface(chp_config=None):
    if chp_config is None:
        chp_config = ChpApiConfig
    return TrapiInterface(
        hosts_filename=chp_config.hosts_filename,
        num_processes_per_host=chp_config.num_processes_per_host,
        bkb_handler=chp_config.bkb_handler,
        joint_reasoner=chp_config.joint_reasoner,
        dynamic_reasoner=chp_config.dynamic_reasoner,
        )

def get_curies():
    interface = get_trapi_interface()
    return interface.get_curies()

def get_meta_knowledge_graph():
    interface = get_trapi_interface()
    return interface.get_meta_knowledge_graph()

def get_response(consistent_queries):
    """ Should return app responses plus app_logs, status, and description information.
    """
    app_logs:list = []
    try:
        interface_dict:defaultdict = setup_queries_based_on_disease_interfaces(consistent_queries)
    except ValueError as ex:
        responses = []
        status = 'Bad request. See description.'
        description = 'Problem during setup. ' + str(ex)
        return responses, app_logs, status, description

    if len(interface_dict.keys()) == 0:
        responses = []
        status = 'Bad request. See description.'
        description = 'Disease Curies are invalid or can not be handled.'
        return responses, app_logs, status, description

    # Setup for CHP inferencing
    try:
        setup_time = time.time()
        for interface, queries in interface_dict.items():
            interface.setup_trapi_queries(queries)
        logger.info('Trapi Interface setup time: {} seconds.'.format(time.time() - setup_time))
    except Exception as ex:
        responses = []
        status = 'Bad request. See description.'
        description = 'Problem during interface setup. ' + str(ex)
        return responses, app_logs, status, description
    
    # Build CHP queries
    try:
        build_time = time.time()
        for interface in interface_dict:
            interface.build_chp_queries()
        logger.info('CHP query build time: {} seconds.'.format(time.time() - build_time))
    except Exception as ex:
        # Add logs from interfaces level
        for interface in interface_dict:
            app_logs.extend(interface.logger.to_dict())
        response = []
        status = 'Bad request. See description.'
        description = 'Problem during CHP query building. '+ str(ex)
        return response, app_logs, status, description

    logger.info('Built Queries.')
    # Run queries
    try:
        reasoning_start_time = time.time()
        for interface in interface_dict:
            interface.run_chp_queries()
        logger.info('Completed Reasoning in {} seconds.'.format(time.time() - reasoning_start_time))
    except Exception as ex:
        # Add logs from interfaces level
        for interface in interface_dict:
            app_logs.extend(interface.logger.to_dict())
        response = []
        status = 'Unexpected error. See description.'
        description = 'Problem during reasoning. ' + str(ex)
        # Report critical error to logs
        logger.critical('Error during reasoning. Check query: {}'.format(query_copy.id))
        return response, app_logs, status, description

    # Construct Response
    responses = []
    for interface in interface_dict:
        responses.extend(interface.construct_trapi_responses())
    
    # Check if any responses came back
    if len(responses) == 0:
        # Add logs from interfaces level
        for interface in interface_dict:
            app_logs.extend(interface.logger.to_dict())
        status = 'No results.'
        return responses, app_logs, status, None

    # Collect app logs from interfaces level
    for interface in interface_dict:
        app_logs.extend(interface.logger.to_dict())

    # Return successful status
    return responses, app_logs, 'Success', None


def get_disease_nodes(query):
    disease_node_ids = query.message.query_graph.find_nodes(categories=[BIOLINK_DISEASE_ENTITY])
    if disease_node_ids is not None:
        disease_nodes = [query.message.query_graph.nodes[_id] for _id in disease_node_ids]
        return disease_nodes
    return None

def get_disease_specific_config(disease_node):
    if disease_node.ids is None:
        return None
        #raise ValueError('Do not support Disease wildcards. Must specify a disease curie.')
    if disease_node.ids[0] == 'MONDO:0005061':
        chp_config = ChpLungApiConfig
    elif disease_node.ids[0] == 'MONDO:0001657':
        chp_config = ChpBrainApiConfig
    elif disease_node.ids[0] == 'MONDO:0007254':
        chp_config = ChpBreastApiConfig
    else:
        return None
        #chp_config = ChpApiConfig
    return chp_config

def setup_queries_based_on_disease_interfaces(consistent_queries):
    config_dict = defaultdict(list)
    for consistent_query in consistent_queries:
        disease_nodes = get_disease_nodes(consistent_query)
        if disease_nodes is not None:
            if len(disease_nodes) > 1:
                consistent_query.info('Using {} config.'.format(ChpApiConfig.name))
                config_dict[ChpApiConfig].append(consistent_query)
                continue
            chp_config = get_disease_specific_config(disease_nodes[0])
            if chp_config is not None:
                consistent_query.info('Using {} config.'.format(chp_config.name))
                config_dict[chp_config].append(consistent_query)
        else:
            config_dict[ChpApiConfig].append(consistent_query)
            continue
    interface_dict = {get_trapi_interface(chp_config): _queries for chp_config, _queries in config_dict.items()}
    return interface_dict
