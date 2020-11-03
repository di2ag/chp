import requests
import json
import unittest
import logging

API_ADDRESS = 'http://localhost:80'

logger = logging.getLogger(__name__)
logger.warning('Make sure your django server is running at: {}.'.format(API_ADDRESS))

class testEndpoint(unittest.TestCase):
    def setUp(self):
        self._default_url = API_ADDRESS
        self._query_endpoint = '/query/'
        self._predicates_endpoint = '/predicates/'
        self._curies_endpoint = '/curies/'

    def test_predicates(self):
        _url = self._default_url + self._predicates_endpoint
        logger.info('Testing endpoint: {}'.format(_url))
        res = requests.get(_url)
        print(res.content)
        ret = res.json()
        print(ret)


    def test_curies(self):
        _url = self._default_url + self._curies_endpoint
        logger.info('Testing endpoint: {}'.format(_url))
        res = requests.get(_url)
        ret = res.json()
        print(ret)
