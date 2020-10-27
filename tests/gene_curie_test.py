import unittest
from chp.gene_curie import curie_match, spinoff_requests, get_gene_list


# Tests chp.gene_curie.py to make sure we get solid request returns. Also
# tests that the program is thread safe meaning all all gene requests are returned
# with their respective curie IDs rather than another gene's curie ID

class testGeneCurie(unittest.TestCase):
    #self.gene_file = '/home/public/data/ncats/data_drop_03-04-2020/wxs.csv'
    #self.gene_list = get_gene_list(self.gene_file)

    # test a good gene request and ensure the ensembl ID matches a known ensemble ID
    def test_curie_request(self):
        raf_curie = 'ENSG00000132155'
        gene = 'RAF1'
        curie_info = curie_match(gene)
        self.assertEqual(curie_info[1], raf_curie)

    # validating that a bunk gene request returns 'not found'
    def test_bunk_curie_request(self):
        bunk_gene = 'Bunk'
        curie_info = curie_match(bunk_gene)
        self.assertEqual(curie_info[1], 'not found')

    # make sure there are no anomolies with linking genes back up with multithreading
    def test_safe_queue_results(self):
        max_threads = 4
        gene_symbol_tuples = {'RAF1': 'ENSG00000132155',
                              'TTN': 'ENSG00000155657',
                              'ZFP2': 'ENSG00000198939',
                              'CMAS': 'ENSG00000111726',
                              'HRAS': 'ENSG00000174775',
                              'DAG1': 'ENSG00000173402',
                              'PPBP': 'ENSG00000163736',
                              'BRAF': 'ENSG00000157764',
                              'BMX': 'ENSG00000102010'}
        gene_list = list(gene_symbol_tuples.keys())
        results = spinoff_requests(gene_list, max_threads)

        for r in results:
            self.assertEqual(r[1], gene_symbol_tuples[r[0]])

if __name__ == '__main__':
    unittest.main()
