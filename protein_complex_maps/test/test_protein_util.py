#!/usr/bin/python

# Unit tests for protein util

import unittest
import requests
import protein_complex_maps.util.protein_util as pu
import numpy as np


#claude: these tests query the live UniProt REST API; skip (rather than error) when it is unreachable
def uniprot_available():
    try:
        requests.head("https://rest.uniprot.org/", timeout=10)
        return True
    except requests.exceptions.RequestException:
        return False


@unittest.skipUnless(uniprot_available(), "UniProt REST API (rest.uniprot.org) is unreachable; skipping network tests")
class ProteinUtilTest(unittest.TestCase):

    def setUp(self,):
        self.protein_id = ["Q6N089","Q53HL2"]

    def testLength(self,):
        length = pu.get_length_uniprot(self.protein_id)
        assert( length["Q6N089"] == 472 )
        assert( length["Q53HL2"] == 280 )

if __name__ == "__main__":
	unittest.main()
