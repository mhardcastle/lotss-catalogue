"""
Tests to be run on the catalogues.
The location of the source catalogue and the component catalogue can be
entered as environment variables called 'SOURCE_CATALOGUE' and 
'COMPONENT_CATALOGUE' or as the last two arguments of the call to the test 
program.
Examples:
$ SOURCE_CATALOGUE=cat1.fits COMPONENT_CATALOGUE=cat2.fits \
 python test_source_catalogue.py
$ python test_source_catalogue.py cat1.fits cat2.fits
TODO: Improve the argument parsing
"""
import unittest
import os
from astropy.table import Table, join
import numpy as np


class TestCatalogues(unittest.TestCase):
    """
    Tests on the catalogues
    """
    SOURCE_CATALOGUE = os.environ.get('SOURCE_CATALOGUE', "")
    COMPONENT_CATALOGUE = os.environ.get('COMPONENT_CATALOGUE', "")

    @classmethod
    def setUpClass(cls):
        """
        Read the catalogues (expensive) once before running the tests.
        """
        cls.source = Table.read(cls.SOURCE_CATALOGUE)
        cls.component = Table.read(cls.COMPONENT_CATALOGUE)

    def test_duplicate_opticalID(self):
        """
        Test for duplicates in the ID_name column of the source catalogue 
        excluding IDs designated as 'Mult'
        """
        oids = self.source['ID_name']
        #oids = np.char.rstrip(oids)
        unique, counts = np.unique(oids, return_counts=True)
        duplicates = unique[(counts > 1) & (unique != "Mult")]
        self.assertEqual(len(duplicates), 0, msg="Duplicated optical IDs found")
    
    def test_duplicate_sourceName(self):
        """
        Test for duplicates in the Source_Name column of the source catalogue
        """
        sns = self.source['Source_Name']
        unique, counts = np.unique(sns, return_counts=True)
        duplicates = unique[(counts > 1)]
        self.assertEqual(len(duplicates), 0, 
                         msg="Duplicated source names found")
    
    def test_duplicated_components(self):
        """
        Test for duplicates in the Source_Name column
        """
        cids = self.component['Component_Name']
        unique, counts = np.unique(cids, return_counts=True)
        duplicates = unique[(counts > 1)]
        self.assertEqual(len(duplicates), 0, 
                         msg="Duplicated component names found")

    def test_match_catalogues_source(self):
        """
        Test the match between catalogues. Check entries that are exclusively 
        on the source catalogue
        """
        x_sources = np.setdiff1d(self.source["Source_Name"], 
                                 self.component["Source_Name"])
        self.assertEqual(len(x_sources), 0, 
                         msg="Unmatched entries found in the source catalogue")
    
    def test_match_catalogues_component(self):
        """
        Test the match between catalogues. Check entries that are exclusively 
        on the component catalogue
        """
        x_components = np.setdiff1d(self.component["Source_Name"], 
                                    self.source["Source_Name"])
        self.assertEqual(len(x_components), 0, 
                         msg="Unmatched entries found in the component catalogue")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        TestCatalogues.COMPONENT_CATALOGUE = sys.argv.pop()
        TestCatalogues.SOURCE_CATALOGUE = sys.argv.pop()
    unittest.main()