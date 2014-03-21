#!/usr/bin/env python

__author__ = "Sam Way"
__copyright__ = "Copyright 2013, The American Gut Project"
__credits__ = ["Sam Way"]
__license__ = "BSD" 
__version__ = "unversioned"
__maintainer__ = "Sam Way"
__email__ = "samuel.way@colorado.edu"

from ..lib.parse import parse_mapping_file_to_dict, get_filtered_taxa_summary
from unittest import TestCase, main

TEST_MAPPING_FILE = './files/test_mapping.txt'

class test_mapping_file_parse(TestCase):
    def setUp(self):
        self.categories = ['TEST_CATEGORY', 'AWESOME_CATEGORY']
        self.sample_ids = ['sample_a', 'sample_b']
        self.metadata_dict = {'sample_a':{'TEST_CATEGORY':'1', 'AWESOME_CATEGORY':'super'},
                             'sample_b':{'TEST_CATEGORY':'2', 'AWESOME_CATEGORY':'totally'}}

    def test_mapping_file(self):
        mapping_dict, comments = parse_mapping_file_to_dict(open(TEST_MAPPING_FILE,'rU'))
        for sample_id, sample_dict in mapping_dict.iteritems():

            # Does the sample dictionary contain all of the metadata categories?
            self.assertEqual(set(sample_dict.keys()), set(self.categories))
    
            # Are all metadata values correct?
            for category in sample_dict.keys():
                self.assertEqual( sample_dict[category], self.metadata_dict[sample_id][category] )

if __name__ == '__main__':
    main()

