'''
Author: Sara Tahir
Date: 04/05/2022
Purpose: Tests for frameshift_detector.y
Example usage: python3 frameshift_detector_test.py 
'''

import unittest
from Bio import SeqIO
import argparse

from frameshift_detector import  GenomeFeature, find_downstream_frameshift, find_upstream_frameshift

genome_features = []
heptamers = [['CTTAGGC', 10]]
stop_codons = ["TAG","TAA","TGA"]
frame = 1
search_limits = 10000
parser = argparse.ArgumentParser(description='Find some frameshifts')

parser.set_defaults(verbose=True)
args = parser.parse_args()

class TestFrameshiftDetector(unittest.TestCase):
    def test_downstream_search_one(self):
        frameshifts = find_downstream_frameshift(genome_features[0],frame,search_limits,stop_codons, heptamers, args)
        print(frameshifts[0].frameshift_seq)
        self.assertEqual(str(frameshifts[0].genome_feature.location), '[0:18](+)')
        self.assertEqual(frameshifts[0].frameshift_seq, 'ATGCTTGGCTGA')

    def test_downstream_search_two(self):
        frameshifts = find_downstream_frameshift(genome_features[1],frame,search_limits,stop_codons, heptamers, args)
        self.assertEqual(0, len(frameshifts))

    def test_downstream_search_three(self):
        frameshifts = find_downstream_frameshift(genome_features[2],frame,search_limits,stop_codons, heptamers, args)
        print(frameshifts[0].frameshift_seq)
        self.assertEqual('atgcctataactaccctaacacagcccaaacttggccattcaaccataccactccgtgataaatga', frameshifts[0].frameshift_seq.lower())


if __name__ == '__main__':
    testing_data_path = 'testing_data/test.gbff'
    nucrecs = SeqIO.parse(testing_data_path, "genbank")
 
    for nucrec in nucrecs:
        nucrec_id = nucrec.id
        print('Processing: ', nucrec_id)
        nucrec_desc = nucrec.description
        #print("Processing: " + nucrec.id)
        # For each feature in the genebank file, create a GenomeFeature object and store in global genome_features array
        for i in range(len(nucrec.features)):
            feat = nucrec.features[i]
            if (feat.type == 'CDS'):
                if (feat.strand == 1):
                    genome_features.append(GenomeFeature(nucrec=nucrec, nuc_acc=nucrec_id, nuc_desc=nucrec_desc, feature=feat, strand='+1', species='testing'))

    unittest.main()
    