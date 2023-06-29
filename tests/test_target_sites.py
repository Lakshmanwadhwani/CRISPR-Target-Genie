import unittest
from target_sites import find_target_sites 

class TestTargetSites(unittest.TestCase):

        def test_find_target_sites(self):
        # Test case 1: standard case
        gene_sequence = "AGTCGATCGTACGAGACGGTGCTAGTCACTCGGCGTATGCATCGTGAGCTACGGATCGACGACGTCGATCGGACGTACGTAGCGTACGTACGTACGGTTT"
        pam_sequence = "GGC"
        # Correcting the expected_output to [(4, 7), (16, 19)]
        expected_output = [(16, 18), (30, 33), (51, 54), (69, 72), (93, 96)]
        self.assertEqual(find_target_sites(gene_sequence, pam_sequence), expected_output)

        # Test case 2: No target sites
        gene_sequence = "ATTGGGATGCTTAGTTCCCTGAACGCATTCGA"
        pam_sequence = "CGG"
        expected_output = []
        self.assertEqual(find_target_sites(gene_sequence, pam_sequence), expected_output)

        # Test case 3: Testing with lowercase
        gene_sequence = "atgggcgtagctaggcgcag"
        pam_sequence = "ggc"
        expected_output = [(3, 6), (13, 16)]
        self.assertEqual(find_target_sites(gene_sequence, pam_sequence), expected_output)

        # Add more test cases as needed

if __name__ == '__main__':
    unittest.main()

