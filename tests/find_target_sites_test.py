import unittest
from find_target_sites import find_target_sites  # replace your_script_name with the actual script name where your functions are


class TestFindTargetSites(unittest.TestCase):

    def test_find_target_sites(self):
        # Test case 1: standard case
        gene_sequence = "ATGGGCGTAGCTAGGCGTAG"
        pam_sequence = "GGC"
        expected_output = [(3, 5), (13, 15)]
        self.assertEqual(find_target_sites(gene_sequence, pam_sequence), expected_output)

        # Test case 2: case with 'N' in PAM sequence representing any nucleotide
        gene_sequence = "ATGGGCGTAGCTTTANGGCGTAG"
        pam_sequence = "NGG"
        expected_output = [(15, 17)]
        self.assertEqual(find_target_sites(gene_sequence, pam_sequence), expected_output)

        # Test case 3: case with no matches
        gene_sequence = "ATGCACTGACTG"
        pam_sequence = "GGC"
        expected_output = []
        self.assertEqual(find_target_sites(gene_sequence, pam_sequence), expected_output)

        # Test case 4: case with lowercase input
        gene_sequence = "aatgggcgtagctaggcgtag"
        pam_sequence = "ggc"
        expected_output = [(4, 6), (14, 16)]
        self.assertEqual(find_target_sites(gene_sequence, pam_sequence), expected_output)


if __name__ == '__main__':
    unittest.main()
