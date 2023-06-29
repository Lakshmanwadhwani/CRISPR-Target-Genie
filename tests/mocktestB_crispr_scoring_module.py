import unittest
from Bcrispr_scoring_module import CRISPRScoringModule
from tempfile import NamedTemporaryFile


class TestCRISPRScoringModule(unittest.TestCase):
    def test_find_target_sites(self):
        # Mock fasta data
        mock_fasta_data = """\
>mock_sequence
GCACGGCTCTAGTCCCAACGGCGTCTTTCTGAAAAGTGGTAAATGTCACATGAT
CAACGGCACACTAGCAACGGGTCAGGCCATGCATGGTACGGGTCTAAGGCACAC
AAACGGCTGTAGTCCCCTGATGATTACGCATTGCCCTTCCGCGCCTAGAGTGGG
GCCCGGGGCTGAATGCTACCGCATCATAGCGCTCCGCGAATGGCGAAGCAAGGC
AAGCGGAACATCGCGTTAGTTGTTTCAGCAACCTTAGCAGAGATGTCTACATAG
CTCCGGCGTACGGGTGTGTACGAGACTAAAAGGAGGCGTGCAGTTTTCTCGGGC
GATCGGTATTTATGGTTGGGTGCAAGACTTGCCGCCTATCGTTCAACGAGTGGC
AAACGGCATCCAGCTAGCTTACTTAGTCGTCGACATATTACTAGTAGTATAAGG"""

        # Create a temporary file with the mock fasta data
        with NamedTemporaryFile(mode="w+", delete=False) as temp_file:
            temp_file.write(mock_fasta_data)
            temp_file.seek(0)
            # Run CRISPRScoringModule using the temporary file
            crispr_tool = CRISPRScoringModule(temp_file.name)
            crispr_tool.find_target_sites()
            crispr_tool.score_target_sites()

            # Add assertions here to check if the output is as expected
            # Check the position and sequence of the first target site
            first_target_site = crispr_tool.target_sites[0]
            self.assertEqual(first_target_site['sequence'], 'CGG')  # Check if position is expected (1-based indexing)
            self.assertEqual(first_target_site['start'], 3)  # Check if position is expected (1-based indexing)
            self.assertEqual(first_target_site['end'], 6)  # Check if position is expected (1-based indexing)
            second_target_site = crispr_tool.target_sites[1]
            self.assertEqual(second_target_site['sequence'], 'CGG')
            self.assertEqual(second_target_site['start'], 18)
            self.assertEqual(second_target_site['end'], 21)
            # ... additional assertions for the second target site ...

            # Check the position and sequence of the third target site
            third_target_site = crispr_tool.target_sites[2]
            self.assertEqual(third_target_site['sequence'], 'CGG')
            self.assertEqual(third_target_site['start'], 57)
            self.assertEqual(third_target_site['end'], 60)
            # Optionally, you can also check the score
            self.assertGreaterEqual(first_target_site['score'], 0)
            self.assertLessEqual(first_target_site['score'], 1)


if __name__ == '__main__':
    unittest.main()

