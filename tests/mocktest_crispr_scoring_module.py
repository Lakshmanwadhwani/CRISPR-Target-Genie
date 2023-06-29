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
            # Check the start and end positions of the first target site
            first_target_site = crispr_tool.target_sites[3]
            self.assertEqual(first_target_site['start'], 3)  # Check if start position is expected
            self.assertEqual(first_target_site['end'], 5)    # Check if end position is expected
            second_target_site = crispr_tool.target_sites[57]
            self.assertEqual(second_target_site['start'], 57)  # Check if start position is expected
            self.assertEqual(second_target_site['end'], 59)    # Check if end position is expected
            self.assertEqual(first_target_site['score'], 0.6666666666666666)  # Check if score is as expected


if __name__ == '__main__':
    unittest.main()

