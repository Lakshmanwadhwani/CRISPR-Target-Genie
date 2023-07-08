import os
import tempfile
import unittest
from .crispr_scoring_moduleA import CRISPRScoringModule, validate_pam_sequence

class TestCRISPRScoringModule(unittest.TestCase):

    def setUp(self):
        # Use a temporary file for testing
        self.fasta_file = tempfile.mktemp(suffix=".fasta")
        with open(self.fasta_file, 'w') as f:
            f.write(">Test\nATGCATGCATGCGGATGCATGCATGC")

        self.pam_sequence = 'CGG'
        self.output_file = tempfile.mktemp(suffix=".txt")
        self.crispr_tool = CRISPRScoringModule(self.fasta_file, self.pam_sequence)

    def test_find_target_sites(self):
        self.crispr_tool.find_target_sites()
        # There should be one target site with 'CGG'
        self.assertEqual(len(self.crispr_tool.target_sites), 1)

    def test_score_target_sites(self):
        self.crispr_tool.find_target_sites()
        self.crispr_tool.score_target_sites()
        # Check if score has been computed
        self.assertIsNotNone(self.crispr_tool.target_sites[0]['score'])

    def test_validate_pam_sequence(self):
        # Test valid PAM sequence
        self.assertTrue(validate_pam_sequence('CGG'))
        # Test invalid PAM sequence
        self.assertFalse(validate_pam_sequence('XYZ'))

    def test_write_results_to_file(self):
        self.crispr_tool.find_target_sites()
        self.crispr_tool.score_target_sites()
        self.crispr_tool.write_results_to_file(self.output_file)

        # Test if the file exists
        self.assertTrue(os.path.exists(self.output_file))

        # Test if the file has correct content
        with open(self.output_file, 'r') as f:
            lines = f.readlines()
        self.assertTrue("Start\tEnd\tSequence\tScore\n" in lines)
        self.assertTrue("\tCGG\t" in "".join(lines))

    def tearDown(self):
        # Clean up temporary files after test
        os.remove(self.fasta_file)
        if os.path.exists(self.output_file):
            os.remove(self.output_file)


if __name__ == "__main__":
    unittest.main()
list
