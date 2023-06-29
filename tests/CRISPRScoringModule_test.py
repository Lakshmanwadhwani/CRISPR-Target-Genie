import unittest
from unittest.mock import patch, Mock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from crispr_scoring_module import CRISPRScoringModule

class TestCRISPRScoringModule(unittest.TestCase):

    def setUp(self):
        self.fasta_sequences = [
            SeqRecord(Seq("AAACCCGGGTTT"), id="seq1"),
            SeqRecord(Seq("GGGCCCAAATTTCGG"), id="seq2")
        ]
        self.pam_sequence = 'CGG'
        
        # Initialize CRISPRScoringModule with the test data
        self.crispr_tool = CRISPRScoringModule(None, self.pam_sequence)
        self.crispr_tool.target_sites = [
            {'start': 5, 'end': 7, 'sequence': 'CGG', 'score': None},
            {'start': 12, 'end': 14, 'sequence': 'CGG', 'score': None}
        ]
    
    @patch('crispr_scoring_module.CRISPRScoringModule.load_fasta_sequences')
    def test_find_target_sites(self, mock_load_fasta_sequences):
        # Mocking the file reading part
        mock_load_fasta_sequences.return_value = self.fasta_sequences
        
        # Call the method you wish to test (in this example, find_target_sites)
        # The mock_load_fasta_sequences method will now return self.fasta_sequences instead of reading a file
        target_sites = self.crispr_tool.find_target_sites('fake_file.fasta')
        
        # Assert your expected output here
        # For example:
        self.assertIsNotNone(target_sites)
        # Further assertions can be added to check contents of target_sites
    
    def test_score_target_sites(self):
        # Test scoring based on GC content
        self.crispr_tool.score_target_sites()
        
        # Check the scores
        for target_site in self.crispr_tool.target_sites:
            sequence = target_site['sequence']
            expected_score = 100.0 * sequence.count('G') + sequence.count('C') / len(sequence)
            self.assertEqual(target_site['score'], expected_score)

if __name__ == '__main__':
    unittest.main()