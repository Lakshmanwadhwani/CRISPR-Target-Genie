import unittest
from pam_input import validate_pam_sequence

class TestPAMInput(unittest.TestCase):

    def test_validate_pam_sequence(self):
        # Test cases for the validate_pam_sequence function

        # Testing valid PAM sequences
        self.assertTrue(validate_pam_sequence("NGG"))
        self.assertTrue(validate_pam_sequence("NAG"))
        self.assertTrue(validate_pam_sequence("NGA"))

        # Testing invalid PAM sequences
        self.assertFalse(validate_pam_sequence("XYZ"))
        self.assertFalse(validate_pam_sequence("N"))
        self.assertFalse(validate_pam_sequence(""))

        # Testing case insensitivity
        self.assertTrue(validate_pam_sequence("ngg"))

if __name__ == '__main__':
    unittest.main()
