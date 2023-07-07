from crispr_scoring_module import CRISPRScoringModule, validate_pam_sequence


def test_validate_pam_sequence():
    assert validate_pam_sequence('CGG')
    assert not validate_pam_sequence('INVALID')


# Add more tests here for CRISPRScoringModule...
