def get_pam_sequence():
    """
    Function to get PAM sequence input from the user
    """
    pam = input("Enter the PAM sequence: ")
    return pam.upper()

def validate_pam_sequence(pam):
    """
    Function to validate if the input PAM sequence is valid
    """
 # Convert the input to uppercase before validating it
    pam = pam.upper()   
    # Check if the PAM sequence contains only 'A', 'T', 'C', 'G', 'N'
    for base in pam:
        if base not in ['A', 'T', 'C', 'G', 'N']:
            return False
    # Check if the length of PAM is appropriate (usually 2-6 bases)
    if 2 <= len(pam) <= 6:
        return True
    else:
        return False

if __name__ == "__main__":
    pam_sequence = get_pam_sequence()
    if validate_pam_sequence(pam_sequence):
        print(f"The PAM sequence '{pam_sequence}' is valid.")
    else:
        print(f"The PAM sequence '{pam_sequence}' is invalid.")
