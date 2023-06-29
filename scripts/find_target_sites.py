from Bio import SeqIO

def validate_pam_sequence(pam):
    """
    Function to validate if the input PAM sequence is valid
    """
    for base in pam.upper():
        if base not in ['A', 'T', 'C', 'G', 'N']:
            return False
    return 2 <= len(pam) <= 6


def find_target_sites(gene_sequence, pam_sequence):
    """
    Function to find the target sites in a gene sequence based on PAM sequence
    """
    target_sites = []
    sequence = str(gene_sequence).upper()
    pam = pam_sequence.upper()

    # Sliding window through the gene_sequence
    for i in range(len(sequence) - len(pam) + 1):
        if sequence[i:i + len(pam)] == pam:
            target_sites.append((i, i + len(pam) -1))

    return target_sites


def main(fasta_file, pam_sequence):
    # Validate PAM sequence
    if not validate_pam_sequence(pam_sequence):
        print("Invalid PAM sequence")
        return

    # Read gene sequences from the FASTA file
    records = SeqIO.parse(fasta_file, "fasta")

    # Find and print the target sites for each gene sequence
    for record in records:
        target_sites = find_target_sites(record.seq, pam_sequence)
        print(f"Target sites for sequence {record.id}: {target_sites}")


if __name__ == "__main__":
    fasta_file = r"C:\Users\Admin\OneDrive\Desktop\p53fasta"  # Replace this with the actual path
    pam_sequence = "NGG"  # or any other PAM sequence you want to use
    main(fasta_file, pam_sequence)
