from Bio import SeqIO

def parse_fasta(file_path):
    """
    This function parses a FASTA file and returns a dictionary
    where the keys are sequence IDs, and the values are sequences.
    """
    sequences = {}
    try:
        for record in SeqIO.parse(file_path, "fasta"):
            sequences[record.id] = str(record.seq)
    except Exception as e:
        print(f"An error occurred while parsing the FASTA file: {e}")
    return sequences

file_path = "C:\\Users\\Admin\\OneDrive\\Desktop\\p53fasta\\p53_seq.fasta.txt"
parsed_sequences = parse_fasta(file_path)
print(parsed_sequences)
