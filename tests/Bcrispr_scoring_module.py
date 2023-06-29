from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class CRISPRScoringModule:
    def __init__(self, fasta_file, pam_sequence='CGG'):
        self.fasta_file = fasta_file
        self.pam_sequence = pam_sequence.upper()
        self.target_sites = []

    def find_target_sites(self):
        records = SeqIO.parse(self.fasta_file, "fasta")
        # Implement the sliding window approach to find target sites
        for record in records:
            sequence = str(record.seq).upper()
            for i in range(len(sequence) - len(self.pam_sequence) + 1):
                if sequence[i:i + len(self.pam_sequence)] == self.pam_sequence:
                    self.target_sites.append({
                        'start': i,
                        'end': i + len(self.pam_sequence),
                        'sequence': sequence[i:i + len(self.pam_sequence)],
                        'score': None
                    })

    def score_target_sites(self):
        # Score based on GC content
        for target_site in self.target_sites:
            target_site['score'] = gc_fraction(target_site['sequence'])

    def display_results(self):
        # Display the results
        for target_site in self.target_sites:
            print(f"Start: {target_site['start']}, End: {target_site['end']}, Sequence: {target_site['sequence']}, Score: {target_site['score']:.2f}")


if __name__ == "__main__":
    fasta_file = r'C:\Users\Admin\OneDrive\Desktop\p53fasta'
    crispr_tool = CRISPRScoringModule(fasta_file)
    crispr_tool.find_target_sites()
    crispr_tool.score_target_sites()
    crispr_tool.display_results()
