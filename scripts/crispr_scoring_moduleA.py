from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class CRISPRScoringModule:
    def __init__(self, fasta_file, pam_sequence='CGG'):
        self.fasta_file = fasta_file
        self.pam_sequence = pam_sequence.upper()
        self.target_sites = []

    def find_target_sites(self):
        records = SeqIO.parse(self.fasta_file, "fasta")
        for record in records:
            sequence = str(record.seq).upper()
            for i in range(len(sequence) - len(self.pam_sequence) + 1):
                if sequence[i:i + len(self.pam_sequence)] == self.pam_sequence and 'N' not in self.pam_sequence:
                    self.target_sites.append({
                        'start': i,
                        'end': i + len(self.pam_sequence),
                        'sequence': sequence[i:i + len(self.pam_sequence)],
                        'score': None
                    })

    def score_target_sites(self):
        for target_site in self.target_sites:
            target_site['score'] = gc_fraction(target_site['sequence'])

    def write_results_to_file(self, output_filename):
        with open(output_filename, 'w', newline='\n') as output_file:
            output_file.write("Start\tEnd\tSequence\tScore\n")
            for target_site in self.target_sites:
                output_file.write(f"{target_site['start']}\t{target_site['end']}\t{target_site['sequence']}\t{target_site['score']:.2f}\n")


def validate_pam_sequence(pam):
    for base in pam.upper():
        if base not in ['A', 'T', 'C', 'G', 'N']:
            return False
    return 2 <= len(pam) <= 6
