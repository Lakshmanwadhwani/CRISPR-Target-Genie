import tkinter as tk
from tkinter import filedialog, messagebox
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


class CRISPRGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("CRISPR Scoring Tool")

        # File path input
        self.label_file = tk.Label(self, text="FASTA File:")
        self.label_file.pack()
        self.file_path = tk.StringVar()
        self.entry_file = tk.Entry(self, textvariable=self.file_path, width=40)
        self.entry_file.pack()
        self.button_browse = tk.Button(self, text="Browse", command=self.browse_file)
        self.button_browse.pack()

        # PAM sequence input
        self.label_pam = tk.Label(self, text="PAM Sequence:")
        self.label_pam.pack()
        self.pam_sequence = tk.StringVar(value='CGG')
        self.entry_pam = tk.Entry(self, textvariable=self.pam_sequence, width=20)
        self.entry_pam.pack()

        # Output file
        self.label_output = tk.Label(self, text="Output File:")
        self.label_output.pack()
        self.output_file = tk.StringVar()
        self.entry_output = tk.Entry(self, textvariable=self.output_file, width=40)
        self.entry_output.pack()
        self.button_output = tk.Button(self, text="Set Output", command=self.set_output)
        self.button_output.pack()

        # Execute button
        self.button_execute = tk.Button(self, text="Execute", command=self.execute)
        self.button_execute.pack()

    def browse_file(self):
        file_path = filedialog.askopenfilename()
        self.file_path.set(file_path)

    def set_output(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".txt")
        self.output_file.set(file_path)

    def execute(self):
        fasta_file = self.file_path.get()
        pam_sequence = self.pam_sequence.get()
        output_file = self.output_file.get()

        # Check if the output file has been set
        if not output_file:
            messagebox.showerror("Error", "Please set the output file")
            return

        if not validate_pam_sequence(pam_sequence):
            messagebox.showerror("Error", "Invalid PAM sequence")
            return

        crispr_tool = CRISPRScoringModule(fasta_file, pam_sequence)
        crispr_tool.find_target_sites()
        crispr_tool.score_target_sites()
        crispr_tool.write_results_to_file(output_file)

        messagebox.showinfo("Success", "Results written to output file.")


if __name__ == "__main__":
    app = CRISPRGUI()
    app.mainloop()
