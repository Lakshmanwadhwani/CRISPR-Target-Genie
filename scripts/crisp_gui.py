import tkinter as tk
from tkinter import filedialog, messagebox
from crispr_scoring_module import CRISPRScoringModule, validate_pam_sequence


class CRISPRGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        # ... rest of the code for GUI ...


if __name__ == "__main__":
    app = CRISPRGUI()
    app.mainloop()
