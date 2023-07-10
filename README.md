![Anaconda](https://img.shields.io/badge/Anaconda-47dfd3)![Visual Studio Code](https://img.shields.io/badge/VisualStudioCode-cbf808)![Powershell](https://img.shields.io/badge/Powershell-c2bdc9)![Biopython](https://img.shields.io/badge/Biopython-b470ff)

![Python application test with Github Actions](https://github.com/lakshmanwadhwani/CRISPR-Target-Genie/actions/workflows/main.yml/badge.svg)

![Target Genie Image](./Target%20Genie%20image.png)

# CRISPR-Target-Genie
CRISPR Target Genie 🧬

CRISPR Target Genie is a software tool I designed to aid researchers in identifying potential target sequences in cancer genes for CRISPR/Cas9 editing. The tool accepts genomic sequences of cancer genes as input and outputs potential CRISPR target sites, which can be used for gene knockout or modification. CRISPR Target Genie utilizes BioPython for parsing various genomic data formats, pattern finding, and basic sequence manipulation.

## Features 🚀
- Ability to input a list of known cancer-related genes or their genomic sequences.
- CRISPR target prediction based on gene sequences.
- A scoring system for each target site based on the likelihood of successful binding and editing.
- User-friendly graphical interface for easy input and output handling.
- Export results in a convenient format for further analysis.

## Getting Started 🏁
### Prerequisites
- Python 3.6 or higher
- BioPython
- Tkinter

### Installation
1. Clone this repository:
   https://github.com/yourusername/CRISPR-Target-Genie.git
2. Navigate to the cloned repository:
   cd CRISPR-Target-Genie
3. Install the required dependencies:
   pip install biopython

### Usage
1. Run the script: python Target_genie_crispr_GUI.py
2. Use the graphical interface to select a FASTA file containing gene sequences, enter a PAM sequence, and choose an output file location.
3. Click on `Execute` and wait for the results to be generated.

## Contributing 🤝
Contributions are welcome! For major changes, please open an issue first to discuss what you would like to change.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/YourFeature`)
3. Commit your Changes (`git commit -m 'Add some YourFeature'`)
4. Push to the Branch (`git push origin feature/YourFeature`)
5. Open a Pull Request

## License 📄
Distributed under the MIT License. See `LICENSE` for more information.

## Acknowledgments 🙌
- [BioPython](https://biopython.org/) for providing the library to work with biological sequences.

## Contact 📬
Your Name - lmarkal.wadhwani@gmail.com

   
   
   

   





