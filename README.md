# Convert_CIF_2_PDB
 
CIF to PDB Converter

cif_2_pdb.py is a Python script that converts CIF (Crystallographic Information File) files to PDB (Protein Data Bank) format. This script is particularly useful for multichain structures and batch CIF downloads. PDB files are simple tab delimited text files that are easier for various computational analyses.

Prerequisites

Before you begin, ensure you have met the following requirements:

Python 3.x
The following Python libraries:
pandas
biopython
biopandas
Installation

To install the required dependencies, you can use pip:

pip install pandas biopython biopandas

Usage

To use this script, follow these steps:

Clone the repository or download the cif_2_pdb.py script to your local machine.
Open a terminal and navigate to the directory containing the cif_2_pdb.py script.

Run the script using the following command:

python cif_2_pdb.py --folder <folder_with_input_cif_files> --output <folder_with_output_pdb_files> [--segi <segments_in_output>]


Arguments

--folder: Path to the folder containing input CIF files.

--output: Path to the folder where output PDB files will be saved.

--segi: (Optional) Segments in output PDB files. this can be a single keyword ("ubiquitin") or multiple comma separated keywords "mRNA,tRNA"


Examples

1. python cif_2_pdb.py --folder input_cif_files/ --output output_pdb_files/ --segi HAEMAGGLUTININ
2. python cif_2_pdb.py --f ./ --o ./ --segi "16S, mRNA"

   
Additional Information

The script uses the BioPython library to handle the conversion process.

Warnings are suppressed to ensure a clean output, but you can remove the warnings.filterwarnings('ignore') line if you want to see the warnings.


License

This project is licensed under the MIT License - see the LICENSE file for details.

Acknowledgements

BioPython
Pandas
BioPandas
Contact

If you have any questions or suggestions, please open an issue or contact the repository owner.
