````markdown name=README.md
# CIF to PDB Converter

`cif_2_pdb.py` is a Python script that converts CIF (Crystallographic Information File) files to PDB (Protein Data Bank) format.

This script converts multi chain CIF files to PDB format and allows for segment keyword search.

CODE WILL ONLY WRITE CHAINS WITH SINGLE CHARACTER IDs. ALL OTHER CHAINS WITH 2-LETTER IDs WILL BE IGNORED (NOT INCLUDED) IN THE OUTPUT. This is because the traditional PDB format accepts a single letter chain ID. 

## Prerequisites

Before you begin, ensure you have met the following requirements:

- Python 3.x
- The following Python libraries:
  - pandas
  - biopython
  - biopandas

## Installation

To install the required dependencies, you can use `pip`:

```bash
pip install pandas biopython biopandas
```

## Usage

To use this script, follow these steps:

1. Clone the repository or download the `cif_2_pdb.py` script to your local machine.
2. Open a terminal and navigate to the directory containing the `cif_2_pdb.py` script.
3. Run the script using the following command:

```bash
python cif_2_pdb.py --folder <folder_with_input_cif_files> --output <folder_with_output_pdb_files> [--segi <segments_in_output>]
```

### Arguments

- `--folder`: Path to the folder containing input CIF files.
- `--output`: Path to the folder where output PDB files will be saved.
- `--segi`: (Optional) Segments to include in output PDB files. This is a keyword search, should be 4 characters or less; if the keyword is longer, it will be truncated to first 4 characters (ubiquitin = ubiq). Segi is listed as the last column in the PDB output.

### Examples

```bash
python cif_2_pdb.py --folder input_cif_files/ --output output_pdb_files/ --segi ubiquitin

python cif_2_pdb.py --folder input_cif_files/ --output output_pdb_files/ --segi "16S, 23S"
```

## Additional Information

- The script uses the BioPython library to handle the conversion process.
- Warnings are suppressed to ensure a clean output, but you can remove the `warnings.filterwarnings('ignore')` line if you want to see the warnings.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

- [BioPython](https://biopython.org/)
- [Pandas](https://pandas.pydata.org/)
- [BioPandas](https://github.com/rasbt/biopandas)

## Contact

If you have any questions or suggestions, please open an issue or contact the repository owner.

````
