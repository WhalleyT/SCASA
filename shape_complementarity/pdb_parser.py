import argparse
from pathlib import Path

def parse_args():
    """
    Parse the command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb', '-P', dest='infile', type=str,
                        help='PDB file of a complex', required=True)
    parser.add_argument('--complex_1', '-C1', dest='complex_1', type=str,
                        help='Chains of first complex, corresponding to the chains in the PDB file. If supplying '
                             'multiple chains they must be a single uninterrupted chain e.g. --complex_1 ABC',
                        required=True)
    parser.add_argument('--complex_2', '-C2', dest='complex_2', type=str,
                        help='Chains of first complex, corresponding to the chains in the PDB file. If not supplied ' 
                        'the default would be all remaining chains',
                        required=False)
    args = parser.parse_args()
    return args


class Complex():
    """
    Object for defining and holding a PDB complex
    """

    def __init__(self, PDBFile, Complex1, Complex2=None):
        self.pdb_file = PDBFile
        self.complex_1 = Complex1
        self.complex_2 = Complex2

        self.verify_chains()


    def verify_chains(self):
        if Path(self.pdb_file).is_file() == False:
            raise FileNotFoundError(f"{self.pdb_file} does not exist.")

        