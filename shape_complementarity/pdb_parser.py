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

class ChainNotFoundException(Exception):
    '''Custom exception raised when there are no variants detected
    '''
    def __init__(self, arg):
        self.arg = arg
        super().__init__(self.arg)

class Complex():
    """
    Object for defining and holding a PDB complex
    """

    def __init__(self, PDBFile, Complex1, Complex2=None, verbose = False):
        self.pdb_ranges = {"ATOM": range(0, 4),
                           "SERIAL": range(6, 11),
                           "ATOM_NAME": range(12, 16),
                           "ALT_LOC": range(16, 17),
                           "RESIDUE_SEQID": range(17, 20),
                           "CHAIN": range(21, 22),
                           "RESIDUE_NUM": range(22, 26),
                           "INSERTION": range(26, 27),
                           "X_COORD": range(30, 38),
                           "Y_COORD": range(38, 46),
                           "Z_COORD": range(46, 54),
                           "OCCUPANCY": range(54, 60),
                           "TEMPERATURE": range(60, 66),
                           "SEGMENT": range(72, 76),
                           "ELEMENT": range(76, 78)}

        self.pdb_file = PDBFile
        self.complex_1 = Complex1
        self.complex_2 = Complex2
        self.chains = None
        self.verbose = verbose

        self.verify_chains()
        self.chain_string_to_list()

        if self.verbose:
            print("PDB successfully validated")
            print(f"Complex 1: {''.join(self.complex_1)}")
            print(f"Complex 2: {''.join(self.complex_2)}")


    def verify_chains(self):
        if Path(self.pdb_file).is_file() == False:
            raise FileNotFoundError(f"{self.pdb_file} does not exist.")

        self.get_all_chains()


        for chain in self.complex_1:
            if chain not in self.chains:
                raise ChainNotFoundException(f" Chain {chain} in complex_1 was not found in PDB file {self.pdb_file}")

        if self.complex_2 is None:
            self.complex_2 = ""
            for chain in self.chains:
                if chain not in list(self.complex_1):
                    self.complex_2 +=chain
            else:
                for chain in self.complex_2:
                    if chain not in self.chains:
                        raise ChainNotFoundException(f" Chain {chain} in complex_2 was not found in PDB file {self.pdb_file}")

    def get_column(self, key):
        list_of_col = []

        with open(self.pdb_file) as f:
            for line in f:
                if line.startswith("ATOM"):
                    indexes = list(self.pdb_ranges[key])
                    contents = "".join([line[x] for x in indexes]).strip()
                    list_of_col.append(contents)

        return list_of_col

    def get_all_chains(self):
        chains = set(self.get_column("CHAIN"))
        self.chains = sorted(list(set(chains)))
        if self.verbose:
            print(f"PDB file contains chains: {self.chains}")

    def chain_string_to_list(self):
        self.complex_1 = list(self.complex_1)
        self.complex_2 = list(self.complex_2)

