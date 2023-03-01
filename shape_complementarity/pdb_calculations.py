import os
from Bio.PDB import PDBParser
from Bio.PDB import ShrakeRupley

class PDBCalculations:
    def __init__(self):
        pass

    def sasa(self, complex_level="R"):
        """
        This function returns a struct obj from biopython that has residue level ASA
        We will use this to compute BSA by taking ASA(complex_1 alone) - ASA(complex_1 bound)
        :return: struct, biopython struct obj
        """
        p = PDBParser(QUIET=1)
        struct = p.get_structure(self.pdb_name, self.pdb_file)
        sr = ShrakeRupley()
        sr.compute(struct, level=complex_level)
        return struct



