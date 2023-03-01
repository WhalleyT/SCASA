import os
from Bio.PDB import PDBParser
from Bio.PDB import ShrakeRupley

class PDBCalculations:
    def __init__(self):
        pass

    def sasa(self):
        p = PDBParser()
        struct = p.get_structure(self.pdb_name, self.pdb_file)
        sr = ShrakeRupley()
        sr.compute(struct, level="S")



