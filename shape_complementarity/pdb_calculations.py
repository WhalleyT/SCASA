import os
from Bio.PDB import PDBParser
from Bio.PDB import ShrakeRupley

class PDBCalculations:
    def __init__(self):
        pass

    def sasa(self, name, file, complex_level="R"):
        """
        This function returns a struct obj from biopython that has residue level ASA
        We will use this to compute BSA by taking ASA(complex_1 alone) - ASA(complex_1 bound)
        :return: struct, biopython struct obj
        """
        p = PDBParser(QUIET=1)
        struct = p.get_structure(name, file)
        sr = ShrakeRupley()
        sr.compute(struct, level=complex_level)
        return struct

    def complex_sasa(self):
        self.complex_sasa = self.sasa(complex_level="S")

    def residue_sasa(self):
        self.complex_sasa = self.sasa(complex_level="R")

    def chain_sasa(self):
        self.complex_sasa = self.sasa(complex_level="C")

    def atom_sasa(self):
        self.complex_sasa = self.sasa(complex_level="A")

    def complex_sasa(self):
        residue_sasa_complex_1 = self.sasa("complex_1", "tmp/complex_1.pdb")
        residue_sasa_complex_2 = self.sasa("complex_2", "tmp/complex_2.pdb")
        residue_sasa_complex_1_2 = self.sasa("complex_1_2", "tmp/complex_1_2.pdb")

        print(residue_sasa_complex_1)



