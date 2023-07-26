import os
from Bio.PDB import PDBParser
from Bio.PDB import ShrakeRupley

import pandas as pd

class PDBLevelNotFound(Exception):
    """
    Custom exception raised the PDB level is wrong
    """

    def __init__(self, arg):
        self.arg = arg
        super().__init__(self.arg)

class SurfaceArea:
    def __init__(self):
        pass

    def sasa(self, name, file, complex_level="R"):
        """
        This function returns a struct obj from biopython that has residue level ASA
        We will use this to compute BSA by taking ASA(complex_1 alone) - ASA(complex_1 bound)
        :return: struct, biopython struct obj
        """

        if complex_level not in ["S", "C", "R", "A"]:
            raise PDBLevelNotFound(f"Complex level {complex_level} in complex_1 was not found in PDB file {self.pdb_file}")

        p = PDBParser(QUIET=1)
        struct = p.get_structure(name, file)
        sr = ShrakeRupley()
        sr.compute(struct, level=complex_level)
        return struct

    def structure_sasa(self):
        self.complex_sasa = self.sasa(complex_level="S")

    def residue_sasa(self):
        self.complex_sasa = self.sasa(complex_level="R")

    def chain_sasa(self):
        self.complex_sasa = self.sasa(complex_level="C")

    def atom_sasa(self):
        self.complex_sasa = self.sasa(complex_level="A")


    def create_residue_asa_df(self, complex):
        df_list = []

        for chain in complex[0]:
                for res in chain:
                    hetflag, resseq, icode = res.get_id()
                    df_list.append([chain.id, resseq, res.get_resname(), res.sasa])

        return pd.DataFrame(df_list, columns=["Chain", "Residue", "Amino_Acid", "ASA"])


    def complex_sasa(self):
        residue_sasa_complex_1 = self.sasa("complex_1", "tmp/complex_1.pdb")
        residue_sasa_complex_2 = self.sasa("complex_2", "tmp/complex_2.pdb")
        residue_sasa_complex_1_2 = self.sasa("complex_1_2", "tmp/complex_1_2.pdb")

        self.complex_1_asa_df = self.create_residue_asa_df(residue_sasa_complex_1)
        self.complex_2_asa_df = self.create_residue_asa_df(residue_sasa_complex_2)
        self.complex_1_2_asa_df = self.create_residue_asa_df(residue_sasa_complex_1_2)


