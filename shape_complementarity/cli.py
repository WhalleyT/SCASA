import argparse
import sys

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
    parser.add_argument("--level", "-L", dest="asa_level", type=str, default="R", required=False,
                        help="Level to calculate ASA and BSA to. They can be 'S' for complex, 'C' for chain"
                             " 'R' for residue, or 'A' for atom")
    args = parser.parse_args()

    if args.asa_level not in ["S", "C", "R", "A"]:
        sys.exit("--level/-L must be one of the following: S, C, R, A")
    else:
        return args