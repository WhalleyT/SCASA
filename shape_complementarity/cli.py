import argparse
import sys

def parse_args():
    """
    Parse the command line arguments
    """
    parser = argparse.ArgumentParser("SCASA: Surface Complementarity and Available Surface Area Calculation")

    subparsers = parser.add_subparsers(dest="command")

    #subparser for asa
    asa = subparsers.add_parser("asa", help="Calcuate Avaible and Buried Surface Area (ASA/BSA) of a PDB complex")
    asa.add_argument('--pdb', '-P', dest='infile', type=str,
                        help='PDB file of a complex', required=True)
    asa.add_argument('--complex_1', '-C1', dest='complex_1', type=str,
                        help='Chains of first complex, corresponding to the chains in the PDB file. If supplying '
                             'multiple chains they must be a single uninterrupted chain e.g. --complex_1 ABC',
                        required=True)
    asa.add_argument('--complex_2', '-C2', dest='complex_2', type=str,
                        help='Chains of first complex, corresponding to the chains in the PDB file. If not supplied '
                             'the default would be all remaining chains',
                        required=False)
    asa.add_argument("--level", "-L", dest="asa_level", type=str, default="R", required=False,
                        help="Level to calculate ASA and BSA to. They can be 'S' for complex, 'C' for chain"
                             " 'R' for residue, or 'A' for atom")

    #subparser for sc
    sc = subparsers.add_parser("sc", help="Calcuate Shape Complementarity (SC) of a PDB complex")
    sc.add_argument('--pdb', '-P', dest='infile', type=str,
                        help='PDB file of a complex', required=True)
    sc.add_argument('--complex_1', '-C1', dest='complex_1', type=str,
                        help='Chains of first complex, corresponding to the chains in the PDB file. If supplying '
                             'multiple chains they must be a single uninterrupted chain e.g. --complex_1 ABC',
                        required=True)
    sc.add_argument('--complex_2', '-C2', dest='complex_2', type=str,
                        help='Chains of first complex, corresponding to the chains in the PDB file. If not supplied '
                             'the default would be all remaining chains',
                        required=False)
    sc.add_argument("--distance", "-D", dest="distance", type=float, default=4.0, required=False,
                        help="Distance parameter used for generating an interface between the two surfaces. Atoms with"
                             " no neighbours within this range are excluded")

    args = parser.parse_args()

    print(args)
    if args.command == "asa":
        if args.asa_level not in ["S", "C", "R", "A"]:
            sys.exit("--level/-L must be one of the following: S, C, R, A")
        else:
            return args
    elif args.command == "sc":
        if args.distance <= 0:
            sys.exit("--distance/-D must be non-negative")
        else:
            return args
    else:
        parser.print_usage()
        parser.print_help()
        sys.exit()

