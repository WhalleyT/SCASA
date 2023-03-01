import shape_complementarity.pdb_parser


def main():
    args =  shape_complementarity.pdb_parser.parse_args()

    #set up class for complex
    complex = shape_complementarity.pdb_parser.Complex(args.infile, args.complex_1, args.complex_2, verbose=True)
    complex.sasa()


if __name__ == "__main__":
    main()