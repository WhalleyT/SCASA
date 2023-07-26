# SCASA: Shape Complementarity and Surface Area

## What is Shape Complementarity?
Shape Complementarity (SC) [Lawrence and Colman, 1993](https://pubmed.ncbi.nlm.nih.gov/8263940/) is a measure of "goodness of fit" between two protein interfaces. 
It ranges from 0 to 1, where 1 is maximum compatibility. It relies on the relative shape of the two interfaces between each other. It has been utilised to explore the
relation ship between antibody/T-cell receptor and antigen.

## What is (buried or available) surface area?
 Buried Surface Area (BSA; reviewed by [Ali et al., 2014](https://pubmed.ncbi.nlm.nih.gov/24678666/)) is a geometric quantity that measures to total surface area (in
Å<sup>2</sup>) of a complex buried within another complex.  Available Surface Area (ASA) is the reciprocal; 
 that is, the total surface area accesible. This calculation of accessibility is generally calcualted in relation
to solvent acessible surface area (SASA). SASA imagines that a ball of a given radius is rolling along the surface ofthe protein of interest.
If there is sufficient geometric space to allow the solvent to pass unimpeded it said to be accessible.

This tool calculates ASA by using the Shrake Ruplley algorithm implemented in Biopython. It then computes BSA
by calculating the ASA of the complex of interest unbound to the other complexes in the structure and comparing the two numbers.

## What is the scope of this tool?
This tools is designed to be used as both a standalone tool and a module to be imported.

## Installation
```pip3 install .```

## Usage

### Standalone Tool
There are two functions ```SCASA sc```, which calculates SC and ```SCASA asa``` which calculates ASA. There are 3 common argumemnts:

```
  --pdb INFILE, -P INFILE
                        PDB file of a complex
  --complex_1 COMPLEX_1, -C1 COMPLEX_1
                        Chains of first complex, corresponding to the chains in the PDB file. If supplying
                        multiple chains they must be a single uninterrupted chain e.g. --complex_1 ABC
  --complex_2 COMPLEX_2, -C2 COMPLEX_2
                        Chains of first complex, corresponding to the chains in the PDB file. If not supplied the
                        default would be all remaining chains
  --level ASA_LEVEL, -L ASA_LEVEL
```

The specific flags for ```sc``` are:
```
  --distance DISTANCE, -D DISTANCE
                        Distance parameter used for generating an interface between the two surfaces. Atoms with
                        no neighbours within this range are excluded
```

The specific flags for ```sc``` are:
```
  --level ASA_LEVEL, -L ASA_LEVEL
                        Level to calculate ASA and BSA to. They can be 'S' for complex, 'C' for chain 'R' for
                        residue, or 'A' for atom
```

So to run sc an example command line argument would be: ```SCASA sc --pdb test/data/1FYT.pdb --complex_1 DE --complex_2 CAB --distance 4.0```. and
an example for asa would be ```SCASA asa --pdb test/data/1FYT.pdb --complex_1 DE --complex_2 CAB --level R```

### Module
SCASA can be imported as a module using ```import scasa```.

## Contact
[Tom Whalley](mailto:whalleyt@cardiff.ac.uk)