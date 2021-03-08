#!/usr/bin/env python

def compare_MD(reference, current, reltol):
    """Compare MD energies.

    Given a reference output and the current output, compare the MD
    energies to within the relative tolerance given by reltol.

    """

    import sys

    fd = open(reference)
    reference_energies = []
    for line in fd:
        result = line.split()
        reference_energies.append(float(result[0]))
    fd.close()
    
    fd = open(current)
    current_energies = []
    for line in fd:
        result = line.split()
        current_energies.append(float(result[0]))
    fd.close()

    for i in range(len(reference_energies)):
        diff = abs(reference_energies[i] - current_energies[i])
        if reference_energies[i] != 0:
            diff = abs(diff/reference_energies[i])
        if diff > reltol:
            print("Error")	
            exit(0)
    print("Ok")

def main():
    """The main function.
    """

    import argparse, os, sys
    
    parser = argparse.ArgumentParser(description="""Script to compare MD results by using the total energy""")
    parser.add_argument("--reference",
                        help="The reference output")
    parser.add_argument("--current",
                        help="The current output")
    parser.add_argument("--reltol",
                        help="Relative tolerance when comparing, default is %(default)s",
                        type=float,
                        default=1e-10)
    options = parser.parse_args()

    compare_MD(options.reference, options.current, options.reltol)

if __name__ == "__main__":
    main()
