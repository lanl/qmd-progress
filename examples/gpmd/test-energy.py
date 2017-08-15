#!/usr/bin/env python


def read_energies(filename):
    """Read a file and return a list of the read energies."""

    energies = []
    with open(filename) as fd:
        for line in fd:
            result = line.split()
            energies.append(float(result[0]))
    return energies


def dump_file(filename):
    """Print a file to screen."""

    print("reading %s" % filename)
    with open(filename) as fd:
        for line in fd:
            print(line.rstrip())


def compare_MD(reference, current, reltol):
    """Compare MD energies.

    Given a reference output and the current output, compare the MD
    energies to within the relative tolerance given by reltol.

    """

    reference_energies = read_energies(reference)
    current_energies = read_energies(current)

    if len(reference_energies) != len(current_energies):
        dump_file(reference)
        dump_file(current)
        raise Exception(
            ("[error] different number of MD steps\n" +
             "  reference ran for %4d steps\n" +
             "  current ran for   %4d steps\n" +
             "  can not compare") %
            (len(reference_energies), len(current_energies)))

    result = True
    for i in range(len(reference_energies)):
        diff = abs(reference_energies[i] - current_energies[i])
        if reference_energies[i] != 0:
            diff = abs(diff / reference_energies[i])
        if diff > reltol:
            print("failure in MD step %d" % (i + 1))
            result = False
    if not result:
        raise Exception(
            ("[error] when comparing '%s' with '%s'" % (reference, current)) +
            "energies do not agree")

    print("Energy test passed without failure ...")


def main():
    """The main function.
    """

    import argparse

    parser = argparse.ArgumentParser(
        description="Script to compare MD results by using the total " +
        "energy")
    parser.add_argument("--reference",
                        help="The reference output")
    parser.add_argument("--current",
                        help="The current output")
    parser.add_argument("--reltol",
                        help="Relative tolerance when comparing, default " +
                        "is %(default)s",
                        type=float,
                        default=1e-10)
    options = parser.parse_args()

    compare_MD(options.reference, options.current, options.reltol)

if __name__ == "__main__":
    main()
