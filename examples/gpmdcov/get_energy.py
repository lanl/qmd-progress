#!/usr/bin/env python

import argparse
import re

energy_re = re.compile(r"Energy Total.*=\s+([0-9-.]+)")

parser = argparse.ArgumentParser()
parser.add_argument("OUT", help="The output")
options = parser.parse_args()

with open(options.OUT, encoding="utf-8") as fd:
    for line in fd:
        result = energy_re.search(line)
        if result:
            print(result.group(1))
