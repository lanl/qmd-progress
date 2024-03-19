#!/usr/bin/env python
import sys
import numpy

MyFileName = str(sys.argv[1])

count = -1
MyFile = open(MyFileName, "r", encoding="utf-8")

for lines in MyFile:
    count = count + 1
Dim = count
datos = numpy.zeros(Dim+1)

count = -1
MyFile = open(MyFileName, "r", encoding="utf-8")
for lines in MyFile:
    lines_split = lines.split()
    count = count + 1
    if len(lines_split) > 2:
        if lines_split[0] == "Energy":
            if lines_split[1] == "Total":
                print(" ", lines_split[4])
