#!/usr/bin/env python3

# modified from https://raw.githubusercontent.com/J35P312/poorpipe/cd0ccd0638a0b182d62d9d75b3237b6087ebcce0/poorpipe/utils/clean_sniffles.py

import sys

for line in open(sys.argv[1]):
    if line[0] == "#":
        if line.startswith('##FILTER'):
            print('##FILTER=<ID=STRANDBIAS,Description="Not sure what this means">')
        print(line.strip())
        continue

    if "SVTYPE=INS" in line or "SVTYPE=BND" in line or "SVTYPE=DUP/INS" in line:
        var=line.strip().split(";END=")
        before=var[0]
        after=var[-1].lstrip("0123456789")
        line=before+after

    if "SVTYPE=BND" in line:
        var=line.strip().split(";SVLEN=1")
        before=var[0]
        after=var[-1]
        line=before+after
 
    print(line.strip())
