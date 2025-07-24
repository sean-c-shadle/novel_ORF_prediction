#!/bin/python3

import sys

workfile = sys.argv[1]

with open(workfile, 'r') as f:
    for line in f:
        line = line.strip()
        fields = line.split("\t")
        
        #extract required fields
        chr_field = fields[7]
        ID_field = fields[0]
        strand = fields[8]
        ranges = fields[-1].split(",")
        
        for r in ranges:
            locations = r.split("-")
            print(f"{chr_field}\t{locations[0]}\t{locations[1]}\t{ID_field}\t0\t{strand}")
