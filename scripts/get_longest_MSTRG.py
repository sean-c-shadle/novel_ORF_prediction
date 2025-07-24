#!/bin/python3

import sys
import re

workfile = sys.argv[1]

data = {}
ID = {}

with open(workfile, 'r') as f:
    for line in f:
        line = line.strip()
        fields = line.split("\t") 
        gene_type = fields[0]
        
        match = re.match(r'^(MSTRG\.\d+\.\d+)(_+\d+_+\d+_+)(\d+)', gene_type)

        if match:
            MSTRG = match.group(1)
            location = match.group(2)
            number = int(match.group(3))

            if MSTRG in data:
                if data[MSTRG] < number:
                    data[MSTRG] = number
                    ID[MSTRG] = location
            else:
                data[MSTRG] = number
                ID[MSTRG] = location

for MSTRG in data:
    print(f"{MSTRG}{ID[MSTRG]}{data[MSTRG]}")
