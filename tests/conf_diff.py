#!/usr/bin/env python3

import numpy as np
import sys

EPS = 1e-6

def line_differ(line1, line2):
    try:
        a1 = np.array([float(x) for x in line1.split()])
        a2 = np.array([float(x) for x in line2.split()])
    except ValueError:
        return line1 != line2
    return np.any(np.abs(a1 - a2) > EPS)

if len(sys.argv) < 3:
    print("Usage is %s conf1 conf2" % sys.argv[0])
    exit(1)
    
try:
    with open(sys.argv[1]) as inp1, open(sys.argv[2]) as inp2:
        # the first three lines are headers
        for _ in range(3):
            l1 = inp1.readline()
            l2 = inp2.readline()
            if l1 != l2:
                print("Found a difference in the header line ('%s' vs '%s')" % (l1.strip(), l2.strip()), file=sys.stderr)
                exit(1)
                
        lines_1 = inp1.readlines()
        lines_2 = inp2.readlines()
        
        if len(lines_1) != len(lines_2):
            print("The two files have different numbers of lines (%d vs %d)" % (len(lines_1), len(lines_2)), file=sys.stderr)
            exit(1)
            
        for i, pair in enumerate(zip(lines_1, lines_2)):
            if line_differ(pair[0], pair[1]):
                print("Found difference on line %d:" % (i + 3), file=sys.stderr)
                print("> %s\n< %s\n" % (pair[0].strip(), pair[1].strip()), file=sys.stderr)
            
except FileNotFoundError as e:
    print(e)
    exit(1)
    
exit(0)
