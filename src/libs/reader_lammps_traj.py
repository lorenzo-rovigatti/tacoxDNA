import numpy as np
import sys

SECTIONS = set([
    'ITEM'
])

class Lammps_parser(object):	

    def __init__(self, filename):
        self.filename = filename

        sects = self.grab_trajfile()
        print(sects)
        self.parse_Atoms(sects['ITEM'])

        self.nstrands = len(np.unique(self.strand))

        # checking the nucleotides have indexes ordered the same way as bonds and are compatible with strands
        for i in range(self.natoms):
            next_bond = self.bonds[i][1]
            
            if next_bond!=-1:
                #check the consecutive bond is on the same strand
                if self.strand[i]!=self.strand[next_bond]:
                    print("Wrong bond arising between two different strands", file=sys.stderr)
                #check the right bond is an higher index (except from the circular closure)
                if next_bond==i-1:
                    print("The bonds should be in incremental order (i i+1) except for strand circularization N 0", file=sys.stderr)
                if next_bond>i+1:
                    print("The bonds should be in incremental order (i i+1) except for strand circularization N 0", file=sys.stderr)
                if next_bond<i+1:
                    if self.bonds[next_bond][1]!=next_bond+1:
                        print("The bonds should be in incremental order (i i+1) except for strand circularization N 0", file=sys.stderr)

                # more check to insert about completely random ordering

    def parse_Atoms(self, datalines):
        if self.natoms != len(datalines):
            raise ValueError("Number of atoms in header %d and in Atoms %d do not coincide" % (self.natoms, len(datalines)))
        # Fields per line
        if len(datalines[1].split()) != 8:
            raise ValueError("Atoms section should be the default one # Atom-ID, type, position, molecule-ID, ellipsoid flag, density with 8 columns and not %d" % len(datalines[1].split()))
        N = self.natoms
        # atom ids aren't necessarily sequential
        self.bases = np.zeros(N, dtype=int)
        self.strand = np.zeros(N, dtype=int) 
        self.xyz = np.zeros((N, 3), dtype=float) 
        for i, line in enumerate(datalines):
            line = line.split()
            index = int(line[0]) - 1
            self.bases[index] = line[1]
            self.strand[index] = line[5]
            self.xyz[index, :] = line[2:5]

    def iterdata(self):
        with open(self.filename) as f:
            for line in f:
                if line:
                    yield line

    def grab_trajfile(self):
        f = list(self.iterdata())
        print(f)
        starts = [i for i, line in enumerate(f)
                  if line.split()[0] in SECTIONS]
        starts += [None]

        # we associate to each section the content (Atoms, Bonds, etc)
        sects = {f[l]:f[l + 1:starts[i + 1]] for i, l in enumerate(starts[:-1])}
        
        return sects

