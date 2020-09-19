import numpy as np
import sys

SECTIONS = set([
    'Atoms',
    'Velocities',
    'Ellipsoids',
    'Bonds',
])

HEADERS = set([
    'atoms',
    'bonds',
    'atom types',
    'bond types',
    'ellipsoids',
    'xlo xhi',
    'ylo yhi',
    'zlo zhi',
])


class Lammps_parser(object):	

    def __init__(self, filename):
        self.filename = filename

        head, sects = self.grab_datafile()

        self.natoms = int(head['atoms'])
        self.nbonds = int(head['bonds'])	
        self.nellipsoids = int(head['ellipsoids'])
        x1, x2 = np.float32(head['xlo xhi'].split())
        self.Lx = x2 - x1
        y1, y2 = np.float32(head['ylo yhi'].split())
        self.Ly = y2 - y1
        z1, z2 = np.float32(head['zlo zhi'].split())
        self.Lz = z2 - z1

        self.parse_Atoms_header(sects['Atoms'])
        self.parse_bonds(sects['Bonds'])
        self.parse_ellipsoids(sects['Ellipsoids'])

        self.parse_velocities(sects['Velocities'])

        self.nstrands = len(np.unique(self.strand))

        # checking the nucleotides have indexes ordered the same way as bonds and are compatible with strands
        for i in range(self.natoms):
            next_bond = self.bonds[i][1]
            
            if next_bond != -1:
                # check the consecutive bond is on the same strand
                if self.strand[i] != self.strand[next_bond]:
                    print >> sys.stderr, "Wrong bond arising between two different strands"
                # check the right bond is an higher index (except from the circular closure)
                if next_bond == i - 1:
                    print >> sys.stderr, "The bonds should be in incremental order (i i+1) except for strand circularization N 0"
                if next_bond > i + 1:
                    print >> sys.stderr, "The bonds should be in incremental order (i i+1) except for strand circularization N 0"
                if next_bond < i + 1:
                    if self.bonds[next_bond][1] != next_bond + 1:
                        print >> sys.stderr, "The bonds should be in incremental order (i i+1) except for strand circularization N 0"

                # more check to insert about completely random ordering

    def parse_Atoms_header(self, datalines):
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

    def parse_velocities(self, datalines):
        if self.natoms != len(datalines):
            raise ValueError("Velocities section do not contain the same amount of particles as number of atoms %d != %d " % self.natoms, len(datalines))
        if len(datalines[1].split()) != 7:
            raise ValueError("Velocities section should be the default one # Atom-ID, translational, rotational velocity with 7 columns and not %d" % len(datalines[1].split()))
        else:
            N = self.natoms
            self.v = np.zeros((N, 3), dtype=float)
            self.Lv = np.zeros((N, 3), dtype=float)
            for i, line in enumerate(datalines):
                line = line.split()
                index = int(line[0]) - 1
                self.v[index] = line[1:4]
                self.Lv[index] = line[4:7]

    def parse_bonds(self, datalines):
        if len(datalines[1].split()) != 4:
            raise ValueError("Bonds section should have 4 columns and not %d" % len(datalines[1].split()))
    
        if self.nbonds != len(datalines):
            raise ValueError("Number of atoms in header %d and in Bonds %d do not coincide" % self.nbonds, len(datalines))	

        # creating a vector indicating for each particle who it is bonded too on its left and right in order of increasing index
        natoms = self.natoms
        self.bonds = np.ones((natoms, 2), dtype=int) * (-1)
        for i, line in enumerate(datalines):
            line = line.split()
            p1 = int(line[2]) - 1
            p2 = int(line[3]) - 1

            self.bonds[p1][1] = p2
            self.bonds[p2][0] = p1

    def parse_ellipsoids(self, datalines):
        if len(datalines[1].split()) != 8:
            raise ValueError("Ellipsoid section should be the default one # Atom-ID, shape, quaternion with 8 columns and not %d" % len(datalines[1].split()))

        if self.nellipsoids != len(datalines):
            raise ValueError("Number of ellipsoids in header %d and in Bonds %d do not coincide" % self.nellipsoids, len(datalines))

        nellipsoids = self.nellipsoids
        self.ellipsoids = np.zeros((nellipsoids, 4), dtype=float)
        for i, line in enumerate(datalines):
            line = line.split()
            index = int(line[0]) - 1
            self.ellipsoids[index, :] = line[4:8]

    def iterdata(self):
        with open(self.filename) as f:
            for line in f:
                line = line.partition('#')[0].strip()
                if line:
                    yield line

    def grab_datafile(self):
        f = list(self.iterdata())
        starts = [i for i, line in enumerate(f)
                  if line.split()[0] in SECTIONS]
        starts += [None]

        # we save here the lammps init header information (mass, N, etc)
        header = {}
        for line in f[:starts[0]]:
            for token in HEADERS:
                if line.endswith(token):
                    header[token] = line.split(token)[0]

        # we associate to each section the content (Atoms, Bonds, etc)
        sects = {f[l]:f[l + 1:starts[i + 1]] for i, l in enumerate(starts[:-1])}
        
        if 'Atoms' not in sects:
            raise ValueError("Data file was missing Atoms section")
        
        return header, sects

