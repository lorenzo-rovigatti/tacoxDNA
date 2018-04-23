import numpy as np
import logging
import string
import functools
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
    def __init__(self, file):
	self.filename=file

        head, sects = self.grab_datafile()

	self.natoms=int(head['atoms'])
	self.nbonds=int(head['bonds'])	
	self.nellipsoids=int(head['ellipsoids'])
        x1, x2 = np.float32(head['xlo xhi'].split())
        self.Lx = x2 - x1
        y1, y2 = np.float32(head['ylo yhi'].split())
        self.Ly = y2 - y1
        z1, z2 = np.float32(head['zlo zhi'].split())
        self.Lz = z2 - z1



	self.parse_Atoms_header(sects['Atoms'])

	self.parse_bonds(sects['Bonds'])

	self.parse_ellipsoids(sects['Ellipsoids'])

        self.nstrands=len(np.unique(self.strand))

    def parse_Atoms_header(self,datalines):

	if self.natoms != len(datalines):
		raise ValueError("Number of atoms in header %d and in Atoms %d do not coincide" % self.natoms,len(datalines))
        # Fields per line
	if len(datalines[1].split())!=8:
            raise ValueError("Atoms section should be the default one # Atom-ID, type, position, molecule-ID, ellipsoid flag, density with 8 columns and not %d" % len(datalines[1].split()))
        else :
		N=self.natoms
		# atom ids aren't necessarily sequential
		self.bases = np.zeros(N, dtype=int)
		self.strand = np.zeros(N, dtype=int) 
		self.xyz = np.zeros((N,3), dtype=float) 
		for i, line in enumerate(datalines):
		    line = line.split()
		    index=int(line[0])-1
		    self.bases[index] = line[1]
		    self.strand[index] = line[5]
		    self.xyz[index,:] = line[2:5]

    def parse_bonds(self,datalines):

	if len(datalines[1].split())!=4:
		raise ValueError("Bonds section should have 4 columns and not %d" %len(datalines[1].split()) )

	if self.nbonds !=len(datalines):
		raise ValueError("Number of atoms in header %d and in Bonds %d do not coincide" % self.nbonds,len(datalines))	
	else:
		nbonds=self.nbonds
		self.bonds=np.zeros((nbonds,2),dtype=int)
		for i, line in enumerate(datalines):
			line = line.split()
			index=int(line[0])-1
			self.bonds[index]= line[2:4]


    def parse_ellipsoids(self,datalines):

	if len(datalines[1].split()) != 8:
		raise ValueError("Ellipsoid section should be the default one # Atom-ID, shape, quaternion with 8 columns and not %d" % len(datalines[1].split()))

        if self.nellipsoids !=len(datalines):
                raise ValueError("Number of ellipsoids in header %d and in Bonds %d do not coincide" % self.nellipsoids,len(datalines))
        else:
        	nellipsoids=self.nellipsoids
                self.ellipsoids=np.zeros((nellipsoids,4),dtype=float)
                for i, line in enumerate(datalines):
                        line = line.split()
                        index=int(line[0])-1
                        self.ellipsoids[index,:]= line[4:8]

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

	#we save here the lammps init header information (mass, N, etc)
        header = {}
        for line in f[:starts[0]]:
	    for token in HEADERS:
                if line.endswith(token):
                    header[token] = line.split(token)[0]
		    continue
        #we associate to each section the content (Atoms, Bonds, etc)
	sects = {f[l]:f[l+1:starts[i+1]]
                 for i, l in enumerate(starts[:-1])}
        
	if 'Atoms' not in sects:
            raise ValueError("Data file was missing Atoms section")
        
	return header, sects

