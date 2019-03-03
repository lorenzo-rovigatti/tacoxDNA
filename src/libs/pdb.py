import numpy as np
import itertools
from math import sqrt
import sys

BASE_SHIFT = 1.13
COM_SHIFT = 0.5
FROM_OXDNA_TO_ANGSTROM = 8.518
FROM_ANGSTROM_TO_OXDNA = 1. / FROM_OXDNA_TO_ANGSTROM

NAME_TO_BASE = {
        "ADE" : "A",
        "CYT" : "C",
        "GUA" : "G",
        "THY" : "T",
    }

BASES = ["A", "T", "G", "C", "U"]

class Nucleotide(object):
    RNA_warning_printed = False
    
    def __init__(self, name, idx):
        object.__init__(self)
        self.name = name.strip()
        if self.name in NAME_TO_BASE.keys():
            self.base = NAME_TO_BASE[self.name]
        elif self.name in BASES:
            if self.name == "U" and not Nucleotide.RNA_warning_printed:
                print >> sys.stderr, "WARNING: uracil detected: since RNA is not supported, U will be treated as a T"
                self.name = "T"
                Nucleotide.RNA_warning_printed = True
                
            self.base = self.name
        else:
            self.base = name[1:]
        self.idx = idx
        self.base_atoms = []
        self.phosphate_atoms = []
        self.sugar_atoms = []
        self.named_atoms = {}
        self.ring_names = ["C2", "C4", "C5", "C6", "N1", "N3"]
        self.chain_id = None

    def get_atoms(self):
        return self.base_atoms + self.phosphate_atoms + self.sugar_atoms

    def add_atom(self, a):
        if 'P' in a.name or a.name == "HO5'": self.phosphate_atoms.append(a)
        elif "'" in a.name: self.sugar_atoms.append(a)
        else: self.base_atoms.append(a)
        
        self.named_atoms[a.name] = a
        if self.chain_id == None: self.chain_id = a.chain_id

    def get_com(self, atoms=None):
        if atoms == None: atoms = self.atoms
        com = np.array([0., 0., 0.])
        for a in atoms:
            com += a.pos

        return com / len(atoms)

    def compute_a3(self):
        base_com = self.get_com(self.base_atoms)
        # the O4' oxygen is always (at least for non pathological configurations, as far as I know) oriented towards the 3' direction
        parallel_to = base_com - self.named_atoms["O4'"].pos
        self.a3 = np.array([0., 0., 0.])
        
        for perm in itertools.permutations(self.ring_names, 3):
            p = self.named_atoms[perm[0]]
            q = self.named_atoms[perm[1]]
            r = self.named_atoms[perm[2]]
            v1 = p.pos - q.pos
            v2 = p.pos - r.pos
            v1 /= sqrt(np.dot(v1, v1))
            v2 /= sqrt(np.dot(v2, v2))
            if abs(np.dot(v1, v2)) > 0.01 or 1:
                a3 = np.cross(v1, v2)
                a3 /= sqrt(np.dot(a3, a3))
                if np.dot(a3, parallel_to) < 0.: a3 = -a3
                self.a3 += a3

        self.a3 /= sqrt(np.dot(self.a3, self.a3))

    def compute_a1(self):
        if "C" in self.name or "T" in self.name:
            pairs = [ ["N3", "C6"], ["C2", "N1"], ["C4", "C5"] ]
        else:
            pairs = [ ["N1", "C4"], ["C2", "N3"], ["C6", "C5"] ]

        self.a1 = np.array([0., 0., 0.])
        for pair in pairs:
            p = self.named_atoms[pair[0]]
            q = self.named_atoms[pair[1]]
            diff = p.pos - q.pos
            self.a1 += diff

        self.a1 /= sqrt(np.dot(self.a1, self.a1))

    def compute_as(self):
        self.compute_a1()
        self.compute_a3()
        self.a2 = np.cross(self.a3, self.a1)
        self.check = abs(np.dot(self.a1, self.a3))
        
    def correct_for_large_boxes(self, box):
        map(lambda x: x.shift(-np.rint(x.pos / box ) * box), self.atoms)

    def to_pdb(self, chain_identifier, print_H, serial_residue, residue_suffix):
        res = []
        for a in self.atoms:
            if not print_H and 'H' in a.name:
                continue
            if residue_suffix == "5" and 'P' in a.name:
                continue
            res.append(a.to_pdb(chain_identifier, serial_residue, residue_suffix))

        return "\n".join(res)

    def to_mgl(self):
        res = []
        for a in self.atoms:
            res.append(a.to_mgl())

        return "\n".join(res)

    def rotate(self, R):
        com = self.get_com()
        for a in self.atoms:
            a.pos = np.dot(R, a.pos - com) + com

        self.compute_as()

    def set_com(self, new_com):
        com = self.get_com()
        for a in self.atoms:
            a.pos += new_com - com - COM_SHIFT * self.a1

    def set_base(self, new_base_com):
        atoms = [v for k, v in self.named_atoms.iteritems() if k in self.ring_names]
        ring_com = self.get_com(atoms)
        for a in self.atoms:
            a.pos += new_base_com - ring_com - BASE_SHIFT * self.a1

        self.compute_as()

    atoms = property(get_atoms)


class Atom(object):
    serial_atom = 1

    def __init__(self, pdb_line):
        object.__init__(self)
        # http://cupnet.net/pdb-format/
        self.name = pdb_line[12:16].strip()
        self.alternate = pdb_line[16]
        self.residue = pdb_line[17:20].strip()
        self.chain_id = pdb_line[21:22].strip()
        self.residue_idx = int(pdb_line[22:26])
        self.pos = np.array([float(pdb_line[31:38]), float(pdb_line[38:46]), float(pdb_line[46:54])])
        
    def is_hydrogen(self):
        return "H" in self.name

    def shift(self, diff):
        self.pos += diff

    def to_pdb(self, chain_identifier, serial_residue, residue_suffix):
        residue = self.residue + residue_suffix
        res = "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%-4s%-2s%-2s" % ("ATOM", Atom.serial_atom, self.name, " ", residue, chain_identifier, serial_residue, " ", self.pos[0], self.pos[1], self.pos[2], 1.00, 0.00, " ", " ", " ")
        Atom.serial_atom += 1
        if Atom.serial_atom > 99999:
            Atom.serial_atom = 1
        return res

    def to_mgl(self):
        colors = {"C" : "0,1,1", "P" : "1,1,0", "O" : "1,0,0", "H" : "0.5,0.5,0.5", "N" : "0,0,1"}
        for c in colors:
            if c in self.name: color = colors[c]
        r = 0.5
        return "%s %s %s @ %f C[%s]" % (self.pos[0], self.pos[1], self.pos[2], r, color)
