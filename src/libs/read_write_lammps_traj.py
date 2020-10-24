import numpy as np
import sys

class Lammps_parser(object):	

    def __init__(self, filename):

        self.filename = filename

        f = open(self.filename,'r')
        t = []

        line = f.readline()

        while line:
 
            if line.startswith('ITEM: TIMESTEP'):
                t = int(f.readline())
                line = f.readline()

            if line.startswith('ITEM: NUMBER OF ATOMS') and t==0:
                self.natoms = int(f.readline())
                line = f.readline()

            if line.startswith('ITEM: BOX BOUNDS') and t==0:
                line = f.readline()
                xlo, xhi = np.float32(line.split()[0]), np.float32(line.split()[1])
                self.Lx = xhi - xlo
                line = f.readline()
                ylo, yhi = np.float32(line.split()[0]), np.float32(line.split()[1])
                self.Ly = yhi - ylo
                line = f.readline()
                zlo, zhi = np.float32(line.split()[0]), np.float32(line.split()[1])
                self.Lz = zhi - zlo

            if line.startswith('ITEM: ATOMS'):
                print(line)
                for n in range(self.natoms):
                    line = f.readline()
                    print(t, line)

            line = f.readline()




        sys.exit(0)

