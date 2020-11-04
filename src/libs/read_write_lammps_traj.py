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
                line = f.readline()
                N = self.natoms

                #self.xyz = np.zeros((N, 3), dtype=float)
                for n in range(self.natoms):

                    position = np.float32(line.split()[3:6])
                    self.xyz = position
                    velocity = np.float32(line.split()[6:9])
                    self.velocity = velocity

                    print(position,velocity)
                    line = f.readline()
              
                    #print(line)







            line = f.readline()
            #print(line)

          


        sys.exit(0)