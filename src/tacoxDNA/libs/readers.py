from . import base
import numpy as np
import os.path

class LorenzoReader:
    def __init__(self, topology, configuration):
        self._conf = False

        if not os.path.isfile(configuration):
            base.Logger.die("Configuration file '%s' is not readable" % configuration)

        if not os.path.isfile(topology):
            base.Logger.die("Topology file '%s' is not readable" % topology)

        self._conf = open(configuration, "r")

        f = open(topology, "r")
        self.N, self.N_strands = [int(x) for x in f.readline().split()]
        self._top_lines = f.readlines()
        if len(self._top_lines) != self.N:
            raise Exception("The number of nucleotides specified in the topology file header (%d) is different from the number of nucleotide lines found in the same file (%d)" % (self.N, len(self._top_lines)))

    def __del__(self):
        if self._conf: self._conf.close()

    def _read(self, only_strand_ends=False, skip=False):
        try:
            timeline = self._conf.readline()
            time = float(timeline.split()[2])
    
            box = np.array([float(x) for x in self._conf.readline().split()[2:]])
            self._conf.readline()
        except Exception as e:
            raise Exception("The header lines of the configuration file are invalid (caught a '%s' exception)" % e)

        if skip:
            for tl in self._top_lines:
                self._conf.readline()

            return False

        system = base.System(box, time=time)
        base.Nucleotide.index = 0
        base.Strand.index = 0

        s = False
        strandid_current = 0
        for i_line, tl in enumerate(self._top_lines):
            tls = tl.split()
            n3 = int(tls[2])
            n5 = int(tls[3])
            strandid = int(tls[0])
            if (len (tls[1]) == 1):
                b = base.base_to_number[tls[1]]
                bb = b
            else:
                try:
                    tmp = int (tls[1])
                except:
                    raise Exception("The line n. %d in the topology file contains an incorrect specific base pairing" % i_line)

                if tmp > 0:
                    b = tmp % 4
                else:
                    b = (3 - ((3 - tmp) % 4))
                bb = tmp

            if strandid != strandid_current:
                # check for circular strand
                if n3 != -1:
                    iscircular = True
                else:
                    iscircular = False

                if s:
                    system.add_strand(s)
                s = base.Strand()
                if iscircular:
                    s.make_circular()
                strandid_current = strandid

            ls = self._conf.readline().split()
            if len(ls) == 0:
                raise Exception("The %d-th nucleotide line in the configuration file is empty" % i_line)
            elif len(ls) != 15:
                raise Exception("The %d-th nucleotide line in the configuration file is invalid" % i_line)
            cm = [float(x) for x in ls[0:3]]
            a1 = [float(x) for x in ls[3:6]]
            a3 = [float(x) for x in ls[6:9]]
            v = [float(x) for x in ls[9:12]]
            L = [float(x) for x in ls[12:15]]
            if not only_strand_ends or n3 == -1 or n5 == -1:
                s.add_nucleotide(base.Nucleotide(cm, a1, a3, b, bb, v, L, n3))

        system.add_strand(s)
        
        return system


    # if only_strand_ends == True then a strand will contain only the first and the last nucleotide
    # useful for some analysis like csd for coaxial interactions
    def get_system(self, only_strand_ends=False, N_skip=0):
        for _ in range(N_skip):
            self._read(skip=True)

        return self._read(only_strand_ends=only_strand_ends, skip=False)
