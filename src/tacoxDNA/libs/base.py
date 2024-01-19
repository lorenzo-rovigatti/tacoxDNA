"""
Utility functions.
base.py includes the classes: System, Strand, Nucleotide
    - Make initial configurations (generate.py)
    - Get detailed energy information (process_data/)
    - If you want to use it with oxRNA, you have to set environment variable OXRNA to 1  (export OXRNA=1) 
"""
import sys
import numpy as np
import os


def partition(s, d):
    if d in s:
        sp = s.split(d, 1)
        return sp[0], d, sp[1]
    else:
        return s, "", ""


number_to_base = {0: 'A', 1: 'G', 2: 'C', 3: 'T'}

base_to_number = {'A': 0, 'a': 0, 'G': 1, 'g': 1,
                  'C': 2, 'c': 2, 'T': 3, 't': 3,
                  'U': 3, 'u': 3, 'D': 4}

try:
    FLT_EPSILON = np.finfo(np.float).eps
except:
    FLT_EPSILON = 2.2204460492503131e-16
    
# oxDNA constants
POS_BACK = -0.4
POS_STACK = 0.34
POS_BASE = 0.4
CM_CENTER_DS = POS_BASE + 0.2
FENE_R0_OXDNA = 0.7525
FENE_EPS = 2.0

POS_MM_BACK1 = -0.3400
POS_MM_BACK2 = 0.3408

# may come in handy later
RNA = False
RNA_POS_BACK_a1 = -0.4
RNA_POS_BACK_a3 = 0.2
RNA_POS_BACK_a2 = 0.0

LENGTH_FACT = 8.518
BASE_BASE = 0.3897628551303122

LENGTH_FACT = 8.518
BASE_BASE = 0.3897628551303122

CREPY_COLOR_TABLE = ['red', 'blue', '0,0.502,0', '1,0.8,0', '0.2,0.8,1']

# INT_POT_TOTAL = 0
INT_HYDR = 4
INT_STACK = 2
INT_CROSS_STACK = 5
INT_COAX_STACK = 6
INT_FENE = 0
INT_EXC_BONDED = 1
INT_EXC_NONBONDED = 3

H_CUTOFF = -0.1

NOGROOVE_ENV_VAR = "OXDNA_NOGROOVE"
MM_GROOVING = os.environ.get(NOGROOVE_ENV_VAR) == '1'

# static class
class Logger(object):
    DEBUG = 0
    INFO = 1
    WARNING = 2
    CRITICAL = 3
    debug_level = INFO

    messages = ("DEBUG", "INFO", "WARNING", "CRITICAL")

    @staticmethod
    def log(msg, level=None, additional=None):
        if level == None: 
            level = Logger.INFO
        if level < Logger.debug_level: 
            return

        if additional != None and Logger.debug_level == Logger.DEBUG:
            print("%s: %s (additional info: '%s')" % (Logger.messages[level], msg, additional), file=sys.stderr)
        else: 
            print("%s: %s" % (Logger.messages[level], msg), file=sys.stderr)

    @staticmethod
    def die(msg):
        Logger.log(msg, Logger.CRITICAL)
        sys.exit()


class Nucleotide():
    """
    Nucleotides compose Strands

    cm_pos --- Center of mass position
        Ex: [0, 0, 0]

    a1 --- Unit vector indicating orientation of backbone with respect to base
        Ex: [1, 0, 0]

    a3 --- Unit vector indicating orientation (tilting )of base with respect to backbone
        Ex: [0, 0, 1]

    base --- Identity of base, which must be designated with either numbers or
        letters (this is called type in the c++ code). Confusingly enough, this
        is similar to Particle.btype in oxDNA.
        
        Number: {0,1,2,3} or any number in between (-inf,-7) and (10, inf)
        To use specific sequences (or an alphabet large than four) one should
        start from the complementary pair 10 and -7. Complementary pairs are
        such that base_1 + base_2 = 3;
        
        Letter: {A,G,T,C} (these will be translated to {0, 1, 3, 2}).
        
        These are set in the dictionaries: number_to_base, base_to_number
    
    btype--- Identity of base. Unused at the moment.

    pair --- Base-paired Nucleotide, used in oxView output
    cluster --- Cluster ID, Number, used in oxView output
    color --- Custom color, Number representing hex value, used in oxView output

    """
    index = 0

    def __init__(self, cm_pos, a1, a3, base, btype=None, v=np.array([0., 0., 0.]), L=np.array([0., 0., 0.]), n3=-1, pair=None, cluster=None, color=None):
        self.index = Nucleotide.index
        Nucleotide.index += 1
        self.cm_pos = np.array(cm_pos)
        self._a1 = np.array(a1)
        self._a3 = np.array(a3)
        self._original_base = base
        # base should be an integer
        if isinstance(base, int) or isinstance(base, np.int_):
            pass
        else:
            try:
                base = base_to_number[base]
            except KeyError:
                Logger.log("Invalid base (%s)" % base, level=Logger.WARNING)
        self._base = base
        if btype is None:
            self._btype = base
        else:
            self._btype = btype
        self._L = L
        self._v = v
        self.n3 = n3
        self.next = -1
        self.pair = pair
        self.cluster = cluster
        self.color = color

    def get_pos_base (self):
        """
        Returns the position of the base centroid
        Note that cm_pos is the centrod of the backbone and base.
        """
        return self.cm_pos + self._a1 * POS_BASE

    pos_base = property(get_pos_base)

    def get_pos_stack (self):
        return self.cm_pos + self._a1 * POS_STACK

    pos_stack = property (get_pos_stack)

    def get_pos_back (self):
        """
        Returns the position of the backbone centroid
        Note that cm_pos is the centrod of the backbone and base.
        """
        if MM_GROOVING:
            return self.cm_pos + self._a1 * POS_MM_BACK1 + self._a2 * POS_MM_BACK2
        elif RNA:
            return self.cm_pos + self._a1 * RNA_POS_BACK_a1 + self._a2 * RNA_POS_BACK_a2 + self._a3 * RNA_POS_BACK_a3
        else:
            return self.cm_pos + self._a1 * POS_BACK

    pos_back = property(get_pos_back)

    def get_pos_back_rel(self):
        """
        Returns the position of the backbone centroid relative to the centre of mass
        i.e. it will be a vector pointing from the c.o.m. to the backbone
        """
        return self.get_pos_back() - self.cm_pos

    def get_a2(self):
        return np.cross (self._a3, self._a1)

    _a2 = property (get_a2)

    def copy(self, disp=None, rot=None):
        copy = Nucleotide(self.cm_pos, self._a1, self._a3, self._base, self._btype, self._L, self._v, self.n3, self.pair, self.cluster, self.color)
        if disp is not None:
            copy.translate(disp)
        if rot is not None:
            copy.rotate(rot)

        return copy

    def translate(self, disp):
        self.cm_pos += disp
        self.cm_pos_box += disp

    def rotate(self, R, origin=None):
        if origin == None: origin = self.cm_pos

        self.cm_pos = np.dot(R, self.cm_pos - origin) + origin
        self._a1 = np.dot(R, self._a1)
        self._a3 = np.dot(R, self._a3)

    def distance(self, other, PBC=True, box=None):
        if PBC and box is None:
            if not (isinstance (box, np.ndarray) and len(box) == 3):
                Logger.die ("distance between nucleotides: if PBC is True, box must be a numpy array of length 3");
        dr = other.pos_back - self.pos_back
        if PBC:
            dr -= box * np.rint (dr / box)
        return dr
    
    def is_uracil(self):
        return isinstance(self._original_base, str) and self._original_base.lower() == "u"

    def get_base(self):
        """
        Returns a number containing base id
        >>> v1 = np.array([0.,0.,0.])
        >>> v2 = np.array([1.,0.,0.])
        >>> v3 = np.array([0.,0.,1.])
        >>> Nucleotide(v1, v2, v3, 'A').get_base()
        'A'
        >>> Nucleotide(v1, v2, v3, "C").get_base()
        'C'
        >>> Nucleotide(v1, v2, v3, "G").get_base()
        'G'
        >>> Nucleotide(v1, v2, v3, "T").get_base()
        'T'
        >>> Nucleotide(v1, v2, v3, 1).get_base()
        'G'
        >>> Nucleotide(v1, v2, v3, 103).get_base()
        '103'
        >>> Nucleotide(v1, v2, v3, -97).get_base()
        '-97'
        """
        if self.is_uracil():
            return "U"
        
        if type(self._base) is not int:
            try:
                number_to_base[self._base]
            except KeyError:
                Logger.log("Nucleotide.get_base(): nucleotide %d: unknown base type '%d', defaulting to 12 (A)" % (self.index, self._base), Logger.WARNING)
        
        if self._base in [0, 1, 2, 3]:
            return number_to_base[self._base]
        else:
            return str(self._base)

    def _get_lorenzo_output(self):
        a = np.concatenate((self.cm_pos, self._a1, self._a3, self._v, self._L))
        return " ".join(str(x) for x in a)

    def _get_ribbon_output(self):
        s1 = self.cm_pos_box + self.get_pos_back_rel()
        return "%lf %lf %lf %lf %lf %lf %lf %lf %lf" % (tuple(s1) + tuple (self._a1) + tuple (self._a2))


class Strand():
    """
    Strand composed of Nucleotides
    Strands can be contained in System
    """
    index = 0

    def __init__(self):
        self.index = Strand.index
        Strand.index += 1
        self._first = -1
        self._last = -1
        self._nucleotides = []
        self._sequence = []
        self.visible = True
        self._circular = False  # bool on circular DNA

    def get_length(self):
        return len(self._nucleotides)

    def get_sequence(self):
        return self._sequence

    def _prepare(self, si, ni):
        self.index = si
        self._first = ni

        for n in range(self.N):
            self._nucleotides[n].index = ni + n

        self._last = ni + n
        return ni + n + 1

    def copy(self):
        copy = Strand()
        for n in self._nucleotides:
            copy.add_nucleotide(n.copy())
        return copy

    def get_cm_pos(self):
        return sum([x.cm_pos for x in self._nucleotides]) / self.N

    def set_cm_pos(self, new_pos):
        diff = new_pos - self.cm_pos
        for n in self._nucleotides: n.translate(diff)

    def translate(self, amount):
        new_pos = self.cm_pos + amount
        self.set_cm_pos(new_pos)

    def rotate(self, R, origin=None):
        if origin == None: 
            origin = self.cm_pos

        for n in self._nucleotides: 
            n.rotate(R, origin)

    def append(self, other):
        if not isinstance (other, Strand):
            raise ValueError

        dr = self._nucleotides[-1].distance (other._nucleotides[0], PBC=False)
        if np.sqrt(np.dot (dr, dr)) > (0.7525 + 0.25):
            print("WARNING: Strand.append(): strands seem too far apart. Assuming you know what you are doing.", file=sys.stderr)

        ret = Strand()

        for n in self._nucleotides:
            ret.add_nucleotide(n)

        for n in other._nucleotides:
            ret.add_nucleotide(n)

        return ret

    def get_slice(self, start=0, end=None):
        if end is None: 
            end = self.get_length()
            
        if end > self.get_length():
            print("The given end parameter is larger than the number of nucleotides of the strand (%d > %d)" % (end, self.get_length()), file=sys.stderr)
            raise ValueError
            
        ret = Strand()
        for i in range(start, end):
            ret.add_nucleotide(self._nucleotides[i].copy())
        return ret

    def set_sequence(self, seq):
        if isinstance(seq, str):
            seq = [base_to_number[x] for x in seq]
        if len(seq) != len(self._nucleotides):
            Logger.log ("Cannot change sequence: lengths don't match", Logger.WARNING)
            return
        i = 0
        for n in self._nucleotides:
            n._base = seq[i]
            i += 1
        self._sequence = seq

    def bring_in_box_nucleotides(self, box):
        diff = np.rint(self.cm_pos / box) * box
        for n in self._nucleotides:
            n.cm_pos_box = n.cm_pos - diff

    def add_nucleotide(self, n):
        if len(self._nucleotides) == 0:
            self._first = n.index
        n.strand = self.index
        self._nucleotides.append(n)
        self._last = n.index
        self.sequence.append(n._base)

    def _get_lorenzo_output(self):
        if not self.visible:
            return ""

        conf = "\n".join(n._get_lorenzo_output() for n in self._nucleotides) + "\n"

        top = ""
        for n in self._nucleotides:
            if self._circular:
                if n.index == self._first:
                    n3 = self._last
                else:
                    n3 = n.index - 1
                if n.index == self._last:
                    n5 = self._first
                else:
                    n5 = n.index + 1
            else:
                if n.index == self._first:
                    n3 = -1
                else:
                    n3 = n.index - 1
                if n.index == self._last:
                    n5 = -1
                else:
                    n5 = n.index + 1
            top += "%d %s %d %d\n" % (self.index + 1, n.get_base(), n3, n5)

        return conf, top

    def get_lammps_N_of_bonds_strand(self):
        N_bonds = 0
        for n in self._nucleotides:
            if n.index != self._last:
                N_bonds += 1
            elif self._circular:
                N_bonds += 1

        return N_bonds

    def get_lammps_bonds(self):
        top = []
        for n in self._nucleotides:
            if n.index != self._last:
                top.append("%d  %d" % (n.index + 1, n.index + 2))
            elif self._circular:
                top.append("%d  %d" % (n.index + 1, self._first + 1))

        return top

    def _get_ribbon_output(self):
        if not self.visible:
            return ""
        v = [n._get_ribbon_output() for n in self._nucleotides]
        v.insert(0, "0. 0. 0. @ 0.3 C[red] RIB")
        return " ".join(v)

    cm_pos = property(get_cm_pos, set_cm_pos)
    N = property(get_length)
    sequence = property(get_sequence)
    
    def make_circular(self, check_join_len=False):
        if check_join_len:
            dr = self._nucleotides[-1].distance (self._nucleotides[0], PBC=False)
            if np.sqrt(np.dot (dr, dr)) > (0.7525 + 0.25):
                Logger.log("Strand.make_circular(): ends of the strand seem too far apart. \
                            Assuming you know what you are doing.", level=Logger.WARNING)
        self._circular = True

    def make_noncircular(self):
        self._circular = False
        
    def is_circular(self):
        return self._circular

    def cut_in_two(self, copy=True):  # cuts a strand into two strands in the middle
        fragment_one = Strand()
        fragment_two = Strand()
        counter = 0
        for n in self._nucleotides:
            if counter < (len(self._nucleotides) / 2):
                fragment_one.add_nucleotide(n.copy() if copy else n)
            else:
                fragment_two.add_nucleotide(n.copy() if copy else n)
            counter += 1
        return fragment_one , fragment_two


def parse_visibility(path):
    try:
        inp = open (path, 'r')
    except:
        Logger.log ("Visibility file `" + path + "' not found. Assuming default visibility", Logger.WARNING)
        return []

    output = []
    for linea in inp.readlines():
        linea = linea.strip().lower()
        # remove everything that comes after '#'
        linea = linea.split('#')[0]
        if len(linea) > 0: output.append(linea)

    return output


class System(object):
    """
    Object representing an oxDNA system
    Contains strands

    Arguments:
    box -- the box size of the system
        Ex: box = [50, 50, 50]

    time --- Time of the system

    E_pot --- Potential energy

    E_kin --- Kinetic energy

    """

    def __init__(self, box, time=0, E_pot=0, E_kin=0):
        self._time = time
        self._ready = False
        self._box = np.array(box, np.float64)
        self._N = 0
        self._N_strands = 0
        self._strands = []
        self._nucleotide_to_strand = []
        self._N_cells = np.array(np.floor (self._box / 3.), np.int_)
        for kk in [0, 1, 2]:
            if self._N_cells[kk] > 100:
                self._N_cells[kk] = 100
        self._cellsides = box / self._N_cells
        self._head = [False, ] * int(self._N_cells[0] * self._N_cells[1] * self._N_cells[2])
        self.E_pot = E_pot
        self.E_kin = E_kin
        self.E_tot = E_pot + E_kin
        self.cells_done = False
        
        Nucleotide.index = 0
        Strand.index = 0

    def get_sequences (self):
        return [x._sequence for x in self._strands]

    _sequences = property (get_sequences)

    def get_N_Nucleotides(self):
        return self._N

    def get_N_strands(self):
        return self._N_strands

    def _prepare(self, visibility):
        sind = 0
        nind = 0
        for sind in range(self._N_strands):
            nind = self._strands[sind]._prepare(sind, nind)

        if visibility != None: self.set_visibility(visibility)

        for s in self._strands:
            s.bring_in_box_nucleotides(self._box)

    def copy (self):
        copy = System (self._box)
        for s in self._strands:
            copy.add_strand (s.copy (), check_overlap=False)
        return copy

    def get_reduced(self, according_to, bbox=None, check_overlap=False):
        visibility_list = self.get_visibility(according_to)

        if bbox == None or bbox == True:
            bbox = self._box
        elif isinstance(bbox, list) and not isinstance(bbox, np.array):
            bbox = np.array(bbox)
        else:
            Logger.die("Cannot reduce system, bbox not correct")

        copy = System(bbox)
        for i in range(self._N_strands):
            if visibility_list[i]:
                copy.add_strand(self._strands[i].copy(), check_overlap)

        return copy

    def join(self, other, box=None):
        if box is None:
            box = np.array([0., 0., 0.])
            for i in range(3):
                if other._box[i] > self._box[i]:
                    box[i] = other._box[i]
                else:
                    box[i] = self._box[i]

        ret = System(np.array(box, np.float64))
        for s in self._strands:
            if s.visible:
                ret.add_strand(s.copy(), check_overlap=False)
        for s in other._strands:
            if s.visible:
                ret.add_strand(s.copy(), check_overlap=False)

        return ret

    def get_visibility(self, arg=None):
        actions = {'vis': True, 'inv': False}
        visibility_list = [True, ] * self._N_strands

        if isinstance (arg, str):
            Logger.log ("Setting visibility with method 'file'", Logger.INFO)
            lines = parse_visibility(arg)

            for line in lines:
                # [uno, due, tre] = [p.strip() for p in line.partition("=")]
                [uno, due, tre] = [p.strip() for p in partition (line, "=")]
                if due != "=" or uno not in ["inv", "vis", "default"]:
                    Logger.log ("Lines in visibility must begin with one of inv=, vis= and default=. Skipping this line: --" + line + "--", Logger.WARNING)
                    continue

                if uno == 'default':
                    if tre not in ['inv', 'vis']:
                        Logger.log ("Wrong default in visibility file. Assuming visible as default", Logger.WARNING)
                        tre = 'vis'
                    if tre == 'inv':
                        visibility_list = [False, ] * self._N_strands
                else:
                    # filter removes all the empty strings
                    arr = [a.strip() for a in [_f for _f in tre.split(',') if _f]]
                    for a in arr:
                        try:
                            ind = int(a)
                        except:
                            Logger.log ("Could not cast '%s' to int. Assuming 0" % a, Logger.WARNING)
                            ind = 0
                        try:
                            visibility_list[ind] = actions[uno]
                        except:
                            Logger.log ("Strand %i does not exist in system, cannot assign visibility. Ignoring" % ind, Logger.WARNING)

        elif isinstance (arg, list):
            Logger.log("Setting visibility with method 'list'", Logger.INFO)
            # first len(arg) elements will be overwritten
            visibility_list[0:len(arg)] = arg
        else:
            if arg is not None:
                Logger.log("Argument of System.set_visibility can be a string or a list. Skipping visibility settings, assuming all strands are visible", Logger.WARNING)
            return visibility_list

        return visibility_list

    def set_visibility(self, arg=None):
        visibility_list = self.get_visibility(arg)

        for i in range(self._N_strands):
            self._strands[i].visible = visibility_list[i]

        return

    def do_cells (self):
        self._N_cells = np.array(np.floor (self._box / 3.), np.int_)
        for kk in [0, 1, 2]:
            if self._N_cells[kk] > 100:
                self._N_cells[kk] = 100
        self._cellsides = self._box / self._N_cells
        for n in self._nucleotides:
            n.next = -1
        self._head = [False, ] * int(self._N_cells[0] * self._N_cells[1] * self._N_cells[2])
        for n in self._nucleotides:
            cs = np.array((np.floor((n.cm_pos / self._box - np.rint(n.cm_pos / self._box) + 0.5) * (1. - FLT_EPSILON) * self._box / self._cellsides)), np.int_)
            cella = cs[0] + self._N_cells[0] * cs[1] + self._N_cells[0] * self._N_cells[1] * cs[2]
            n.next = self._head[cella]
            self._head[cella] = n
        self.cells_done = True
        return

    def add_strand(self, s, check_overlap=True):
        """
        Add a Strand to the System
        """
        '''
        # we now make cells off-line to save time when loading
        # configurations; interactions are computed with h_bonds.py
        # most of the time anyways
        for n in s._nucleotides:
            cs = np.array((np.floor((n.cm_pos/self._box - np.rint(n.cm_pos / self._box ) + 0.5) * (1. - FLT_EPSILON) * self._box / self._cellsides)), np.int_)
            cella = cs[0] + self._N_cells[0] * cs[1] + self._N_cells[0] * self._N_cells[1] * cs[2]
            n.next = self._head[cella]
            self._head[cella] = n
        '''
        self._strands.append(s)
        self._N += s.N
        self._N_strands += 1
        self.cells_done = False
        return True

    def add_strands(self, ss, check_overlap=True):
        if isinstance(ss, tuple) or isinstance(ss, list):
            added = []
            for s in ss:
                if self.add_strand(s, check_overlap):
                    added.append(s)
            if len(added) == len(ss):
                return True
            else:
                for s in added:
                    Nucleotide.index -= s.N
                    Strand.index -= 1
                    self._strands.pop()
                    self._N -= s.N
                    self._N_strands -= 1
                    self._sequences.pop()
                return False

        elif not self.add_strand(ss, check_overlap): 
            return False

        return True

    def get_unique_seq(self):
        # we need only the unique sequences of the system
        # see http://stackoverflow.com/questions/1143379/removing-duplicates-from-list-of-lists-in-python
        unique_seq = list(dict((str(x), x) for x in self._sequences).values())
        return unique_seq

    def rotate (self, amount, origin=None):
        for s in self._strands:
            s.rotate (amount, origin)

    def translate (self, amount):
        for s in self._strands:
            s.translate (amount)

    def print_ribbon_output(self, name, same_colors=False, visibility=None, constr_size=None):
        self._prepare(visibility)
        if constr_size != None:
            for i in range(self.N_strands / constr_size):
                cind = constr_size * i
                s1 = self._strands[cind]
                s2 = self._strands[cind + 1]
                s3 = self._strands[cind + 2]
                s4 = self._strands[cind + 3]

                diff1 = np.rint(s1.cm_pos / self._box) * self._box
                diff2 = np.rint(s2.cm_pos / self._box) * self._box
                diff3 = np.rint(s3.cm_pos / self._box) * self._box
                diff4 = np.rint(s4.cm_pos / self._box) * self._box

                s2.translate(np.rint((s1.cm_pos - s2.cm_pos - diff1 + diff2) / self._box) * self._box)
                s3.translate(np.rint((s1.cm_pos - s3.cm_pos - diff1 + diff3) / self._box) * self._box)
                s4.translate(np.rint((s1.cm_pos - s4.cm_pos - diff1 + diff4) / self._box) * self._box)

        unique_seq = list(self.get_unique_seq())

        if same_colors:
            n = len(unique_seq)
            colors = CREPY_COLOR_TABLE
            while len(colors) < n: colors *= 2

        f = open(name, "w")
        f.write(".Box:%lf,%lf,%lf\n" % tuple(self._box))
        for s in self._strands:
            out = s._get_ribbon_output() + "\n"
            if same_colors:
                color = colors[unique_seq.index(s.sequence)]
                out = out.replace("C[red]", "C[%s]" % color)
            f.write(out)
        f.close()

    def print_lorenzo_output(self, conf_name, top_name, visibility=None):
        self._prepare(visibility)
        conf = "t = %lu\nb = %f %f %f\nE = %lf %lf %lf\n" % (int(self._time), self._box[0], self._box[1], self._box[2], self.E_tot, self.E_pot, self.E_kin)

        visible_strands = 0
        visible_nucleotides = 0
        for s in self._strands:
            if s.visible:
                visible_strands += 1
                visible_nucleotides += s.N

        topology = "%d %d\n" % (visible_nucleotides, visible_strands)
        for s in self._strands:
            sc, st = s._get_lorenzo_output()
            topology += st
            conf += sc
            
        with open(conf_name, "w") as f:
            f.write(conf)
            f.close()
        with open(top_name, "w") as f:
            f.write(topology)
            f.close()

    def print_oxview_output(self, name):
        import json

        out = {
            'box': self._box.tolist(),
            'systems': [{'id':0, 'strands': []}]
        }

        for s in self._strands:
            strand = {
                'id': s.index, 'end3': s._nucleotides[-1].index, 'end5': s._nucleotides[0].index,
                'class': 'NucleicAcidStrand', 'monomers': []
            }
            for i, n in enumerate(s._nucleotides):
                if s._circular:
                    if i == 0:
                        n5 = s._nucleotides[-1].index
                    else:
                        n5 = s._nucleotides[i - 1].index
                    if i == len(s._nucleotides) - 1:
                        n3 = s._nucleotides[0].index
                    else:
                        n3 = s._nucleotides[i + 1].index
                else:
                    if i == 0:
                        n5 = -1
                    else:
                        n5 = s._nucleotides[i - 1].index
                    if i == len(s._nucleotides) - 1:
                        n3 = -1
                    else:
                        n3 = s._nucleotides[i + 1].index
                nucleotide = {
                    'id': n.index,
                    'type': n.get_base(),
                    'class': 'DNA',
                    'p': n.cm_pos.tolist(),
                    'a1': n._a1.tolist(),
                    'a3': n._a3.tolist()
                }
                if n3 >= 0:
                    nucleotide['n3'] = n3
                if n5 >= 0:
                    nucleotide['n5'] = n5
                if n.pair is not None:
                    nucleotide['bp'] = n.pair.index
                if n.cluster is not None:
                    nucleotide['cluster'] = n.cluster
                if n.color is not None:
                    nucleotide['color'] = n.color
                strand['monomers'].append(nucleotide)
            out['systems'][0]['strands'].append(strand)

        with open(name, "w") as f:
            f.write(json.dumps(out))

    N = property(get_N_Nucleotides)
    N_strands = property (get_N_strands)

    def get_nucleotide_list (self):
        ret = []
        for s in self._strands:
            ret += s._nucleotides
        return ret

    _nucleotides = property (get_nucleotide_list)

    def map_nucleotides_to_strands(self):
        # this function creates nucl_id -> strand_id array
        for i in range(len(self._strands)):
            for _ in range(self._strands[i].get_length()):
                self._nucleotide_to_strand.append(i)

    def print_dot_bracket_output(self, filename):
        # assumes each nucleotide has at most 1 hydrogen bond, requires interactions already to be filled for nucleotide objects
        nupack_string = ""
        for n1 in range(self.get_N_Nucleotides()):
            interactions = self._nucleotides[n1].interactions
            if len(interactions) > 1:
                Logger.log ("more than 1 HB for a nucleotide", Logger.WARNING)
            if len(interactions) == 0:
                nupack_string += "."
            elif interactions[0] > n1:
                nupack_string += "("
            elif interactions[0] < n1:
                nupack_string += ")"
            else:
                Logger.log("unexpected interaction detected while building nupack string", Logger.CRITICAL)

        f = open(filename, "w")
        f.write(nupack_string)
        f.close()
        
