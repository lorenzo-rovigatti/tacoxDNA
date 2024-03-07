#!/usr/bin/env python3

# This script converts cando format to oxDNA2 files
import sys
import numpy as np

from .libs.readers import LorenzoReader

BASE_SHIFT = 1.13
COM_SHIFT = 0.54
# COM_SHIFT = 0.545

SCALE = (1 / 0.85) * (1 / 10.)


class Base:

    def __init__(self, base_id, across, down, up, seq):
        self.id = base_id
        self.oxid = -1
        self.across = across
        self.up = up
        self.down = down
        self.seq = seq

        self.cm_pos = None
        self.a1 = np.array([1., 0, 0])
        self.a3 = np.array([0, 1., 0])

    def set_position(self, cm_bp, e1_bp, e2_bp, e3_bp, preferred=True):
        if  not preferred:
            self.cm_pos = SCALE * cm_bp - COM_SHIFT * e2_bp
            self.a1 = e2_bp
            self.a3 = e3_bp
        else:
            self.cm_pos = SCALE * cm_bp + COM_SHIFT * e2_bp
            self.a1 = -e2_bp
            self.a3 = -e3_bp

    def get_5prime(self):
        return self.up

    def pos_str(self):
        if self.cm_pos is None:
            print("Warning, base ", self.id, self.oxid, 'was not assigned cm', file=sys.stderr)
        if self.cm_pos[0] == 0 and self.cm_pos[1] == 0 and self.cm_pos[2] == 0:
            print("Warning, cm of base ", self.id, self.oxid, 'was set to 0 0 0', file=sys.stderr)

        cm = self.cm_pos
        a1 = self.a1
        a3 = self.a3
        s = '%f %f %f %f %f %f %f %f %f 0.0 0.0 0.0 0.0 0.0 0.0' % (
            cm[0], cm[1], cm[2], a1[0], a1[1], a1[2], a3[0], a3[1], a3[2])

        return s


def normalize(v):
    return v / np.sqrt(np.dot(v, v))


class Strand:
    base_counter = 0

    def __init__(self, s_id):
        self.id = s_id
        self.bases = []

    def add_base_at_5_prime(self, b):
        b.oxid = Strand.base_counter
        Strand.base_counter += 1
        self.bases.append(b)

    def __len__(self):
        return len(self.bases)

    def fix_nans(self):
        # some bases are not assigned to any bp, in which case their center of mass needs to be interpolated
        last_assigned = None
        unassigned = []
        for b in self.bases:
            if not b.cm_pos is None:
                if len(unassigned) == 0:
                    last_assigned = b.cm_pos
                else:
                    new_assigned = b.cm_pos
                    N = len(unassigned)
                    v = (np.array(new_assigned) - 
                         np.array(last_assigned)) / float(N)
                    spos = np.array(last_assigned)
                    for ub in unassigned:
                        # print 'Assigning to ',ub.oxid, last_assigned + v
                        ub.cm_pos = spos
                        spos = spos + v
                    unassigned = []
            else:
                unassigned.append(b)

    def better_fix_nans(self):
        # some bases are not assigned to any bp, in which case their center of mass needs to be interpolated
        # this version deals with overlaps by pointing the unpaired sections away from the c.of.m.
        last_assigned = None
        unassigned = []
        cm = np.zeros(3)
        counter = 0.
        for b in self.bases:
            if not b.cm_pos is None:
                cm += np.array(b.cm_pos)
                counter += 1
        cm /= counter

        for b in self.bases:
            if not b.cm_pos is None:
                if len(unassigned) == 0:
                    last_assigned = b.cm_pos
                else:
                    new_assigned = b.cm_pos
                    v = np.array(new_assigned) - cm
                    v = normalize(v)

                    # v = (np.array(new_assigned) -
                    #     np.array(last_assigned))/float(N)
                    spos = np.array(last_assigned) 
                    for ub in unassigned:
                        # print 'Assigning to ',ub.oxid, last_assigned + v
                        ub.cm_pos = spos
                        spos = spos + 0.8 * v 
                    unassigned = []
            else:
                unassigned.append(b)


def reconstruct_strands(three_ends, bases):
    strands = []
    strand_id = 0
    added = []
    for s_id in three_ends:
        s = Strand(strand_id)
        strand_id += 1
        b = bases[s_id]
        added.append(s_id)
        s.add_base_at_5_prime(b)
        s_next = b.get_5prime()
        while s_next != -1:
            b = bases[s_next]
            s.add_base_at_5_prime(b)
            added.append(s_next)
            s_next = b.get_5prime()
        strands.append(s)
    
    if(len(bases) != len(added)):  # maybe circular strands?
        s = None
        for b_id, b in list(bases.items()):
            if b_id not in added:  # maybe new strand:
                s = Strand(strand_id)
                start = b_id
                strand_id += 1
                s.add_base_at_5_prime(b)
                added.append(b_id)
                b_next = b.get_5prime()
                while b_next != -1 and b_next != start and b_next not in added:
                    nb = bases[b_next]
                    s.add_base_at_5_prime(nb)
                    added.append(b_next)
                    b_next = nb.get_5prime()
                    
                    if len(added) > len(list(bases.keys())):
                        raise Exception("error while processing circular strand")
                strands.append(s)

    return strands


def write_topology(strands, outfile):
    handle = open(outfile, 'w')
    total_p = int(np.sum([len(s) for s in strands]))
    total_s = len(strands)
    print(total_p, total_s, file=handle)
    for s in strands:
        for i, b in enumerate(s.bases):
            if i == 0:
                three = -1
            else:
                three = s.bases[i - 1].oxid
            if i == len(s.bases) - 1:
                five = -1
            else:
                five = s.bases[i + 1].oxid
            print(s.id + 1, b.seq, three, five, file=handle)


def write_conf(bases, fname, box):
    handle = open(fname, 'w')
    print("t = 0", file=handle)
    print("b = %f %f %f" % (box, box, box), file=handle)
    print("E = 0 0 0", file=handle)
    ox_bases = {}
    for i, b in list(bases.items()):
        ox_bases[b.oxid] = b

    for i in range(max(ox_bases.keys()) + 1):
        handle.write(ox_bases[i].pos_str() + '\n')

    handle.close()

def write_force(bases, pairs, fname):
    handle = open(fname, 'w')
    forces = ""
    for p in list(pairs.values()):
        a = bases[p[0]].oxid
        b = bases[p[1]].oxid
        forces = forces + get_mutual_force(a, b)   
    handle.write(forces)
    handle.close()


def assign_base_positions(bases, cm_pos, triads, pair_ids):
    for pid, pair in list(pair_ids.items()):
        aid = pair[0]
        bid = pair[1]
        cm = cm_pos[pid]
        e1 = triads[pid][0]
        e2 = triads[pid][1]
        e3 = triads[pid][2]
        bases[aid].set_position(cm, e1, e2, e3, True)
        bases[bid].set_position(cm, e1, e2, e3, False)


def get_mutual_force(id1, id2, stiffness=0.1):
    s = "{\n type = mutual_trap\n particle = %d\n ref_particle = %d\n stiff = %f \n r0 = 1.2 \n}\n" % (id1, id2, stiffness)
    s = s + "{\n type = mutual_trap\n particle = %d\n ref_particle = %d\n stiff = %f \n r0 = 1.2 \n}\n" % (id2, id1, stiffness)
    return s


def load_and_convert(opts, invert_preference=False):
    cando_file = opts['cando_file']
    box = opts['box']
    # loads cando file into oxDNA system
    infile = open(cando_file, 'r')
    lines = infile.readlines()
    top_start = -1
    top_end = -1
    for i, line in enumerate(lines):
        if 'dnaTop' in line:
            top_start = i + 1

        if 'dNode' in line:
            top_end = i - 1
            break

    node_cm_pos = {}
    node_start = -1
    node_end = -1

    # do nodes
    for i, line in enumerate(lines):
        if 'dNode' in line:
            node_start = i + 1

        if 'triad' in line:
            node_end = i - 1
            break
    for line in lines[node_start:node_end]:
        nodes_info = line.strip().split(',')
        # print nodes_info
        nid = int(nodes_info[0])
        node_pos = [float(v) for v in nodes_info[1:]]
        node_pos = np.array(node_pos)
        node_cm_pos[nid] = node_pos

    # do triads
    triads = {}
    triad_start = -1
    triad_end = -1
    for i, line in enumerate(lines):
        if "triad," in line:
            triad_start = i + 1
        elif "id_nt," in line:  # or ',' not in line:
            triad_end = i - 1
            break

    # print triad_start,triad_end
    for line in lines[triad_start:(triad_end)]:
        vals = line.strip().split(',')

        # print vals
        t_id = int(vals[0])
        vectors = [float(v) for v in vals[1:]]
        e1 = np.array(vectors[0:3])
        e2 = np.array(vectors[3:6])
        e3 = np.array(vectors[6:9])
        triads[t_id] = [e1, e2, e3]

    # find paired nucletides
    pair_ids = {}
    pairs_start = -1
    pairs_end = -1
    for i, line in enumerate(lines):
        if "id_nt," in line:
            pairs_start = i + 1
        elif pairs_start > -1 and len(line.strip().split(',')) < 3:
            pairs_end = i
            break

    if pairs_end == -1:
        pairs_end = len(lines)
    # print pairs_start, pairs_end

    for line in lines[pairs_start:pairs_end]:
        vals = line.strip().split(',')
        # print vals
        pid = int(vals[0])
        nuca = int(vals[1])
        nucb = int(vals[2])
        if invert_preference:
            pair_ids[pid] = [nucb, nuca]
        else:
            pair_ids[pid] = [nuca, nucb]

    bases = {}
    three_primes = []
    for topline in lines[top_start:top_end]:
        vals = topline.strip().split(',')
        base_id = int(vals[1])
        up = int(vals[2])  # up is the 5'
        down = int(vals[3])  # down is the 3'
        across = int(vals[4])
        seq = vals[5].strip()
        b = Base(base_id, across, down, up, seq)
        if down == -1:
            three_primes.append(base_id)
        bases[base_id] = b

    # print three_primes
    # print bases.keys()
    strands = reconstruct_strands(three_primes, bases)
    write_topology(strands, cando_file + '.top')

    assign_base_positions(bases, node_cm_pos, triads, pair_ids)
    for s in strands:
        s.better_fix_nans()

    write_conf(bases, cando_file + '.oxdna', box)

    if opts['print_force_file']:
        force_file = cando_file + '.forces.txt'
        print("## Printing forces to the '%s' file" % force_file, file=sys.stderr)
        write_force(bases, pair_ids, force_file)
        
    return strands


def print_usage():
        print("USAGE:", file=sys.stderr)
        print("\t%s CanDo_file" % sys.argv[0], file=sys.stderr)
        print("\t[-b\--box=100]", file=sys.stderr)
        print("\t[-f\--print-force-file]\n", file=sys.stderr)
        exit(1)


def parse_options():
    shortArgs = 'b:f'
    longArgs = ['box=', 'print-force-file']
    
    # for some reason, files originally made in T1 have a different .json form than T2
    # it would be possible to rewrite all the parameters to fix it, but tossing a factor
    # of 1.2 for DNA and 1.6 for RNA on base_vector seems to be good enough
    opts = {
        "box" : 100.,
        "print_force_file" : False
    }
    
    try:
        import getopt
        args, positional_args = getopt.gnu_getopt(sys.argv[1:], shortArgs, longArgs)
        for k in args:
            if k[0] == '-b' or k[0] == '--box':
                try:
                    opts['box'] = float(k[1])
                    print("## Setting the box size to %f" % opts['box'], file=sys.stderr)
                except ValueError:
                    print("The argument of '%s' should be a number (got '%s' instead)" % (k[0], k[1]), file=sys.stderr)
                    exit(1)
            elif k[0] == '-f' or k[0] == '--print-force-file':
                opts['print_force_file'] = True
            
        opts['cando_file'] = positional_args[0]
        
    except Exception:
        print_usage()
        
    return opts

def main():
    opts = parse_options()
    
    s = load_and_convert(opts)
    r = LorenzoReader(opts['cando_file'] + '.top', opts['cando_file'] + '.oxdna')
    mys = r.get_system()
    refdir = mys._strands[0]._nucleotides[1].get_pos_base() - mys._strands[0]._nucleotides[0].get_pos_base()
    refdir /= np.sqrt(np.dot(refdir, refdir))
    ref_angle = np.dot(mys._strands[0]._nucleotides[0]._a3, refdir) 
    if ref_angle < 0 and abs(ref_angle) > 0.8:  # likely selected the wrong orientation for preferred nucleotide!
        Strand.base_counter = 0
        s = load_and_convert(opts, invert_preference=True)

    print("## Wrote data to '%s' / '%s'" % (opts['cando_file'] + '.oxdna', opts['cando_file'] + '.top'), file=sys.stderr)
    print("## DONE", file=sys.stderr)

if __name__ == '__main__':
    main()

