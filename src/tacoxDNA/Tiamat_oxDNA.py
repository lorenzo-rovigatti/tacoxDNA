#!/usr/bin/env python3

import sys
from json import load
import numpy as np
import random
import json

DNA_BASES = ['A', 'C', 'G', 'T']

# TODO: use the stuff in base.py
class Base(object):
    scale = 1 / 0.85
    default_val = 'R'
    
    def __init__(self, local_id, base_info):
        self.base_info = base_info
        self.global_id = base_info['id']
        self.local_id = local_id
        self.up = None
        self.down = None
        self.across = None
        try:
            self.val = base_info['type'][0]
        except:
            if Base.default_val == 'R':
                self.val = random.choice(DNA_BASES)
            else:
                self.val = Base.default_val
            print("WARNING: base n.%s has no associated type, setting it to '%s'" % (self.global_id, self.val), file=sys.stderr)
        self.pos = np.array(base_info['position'])

    def __get_connected_id(self, neighbor):
        if  neighbor == None:
            return -1
        return neighbor.local_id

    def get_up_id(self):
        return self.__get_connected_id(self.up)

    def get_down_id(self):
        return self.__get_connected_id(self.down)

    def get_across_id(self):
        return self.__get_connected_id(self.across)

    def __str__(self):
        return "%d %s " % (self.local_id, self.val)

    def get_pos(self):
        return Base.scale * self.pos

    def get_up(self):
        return self.up

    def get_down(self):
        return self.down
    
    def get_across(self):
        return self.across
    
    
class Strand(object):

    def __init__(self, strand_id, base_info_lst, local_id_strart=0):
        self.strand_id = strand_id
        self.bases = [Base(local_id, base_info) 
                        for local_id, base_info in
                         enumerate(base_info_lst, local_id_strart)]

    def __len__(self):
        return len (self.bases)

    def __str__(self):
        return "Strand %d: %s" % (self.strand_id,
                                    [str(b) for b in  self.bases])
        

def normalize(v):
    return v / np.sqrt(np.dot(v, v))


def get_5primes(bases):
    """  5' ends don't have bases downstream """
    return [b for b in bases if b['down'] is None]


def get_circular_strands(bases, already_included, start_id, loc):
    cs = [b for b in bases if not b['id'] in already_included]
    bases_by_id = {base['id'] : base for base in bases}
    idx = start_id
    strands = []
    local_id = loc

    while len(cs) > 0:
        base_lst = [cs[0]]
        base_lst_ids = [cs[0]['id']]
        next_base = bases_by_id[cs[0]['up']]
        while next_base['id'] not in base_lst_ids:
            base_lst.append(next_base)
            base_lst_ids.append(next_base['id'])
            next_base = bases_by_id[next_base['up']]
        strand = Strand(idx, base_lst, local_id)
        strands.append(strand)
        local_id = strand.bases[-1].local_id + 1
        idx += 1
        cs = [b for b in cs if b['id'] not in base_lst_ids]
    
    return strands


def locate_base_on_strands(strands, global_id):
    for strand in strands:
        for base in strand.bases:
            if base.global_id == global_id:
                return base  
    return NoBase()


# creates a list of all strands in the respective tiamat file 
def split_bases_to_strands(bases):
    # lookup bases by id 
    bases_by_id = {base['id'] : base for base in bases}
    strands = []
    local_id = 0 
    # for #5'
    added = 0
    included = []
    for idx, start_base in enumerate(get_5primes(bases), 1):
        base_lst = [start_base]
        base = start_base   
        # follows every strand up the backbone (5'->3')
        while not base['up'] is None:
            base = bases_by_id[base['up']]
            base_lst.append(base)
            added += 1
        strand = Strand(idx, base_lst, local_id)
        strands.append(strand)
        included = included + [b['id'] for b in base_lst]
        local_id = strand.bases[-1].local_id + 1
    
    extra_strands = get_circular_strands(bases, included, len(strands) + 1, local_id)
    return strands + extra_strands


# define the base connections across strands 
def define_connections(strands):
    # lookup bases by global id 
    for strand in strands:
        for base in strand.bases:
            base.across = locate_base_on_strands(strands, base.base_info['across'])
            base.up = locate_base_on_strands(strands, base.base_info['up'])
            base.down = locate_base_on_strands(strands, base.base_info['down'])


# used to calculate the center of mass 
def neighbor3_cal_vector(AB, AA3, AB5):
    centrmas = normalize(-0.13079674 * AB - 0.22543211 * AA3 + 0.62949112 * AB5)
    a3 = normalize(2.69498211 * AB - 1.04531113 * AA3 - 2.30531223 * AB5)
    a1 = AB
    return [a1, a3, centrmas]


# used to calculate the center of mass 
def neighbor5_cal_vector(AB, AA5, AB3):
    centrmas = normalize(0.81079674 * AB + 0.22543211 * AA5 - 0.50262804 * AB3)
    a3 = normalize(-2.12846367 * AB + 0.82557385 * AA5 + 2.33064701 * AB3)
    a1 = AB
    return [a1, a3, centrmas]


def neighbor3_cal_vector_RNA(AB, AA3, AB5):
    centrmas = normalize(-0.28102082 * AB - 0.25891019 * AA3 + 0.84990909 * AB5)
    a3 = normalize(2.34763359 * AB - 1.1627428 * AA3 - 1.63537381 * AB5)
    a1 = AB
    return [a1, a3, centrmas]


def neighbor5_cal_vector_RNA(AB, AA5, AB3):
    centrmas = normalize(0.85540635 * AB + 0.30569283 * AA5 - 0.44567833 * AB3)
    a3 = normalize(-1.60523423 * AB + 0.58820649 * AA5 + 2.00150202 * AB3)
    a1 = AB
    return [a1, a3, centrmas]


class NoBase(Base):

    def __init__(self):
        self.local_id = -1
        self.val = None 

    def __str__(self):
        return "%d %s " % (self.local_id, self.val)

    def get_pos(self):
        return Base.scale * np.array([0, 0, 1])

    def get_across(self):
        return self  # Returns itself to reduce memory consumtion 

def get_base_config(base):
    base_vector = base.get_pos()

    pairing_base = base.get_across()
    paring_base_vector = pairing_base.get_pos()

    up_base = base.get_up()
    up_base_vector = up_base.get_pos()

    down_base = base.get_down()
    down_base_vector = down_base.get_pos()

    paring_base3 = up_base.get_across()
    paring_base3_vector = paring_base3.get_pos()

    paring_base5 = down_base.get_across()
    paring_base5_vector = paring_base5.get_pos()

    # three backbone vectors, A-base_vector, B-paring_base_vector,
    # A5up_base_vector, A3down_base_vector
    backbone_A_to_backbone_B = normalize(-base_vector + paring_base_vector)
    backbone_A_to_posA5_neibor = normalize(-base_vector + up_base_vector)
    backbone_A_to_posA3_neibor = normalize(-base_vector + down_base_vector)
    backboneA_to_backbone_B3 = normalize(-base_vector + paring_base3_vector)
    backboneA_to_backbone_B5 = normalize(-base_vector + paring_base5_vector)

    # do we have a duplex across the 3' end
    if not(type(paring_base5) is NoBase):
        if opts['isDNA']:
            a1_vector, a3_vector, cm_pos = neighbor3_cal_vector(backbone_A_to_backbone_B, backbone_A_to_posA3_neibor, backboneA_to_backbone_B5)
        else:
            a1_vector, a3_vector, cm_pos = neighbor3_cal_vector_RNA(backbone_A_to_backbone_B, backbone_A_to_posA3_neibor, backboneA_to_backbone_B5)
        cm_pos = cm_pos + (base_vector * opts['tiamat_version_fudge'])
    # do we have a duplex across the 5' end
    elif not(type(paring_base3) is NoBase):
        if opts['isDNA']:
            a1_vector, a3_vector, cm_pos = neighbor5_cal_vector(backbone_A_to_backbone_B, backbone_A_to_posA5_neibor, backboneA_to_backbone_B3)
        else:
            a1_vector, a3_vector, cm_pos = neighbor5_cal_vector_RNA(backbone_A_to_backbone_B, backbone_A_to_posA5_neibor, backboneA_to_backbone_B3)
        cm_pos = cm_pos + (base_vector * opts['tiamat_version_fudge'])
    else:
        # single stranded case, not really treated (could randomize orientations?)
        cm_pos = base_vector
        a3_vector = np.array([0., 0., 1.])
        a1_vector = np.array([0., 1., 0.])

    # assemble a new line in the configuration file
    return cm_pos, a1_vector, a3_vector

def write_topology_file(strands, topology_file):
    # setup the topology header 
    top_lines = [
        '%d %d' % (len(bases), len(strands))  # number of bases, number of strands
    ]

    # go through all the strands to generate topology 
    for strand in strands:
        for base in strand.bases:
            top_lines.append(
                '%d %s %d %d' % (strand.strand_id, base.val, base.get_down_id(), base.get_up_id())
            )

    # spit out topology file
    with open(topology_file, "w") as f_out:
        f_out.write('\n'.join(top_lines))


def write_force_file(strands, force_file_name):
    mutual_trap_template = '{ \ntype = mutual_trap\nparticle = %d\nstiff = 0.9\nr0 = 1.2\nref_particle = %d\nPBC=1\n}\n'
    lines = []
    for strand in strands:
        for base in strand.bases:
            from_particle_id = base.local_id
            to_particle_id = base.get_across().local_id
            if to_particle_id >= 0:
                lines.append(mutual_trap_template % (from_particle_id, to_particle_id))
    
    with open(force_file_name, "w") as f_out:
        f_out.write("\n".join(lines))

def write_oxview_file(strands,file_name):
    box_size = 1
    out = {
        'systems': [{'id':0, 'strands': []}]
    }

    for s in strands:
        strand = {
            'id': s.strand_id, 'end3': s.bases[0].local_id, 'end5': s.bases[-1].local_id,
            'class': 'NucleicAcidStrand',
            'monomers': []
        }
        for n in s.bases:
            cm_pos, a1_vector, a3_vector = get_base_config(n)

            # Make sure box fits base
            for v in cm_pos.tolist():
                if abs(v) > box_size/2:
                    box_size = abs(v)*2

            n3 = n.get_down_id()
            n5 = n.get_up_id()
            pair = n.get_across_id()

            nucleotide = {
                'id': n.local_id,
                'type': n.val,
                'class': 'DNA' if opts['isDNA'] else 'RNA',
                'p': cm_pos.tolist(),
                'a1': a1_vector.tolist(),
                'a3': a3_vector.tolist()
            }
            if n3 >= 0:
                nucleotide['n3'] = n3
            if n5 >= 0:
                nucleotide['n5'] = n5
            if pair >= 0:
                nucleotide['bp'] = pair
            strand['monomers'].append(nucleotide)
        out['systems'][0]['strands'].append(strand)

    box_size *= 1.5 # Add some margin
    box_size = round(box_size)
    out['box'] = [box_size, box_size, box_size]

    with open(file_name, "w") as f:
        f.write(json.dumps(out))


def write_configuration_file(strands, configuration_file, opts):
    # Decide on a box size based on strand length and number of strands - there's not a good way to do this blind
    # Need to know strand lengths, organization, shape
    # strand_lengths = map(len, strands)
    
    # TODO: choose the size of the box according to the configuration we are converting
    box_size = 250  
    # box_size = max(strand_lengths) * 2
    
    # setup the configuration file header
    configuration_lines = [
        't = %d' % 0,  # set time 
        'b = %d %d %d' % (box_size, box_size, box_size),  # set box size
        'E = 0 0 0'  # set energy
    ]
    
    for strand in strands:
        for base in strand.bases:
            cm_pos, a1_vector, a3_vector = get_base_config(base)
            configuration_lines.append(
                "".join([
                    "%f %f %f " % (cm_pos[0], cm_pos[1], cm_pos[2]),
                    "%f %f %f " % (a1_vector[0], a1_vector[1], a1_vector[2]),
                    "%f %f %f " % (a3_vector[0], a3_vector[1], a3_vector[2]),
                    "0 0 0 ",
                    "0 0 0" 
                ])
            )
    # spit out the configuration file 
    with open(configuration_file, 'w') as f_out: 
        f_out.write('\n'.join(configuration_lines) + "\n")
    

def print_usage():
        print("USAGE:", file=sys.stderr)
        print("\t%s Tiamat_json_file" % sys.argv[0], file=sys.stderr)
        print("\t[-m\--molecule=DNA|RNA]", file=sys.stderr)
        print("\t[-t\--tiamat-version=1|2]", file=sys.stderr)
        print("\t[-t\--default-base=A|C|G|T|R|i (R = random, i = any integer)]", file=sys.stderr)
        print("\t[-f\--print-force-file]", file=sys.stderr)
        print("\t[-o\--print-oxview]\n", file=sys.stderr)
        print("\tThe defaults options are --molecule=DNA, --tiamat-version=2, --default-base=R", file=sys.stderr)
        exit(1)


def parse_options():
    shortArgs = 'd:m:t:f:o'
    longArgs = ['default-base=', 'molecule=', 'tiamat-version=', 'print-force-file', 'print-oxview']
    
    # for some reason, files originally made in T1 have a different .json form than T2
    # it would be possible to rewrite all the parameters to fix it, but tossing a factor
    # of 1.2 for DNA and 1.6 for RNA on base_vector seems to be good enough
    opts = {
        "isDNA" : True,
        "tiamat_version_fudge" : 1,
        "print_force_file" : False,
        "print_oxview" : False
    }
    
    tiamat_version = 2
    try:
        import getopt
        args, positional_args = getopt.gnu_getopt(sys.argv[1:], shortArgs, longArgs)
        for k in args:
            if k[0] == '-m' or k[0] == '--molecule':
                k_arg = k[1].upper()
                if k_arg == "DNA":
                    print("## Assuming DNA", file=sys.stderr)
                    opts["isDNA"] = True
                elif k_arg == "RNA":
                    print("## Assuming RNA", file=sys.stderr)
                    opts["isDNA"] = False
                else:
                    print("The argument of '%s' should be either 'DNA' or 'RNA' (got '%s' instead)" % (k[0], k[1]), file=sys.stderr)
                    exit(1)
            elif k[0] == '-t' or k[0] == '--tiamat-version':
                tiamat_version = int(k[1])
            elif k[0] == '-d' or k[0] == '--default-base':
                db = k[1].upper()
                if db != 'R' and db not in DNA_BASES:
                    try:
                        int(k[1])
                    except Exception:
                        print("Invalid default base value '%s'. The only accepted values are %s, R or an integer)" % (k[1], ", ".join(DNA_BASES)), file=sys.stderr)
                        exit(1)
                        
                Base.default_val = db
            elif k[0] == '-f' or k[0] == '--print-force-file':
                opts['print_force_file'] = True
            elif k[0] == '-o' or k[0] == '--print-oxview':
                opts['print_oxview'] = True
            
        if tiamat_version == 1:
            if opts['isDNA']:
                opts['tiamat_version_fudge'] = 1.2
            else: 
                opts['tiamat_version_fudge'] = 1.6
        elif tiamat_version == 2:
            opts['tiamat_version_fudge'] = 1
        else:
            print("The argument of '%s' should be either '1' or '2' (got '%s' instead)" % (k[0], k[1]), file=sys.stderr)
            exit(1)

        opts['tiamat_file'] = positional_args[0]
            
        print("## Assuming Tiamat version %d" % tiamat_version, file=sys.stderr)
        
    except Exception:
        print_usage()
    
    return opts


def main():
    opts = parse_options()

    topology_file = opts['tiamat_file'] + ".top"
    configuration_file = opts['tiamat_file'] + ".oxdna"
    force_file_name = opts['tiamat_file'] + ".forces.txt"
    oxview_file_name = opts['tiamat_file'] + ".oxview"
    
    # read and parse json
    with open(opts['tiamat_file'], 'r') as f_out:
        data = load(f_out)
    # our main interest are the bases
    bases = data['bases']
    
    # get a list of the strands in the system
    strands = split_bases_to_strands(bases)
    # setup connections with new local indices within the strands  
    define_connections(strands)
    
    # # build topology file first 
    write_topology_file(strands, topology_file)

    # write forces file
    if opts['print_force_file']:
        write_force_file(strands, force_file_name)
    
    # work on the configuration file 
    write_configuration_file(strands, configuration_file, opts)

    # write oxview file
    if opts['print_oxview']:
        write_oxview_file(strands, oxview_file_name)
    
    print("## Wrote data to '%s' / '%s'" % (configuration_file, topology_file), file=sys.stderr)
    print("## DONE", file=sys.stderr)



if __name__ == '__main__':
    main()