import numpy as np
import sys, os, subprocess, math, re
from collections import OrderedDict


def base_identify(X):
    base_id = ["A", "T", "C", "G"]
    complementary_base_id = ["T", "A", "G", "C"]
    return base_id[X], complementary_base_id[X]


def GetMatrix(gamma, beta, alpha):
    gamma = math.radians(gamma)
    beta = math.radians(beta)
    alpha = -math.radians(alpha)
    s1 = math.sin(alpha)
    s2 = math.sin(beta)
    s3 = math.sin(gamma)
    c1 = math.cos(alpha)
    c2 = math.cos(beta)
    c3 = math.cos(gamma)
    a = (c1 * c2)
    b = (c1 * s3 * s2) + (s1 * c3)
    c = (c1 * s2 * c3) - (s1 * s3)
    d = 0 - (s1 * c2)
    e = 0 - (s1 * s3 * s2) + (c1 * c3)
    f = 0 - (s1 * s2 * c3) - (c1 * s3)
    g = 0 - s2
    h = (s3 * c2)
    i = (c2 * c3)
    R = np.array([[a, b, c], [d, e, f], [g, h, i]])
    return R


def export_oxDNA(ma_file,box_size):
    base_parents = OrderedDict()
    base_names = []
    base_types = []
    vHelix_names = []
    first_base_parents = OrderedDict()
    dups = []
    helix = 0
    helix_number = 0
    helix_t = np.array([0, 0, 0])
    helix_r = np.array([0, 0, 0])
    helix_trans = {}
    helix_rot = {}
    base = 0
    first_base_types = {}
    first_base_trans = {}
    base_types = {}
    base_trans = {}
    base_translation = np.array([1, 1, 1])
    base_type = 5
    base_number = 0


    first_last_base = np.zeros((5,), dtype=np.int)  # ERIK

    vHelixFile = open(ma_file, "r").readlines()
    for line in vHelixFile:
        if "createNode" in line:
            if "vHelix" in line:
                helix_number += 1
                vHelix_names.append(line.split(" ")[3].strip('"').rstrip('";\n'))
                base = 0  # ERIK seems to be used to keep track if the current item is a helix or a base
                helix = 1  # ERIK seems to be used to keep track if the current item is a helix or a base
                if helix_number > 1:
                    helix_trans["%s" % helix_name] = helix_t
                    helix_rot["%s" % helix_name] = helix_r
                    helix_t = np.array([0, 0, 0])
                    helix_r = np.array([0, 0, 0])
                helix_name = line.split(" ")[3].strip('"').rstrip('";\n')
            if "HelixBase" in line:
                base_number += 1
                helix = 0  # ERIK seems to be used to keep track if the current item is a helix or a base
                base = 1  # ERIK seems to be used to keep track if the current item is a helix or a base
                # if '"forw_2" -p "helix_' in line:
                #     temp_helix_nr = str(re.findall('helix_\d+',line)[0])
                #     print 'found forward '+ str(base_number)+' '+temp_helix_nr
                #
                # if '"backw_2" -p "helix_' in line:
                #     temp_helix_nr = str(re.findall('helix_\d+',line)[0])
                #     print 'found backward ' + str(base_number)+' '+temp_helix_nr

                if base_number > 1:
                    if float(base_translation[0]) == 0. and float(base_translation[1]) == 0. and float(base_translation[2]) == 0.:
                        print("defaulting")
                    first_base_trans["|%s|%s" % (base_parent, base_name)] = base_translation
                    first_base_types["|%s|%s" % (base_parent, base_name)] = base_type
                    # base_translation=np.array([0,0,0])
                    base_type = 5
                try:
                    base_parent = line.split(" ")[5].strip('"').rstrip('";\n')
                    base_name = line.split(" ")[3].strip('"').rstrip('"')
                    first_base_parents["|%s|%s" % (base_parent, base_name)] = base_parent
                except:
                    pass

        if '".t"' in line and helix == 1:
            try:
                data = line.split(" ")[4:7]
                helix_t = np.array([float(data[0]), float(data[1]), float(data[2])])
            except:
                print(line)
                sys.exit()
        elif '".t"' in line and base == 1:
            # print "using base"
            data = line.split(" ")[4:7]
            base_translation = np.array([float(data[0]), float(data[1]), float(data[2])])
        if '".r"' in line and helix == 1:
            try:
                data = line.split(" ")[4:7]
                helix_r = np.array([float(data[0]), float(data[1]), float(data[2])])
            except:
                print(line)
                sys.exit()
        if '".lb"' in line and base == 1:
            try:
                base_type = int(line.split(" ")[2].rstrip(";\n"))
            except:
                print(line.split(" ")[2].rstrip(";\n"))
                sys.exit()

        if "createNode" in line and "lightLinker" in line:
            # print "resetting values"
            helix_trans["%s" % helix_name] = helix_t
            helix_rot["%s" % helix_name] = helix_r
            helix_t = np.array([0, 0, 0])
            helix_r = np.array([0, 0, 0])
            first_base_trans["|%s|%s" % (base_parent, base_name)] = base_translation
            first_base_types["|%s|%s" % (base_parent, base_name)] = base_type
            # base_translation=np.array([0,0,0])
            base_type = 5
            break

    keys = first_base_parents.keys()

    x = []
    for line in keys:
        if line.split("|")[2].strip("'").rstrip("'") in x:
            if line.split("|")[2].strip("'").rstrip("'") not in dups:
                dups.append(line.split("|")[2].strip("'").rstrip("'"))
        else:
            x.append(line.split("|")[2].strip("'").rstrip("'"))

    for key in first_base_parents:
        if key.split("|")[2].strip("'").rstrip("'") in dups:
            base_parents[key] = first_base_parents[key]
            base_names.append(key)
        else:
            base_parents[key.split("|")[2].strip("'").rstrip("'")] = first_base_parents[key]
            base_names.append(key.split("|")[2].strip("'").rstrip("'"))
    keys = first_base_trans.keys()

    x = []
    for line in keys:
        if line.split("|")[2].strip("'").rstrip("'") in x:
            if line.split("|")[2].strip("'").rstrip("'") not in dups:
                dups.append(line.split("|")[2].strip("'").rstrip("'"))
        else:
            x.append(line.split("|")[2].strip("'").rstrip("'"))

    for key in first_base_trans:
        if key.split("|")[2].strip("'").rstrip("'") in dups:
            base_trans[key] = first_base_trans[key]
        else:
            base_trans[key.split("|")[2].strip("'").rstrip("'")] = first_base_trans[key]
    keys = first_base_types.keys()
    x = []
    for line in keys:
        if line.split("|")[2].strip("'").rstrip("'") in x:
            if line.split("|")[2].strip("'").rstrip("'") not in dups:
                dups.append(line.split("|")[2].strip("'").rstrip("'"))
        else:
            x.append(line.split("|")[2].strip("'").rstrip("'"))

    for key in first_base_types:
        if key.split("|")[2].strip("'").rstrip("'") in dups:
            base_types[key] = first_base_types[key]
        else:
            base_types[key.split("|")[2].strip("'").rstrip("'")] = first_base_types[key]

    complementary_pair_list = {}
    # abc=(subprocess.Popen('cat %s |grep "lb"| grep "connectAttr"' % "rod3_10_10_10.ma",
    # 	shell = True,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0].split("\n")[:-1])
    new_abc = []
    new_bcd = []
    for i in vHelixFile:
        if re.search('connectAttr.+lb', i):
            new_abc.append(i[:-1])
        if re.search('connectAttr.+\.bw', i):
            new_bcd.append(i)

    for line in new_abc:
        data = line.split(" ")
        base1 = data[1].strip('"').rstrip('.lb"')
        base2 = data[2].strip('"').rstrip('.lb";')
        complementary_pair_list[base1] = base2
        complementary_pair_list[base2] = base1

    neighbour_list = []
    # bcd=(subprocess.Popen('cat %s | grep "backw\|forw" |grep ".bw"| grep "connectAttr"| grep ".fw"' % "rod3_10_10_10.ma",
    # 	shell = True,stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()[0].split("\n")[:-1])
    for line in new_bcd:
        data = line.split(" ")
        base1 = data[1].strip('"').rstrip('.bw"')
        base2 = data[2].strip('"').rstrip('.fw";\n')
        neighbour_list.append([base1, base2])

    # print len(neighbour_list)
    # sys.exit()
    strand_list = []
    counter = 0

    while counter < len(neighbour_list):
        found = "no"
        base1 = neighbour_list[len(neighbour_list) - counter - 1][0]
        base2 = neighbour_list[len(neighbour_list) - counter - 1][1]
        for line in strand_list:
            check = line[-1]
            if check == base1:
                found = "yes"
                line.append(base2)
        if found == "no":
            strand_list.append([base1, base2])
        counter += 1


    moves = 1
    while moves != 0:
        moves = 0
        i = 0
        for line in strand_list:
            i += 1
            j = 0
            check = line[-1]
            check2 = line[-2]
            checkend = line[1]
            for other_line in strand_list:
                j += 1
                other_check = other_line[0]
                other_check2 = other_line[1]
                other_checkend = other_line[-1]
                if check == other_check and check2 == other_check2 and checkend == other_checkend and (i != j):
                    del strand_list[i - 1]
                elif check == other_check and check2 == other_check2 and (i != j):
                    if len(other_line) > len(line):
                        del strand_list[i - 1]
                    else:
                        del strand_list[j - 1]
                    moves += 1
                elif other_check == check and (i != j):
                    for element in other_line[1:]:
                        line.append(element)
                    del strand_list[j - 1]
                    moves += 1

    base_a1s = {}
    base_CoMs = {}
    new_base_trans = {}
    base_array = [1, 0, 3, 2]
    for name_line in base_names:
        if "forw" in name_line:
            for key in complementary_pair_list:
                base1 = key
                base2 = complementary_pair_list[key]
                if base1 == name_line:
                    partner_base = base2
                elif base2 == name_line:
                    partner_base = base1
            # if base1!=base2:
            try:
                own_index = int(base_names.index(name_line))
            except:
                print(name_line)
                print('error at 265')
                sys.exit()
            try:
                own_translation = base_trans[name_line]
            except:
                print(own_index)
                print("failing translation")
                print(name_line)
                sys.exit()
            own_base_type = base_types[name_line]

            own_parent = base_parents[name_line]
            own_parent_index = vHelix_names.index(own_parent)
            own_parent_translation = helix_trans[own_parent]
            own_parent_rotation = helix_rot[own_parent]
            try:
                partner_index = int(base_names.index(partner_base))
            except:
                print(partner_base)
                print('error at 284')
                sys.exit()
            try:
                partner_translation = base_trans[partner_base]
            except:
                print(partner_index)

                sys.exit()
            partner_base_type = base_types[partner_base]

            partner_parent = base_parents[partner_base]
            partner_parent_index = vHelix_names.index(partner_parent)
            partner_parent_translation = helix_trans[partner_parent]
            partner_parent_rotation = helix_rot[partner_parent]

            link_versor = (partner_translation - own_translation) / np.linalg.norm(
                partner_translation - own_translation)

            perp_versor = np.array([-link_versor[1], link_versor[0], 0.])
            unrotated_CoM_element = (own_translation + (0.8494 * link_versor) - (0.1883 * perp_versor))
            new_forw_pos = unrotated_CoM_element - ((0.6 * 0.8518) * link_versor)
            new_backw_pos = unrotated_CoM_element + ((0.6 * 0.8518) * link_versor)
            rot_mat = GetMatrix(own_parent_rotation[0], own_parent_rotation[1], own_parent_rotation[2])
            forw_CoM_element = np.dot(rot_mat, unrotated_CoM_element) + own_parent_translation

            rotated_forw_pos = np.dot(rot_mat, new_forw_pos)
            back_rot_mat = GetMatrix(partner_parent_rotation[0], partner_parent_rotation[1], partner_parent_rotation[2])
            backw_CoM_element = np.dot(back_rot_mat, unrotated_CoM_element) + partner_parent_translation
            rotated_backw_pos = np.dot(back_rot_mat, new_backw_pos)
            base_CoMs[partner_base] = backw_CoM_element
            base_CoMs[name_line] = forw_CoM_element

            correct_forw_pos = (rotated_forw_pos + own_parent_translation) / 0.8518
            correct_backw_pos = (rotated_backw_pos + partner_parent_translation) / 0.8518
            if (correct_forw_pos[0] == correct_backw_pos[0]) and (correct_forw_pos[1] == correct_backw_pos[1]) and (
                correct_forw_pos[2] == correct_backw_pos[2]):
                print("error in translation")
                print(partner_base, name_line, correct_forw_pos)
                sys.exit()
            a1 = (correct_backw_pos - correct_forw_pos) / np.linalg.norm(correct_backw_pos - correct_forw_pos)
            back_a1 = -a1
            base_a1s[partner_base] = back_a1
            base_a1s[name_line] = a1
            new_base_trans[name_line] = correct_forw_pos
            new_base_trans[partner_base] = correct_backw_pos
            if (new_base_trans[name_line][0] == new_base_trans[partner_base][0]) and (
                new_base_trans[name_line][1] == new_base_trans[partner_base][1]) and (
                new_base_trans[name_line][2] == new_base_trans[partner_base][2]):
                print("error in storing translation")
                print(partner_base, name_line, correct_forw_pos)
                sys.exit()
            if partner_base_type == 5 and own_base_type == 5:
                base_types[name_line] = 0
                base_types[partner_base] = 1
            elif partner_base_type == 5:
                comp_base_type = base_array[own_base_type]
                base_types[partner_base] = comp_base_type
            elif own_base_type == 5:
                comp_base_type = base_array[partner_base_type]
                base_types[name_line] = comp_base_type
            if (new_forw_pos[0] == 0 and new_forw_pos[1] == 0 and new_forw_pos[2] == 0) or (
                        new_backw_pos[0] == 0 and new_backw_pos[1] == 0 and new_backw_pos[2] == 0):
                print("total_is_zero", name_line, partner_base)


    i = -1
    for line in base_names:
        i += 1
        if line not in base_a1s.keys():
            base_a1s[line] = [0, 1, 0]
            base_pos = base_trans[line]
            base_CoMs[line] = base_pos
            parent = base_parents[line]
            parent_index = vHelix_names.index(parent)
            parent_translation = helix_trans[parent]
            parent_rotation = helix_rot[parent]
            rot_mat = GetMatrix(parent_rotation[0], parent_rotation[1], parent_rotation[2])
            new_base_trans[line] = ((np.dot(rot_mat, base_pos) + parent_translation) / 0.8518)
            if base_types[line] == 5:
                base_types[line] = 0
                # print new_base_trans[line]
    # sys.exit()
    base_numbers = []
    a1 = []
    pos = []
    a3 = []
    base_letters = []
    base_letter_array = ["A", "T", "C", "G"]
    strand_number = 0
    strand_length = 0
    nucs_to_current_strand = 0
    topology = []
    temp_strand = []
    base_a3s = {}
    # print strand_list[0]
    # print strand_list[0]
    for strand in strand_list:
        # print strand
        '''for element in reversed(strand):
            temp_strand.append(element)
        strand=temp_strand
        temp_strand=[]'''

        nucs_to_current_strand += strand_length
        # print nucs_to_current_strand, strand_number+1
        strand_number += 1
        first_comp_pos = base_CoMs[strand[1]]

        upper_neighbour_list = []
        lower_neighbour_list = []
        strand_length = len(strand)
        circular = 0
        if strand[0] == strand[
            -1]:  # find out if the strand is circlar i.e scaffold and remove last base (duplicate with first base)
            circular = 1
            temp_strand = strand[:-1]
            strand = temp_strand
        if circular == 1 and len(strand) > 1:
            strand_length = len(strand)
            lower_neighbour_list.append(nucs_to_current_strand - 1 + strand_length)
        else:
            lower_neighbour_list.append(-1)
        i = 0
        while i < strand_length - 1:
            lower_neighbour_list.append(i + nucs_to_current_strand)
            upper_neighbour_list.append(i + nucs_to_current_strand + 1)
            i += 1
        if circular == 1 and len(strand) > 1:
            upper_neighbour_list.append(nucs_to_current_strand)
        else:
            upper_neighbour_list.append(-1)
        base_numbers = []
        base_letters = []
        strand_a3 = []
        first_CoM = base_CoMs[strand[0]]
        if strand[1] not in base_a3s.keys():
            strand_a3.append((first_comp_pos - first_CoM) / np.linalg.norm(first_comp_pos - first_CoM))
        else:
            strand_a3.append(base_a3s[strand[1]])

        if first_comp_pos[0] == first_CoM[0] and first_comp_pos[1] == first_CoM[1] and first_comp_pos[2] == first_CoM[
            2]:
            print("a1 incalculable strand_end")
            sys.exit()
        i = 1
        for base in strand[1:]:
            if base not in base_a3s.keys():

                current_pos = base_CoMs[base]
                previous_pos = base_CoMs[strand[i - 1]]
                strand_a3.append((current_pos - previous_pos) / np.linalg.norm(current_pos - previous_pos))
                if base in complementary_pair_list:
                    complementary_base = complementary_pair_list[base]
                else:
                    complementary_base = base
                base_a3s[complementary_base] = (0 - (current_pos - previous_pos)) / np.linalg.norm(
                    current_pos - previous_pos)
            else:
                strand_a3.append(base_a3s[base])
            if (current_pos[0] == previous_pos[0]) and (current_pos[1] == previous_pos[1]) and (
                current_pos[2] == previous_pos[2]):
                print("a1 incalculable")
                print(base)
                print(strand[i - 1])
                # sys.exit()
            i += 1
        i = 0
        for base in strand:
            i += 1
            base_index = base_names.index(base)
            base_numbers.append(base_types[base])
            a1.append(base_a1s[base])
            pos.append(new_base_trans[base])
        # print len(strand), len(strand_a3)
        for line in strand_a3:
            a3.append(line)
        for line in base_numbers:
            index = int(line)
            try:
                base_letters.append(base_letter_array[index])
            except:
                print(index)
        for base, low, up in zip(base_letters, lower_neighbour_list, upper_neighbour_list):
            topology.append("%s %s %s %s" % (strand_number, base, low, up))
            # print len (base_letters), len(lower_neighbour_list), len(upper_neighbour_list), len (strand_a3)

    # Set condition to true if you want to keep the .ma file's name for the .conf and .top files.
    if True:
        if ma_file.endswith('.ma'):
            ma_file = ma_file[:-3]
            conf_file = ma_file + '.conf'
            top_file = ma_file + '.top'
            end_base_file = ma_file + 'end_bases.txt'
    else:
        conf_file = ma_file + 'prova.conf'
        top_file = ma_file + 'prova.top'
        end_base_file = ma_file + 'end_bases.txt'

    if os.path.isfile(conf_file):
        os.remove(conf_file)

    if os.path.isfile(top_file):
        os.remove(top_file)

    if os.path.isfile(end_base_file):
        os.remove(end_base_file)
    topology_file = open(top_file, "w")
    topology_file.close()

    configuration_file = open(conf_file, "w")
    configuration_file.close()
    total_nucs = len(base_names)
    total_strands = (strand_number)

    string = "%s %s \n" % (total_nucs, total_strands)

    with open(top_file, "a+") as top:
        top.write("%s" % string)
        for line in topology:
            top.write("%s \n" % line)

    with open(conf_file, "a+") as conf:
        conf.write("t = 0 \nb = %s %s %s \nE = 0 0 0 \n" % (box_size, box_size, box_size))

        for nuc, versor, a3 in zip(pos, a1, a3):
            i += 1
            pos = "%s %s %s" % (nuc[0], nuc[1], nuc[2])
            versora = "%s %s %s" % (versor[0], versor[1], versor[2])
            a3_versor = "%s %s %s" % (a3[0], a3[1], a3[2])
            string = "%s %s %s 0 0 0 0 0 0" % (pos, versora, a3_versor)
            conf.write("%s \n" % string)

    print("written conf_file and top_file")


def parse_options():
    shortArgs = 'b:'
    longArgs = ['box=']


    opts = {
        "box": 100.,

    }

    try:
        import getopt
        args, positional_args = getopt.gnu_getopt(sys.argv[1:], shortArgs, longArgs)
        for k in args:
            if k[0] == '-b' or k[0] == '--box':
                try:
                    opts['box'] = float(k[1])
                    print >> sys.stderr, "## Setting the box size to %f" % opts['box']
                except ValueError:
                    print >> sys.stderr, "The argument of '%s' should be a number (got '%s' instead)" % (k[0], k[1])
                    exit(1)


    except Exception:
        print_usage()

    return opts

def print_usage():
    print >> sys.stderr, "USAGE:"
    print >> sys.stderr, "\t%s vHelix file in Maya .Ma format" % sys.argv[0]
    print >> sys.stderr, "\t[-b\--box=100]"
    exit(1)



if __name__ == '__main__':




    if len(sys.argv) < 1:
        print_usage()


    source_file = sys.argv[1]





    opts = parse_options()

    box_size = opts['box']
    export_oxDNA(source_file,box_size)
