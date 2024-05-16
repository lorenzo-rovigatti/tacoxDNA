#!/usr/bin/env python3

import os
import pickle
import re
import sys

import numpy as np

from .libs import base
from .libs import cadnano_utils as cu
from .libs import utils

DEBUG = 0
DIST_HEXAGONAL = 2.55  # distance between centres of virtual helices (hexagonal array)
DIST_SQUARE = 2.60  # distance between centres of virtual helices (square array)
BOX_FACTOR = 2  # factor by which to expand the box (linear dimension)


class vh_nodes(object):
    def __init__(self):
        self.begin = []
        self.end = []

    def __str__(self):
        return str([self.begin, self.end])

    def add_begin(self, begin_index):
        if begin_index not in self.begin:
            self.begin.append(begin_index)

    def add_end(self, end_index):
        if end_index not in self.end:
            self.end.append(end_index)


def vhelix_rotation_origami_sq(direction, perp):
    R = utils.get_rotation_matrix(direction, np.pi * 15. / 180)
    return np.dot(R, perp)


def vhelix_rotation_origami_he(direction, perp):
    R = utils.get_rotation_matrix(direction, np.pi * 160. / 180)
    return np.dot(R, perp)


def insert_loop_skip(strands, start_pos, direction, perp, rot, helix_angles, vhelix, nodes, use_seq, seqs):
    # return a double strand which is a copy of the double strand in the first argument, but with skips and loops

    # strand is generated right to left i.e. opposite direction to even vhelix
    length_change = []
    length_change_total = 0
    new_nodes = vh_nodes()
    new_angle = []
    helix_angles_new = np.copy(helix_angles)

    if vhelix.num % 2 == 1:
        reverse_nodes = vh_nodes()
        reverse_nodes.begin = list(reversed(nodes.begin))
        reverse_nodes.end = list(reversed(nodes.end))

    for i in range(len(nodes.begin)):
        # ltr: left to right; looking at the strand left to right (low to high square index), the beginning/end of the effective strand is here (before skips/loops)
        # gs: generated strand; the index of the nucleotide (BEFORE skips/loops are applied) on the generated strand corresponding to the beginning/end of the effective strand
        if vhelix.num % 2 == 0:
            begin_ltr = nodes.begin[i]
            end_ltr = nodes.end[i]
            begin_gs = nodes.begin[i]
            end_gs = nodes.end[i]
        else:
            begin_ltr = reverse_nodes.end[i]
            end_ltr = reverse_nodes.begin[i]
            begin_gs = vhelix.len - reverse_nodes.begin[i] - 1
            end_gs = vhelix.len - reverse_nodes.end[i] - 1

        # check for zero length effective strand
        if end_gs - begin_gs != 0:
            # get length change for this effective strand
            length_change.append(0)
            for j in vhelix.skip[begin_ltr:end_ltr + 1]:
                length_change[i] -= int(j)
            for j in vhelix.loop[begin_ltr:end_ltr + 1]:
                length_change[i] += int(j)
            # get new pitch angles for this effective strand
            new_angle.append(sum(helix_angles[begin_gs:end_gs]) / (end_gs - begin_gs + length_change[i]))
            helix_angles_new[begin_gs:end_gs] = new_angle[i]
            # adjust beginning/end indices according to length change
            begin_gs += length_change_total
            end_gs += length_change_total + length_change[i]
            new_nodes.add_begin(begin_gs)
            new_nodes.add_end(end_gs)  # begin_gs > end_gs.....
        else:
            length_change.append(0)
            new_angle.append(sum(helix_angles) / len(helix_angles))  # append an average angle
            # adjust beginning/end indices according to length change
            begin_gs += length_change_total
            end_gs += length_change_total + length_change[i]
            new_nodes.add_begin(begin_gs)
            new_nodes.add_end(end_gs)
        length_change_total += length_change[i]

    # adjust the new helix angle array according to skips/loops
    deleted = 0
    inserted = 0
    deleted_this_iteration = 0
    inserted_this_iteration = 0
    for i in range(len(nodes.begin)):
        deleted += deleted_this_iteration
        inserted += inserted_this_iteration
        deleted_this_iteration = 0
        inserted_this_iteration = 0
        if vhelix.num % 2 == 0:
            begin_ltr = nodes.begin[i]
            end_ltr = nodes.end[i]
            begin_gs = nodes.begin[i]
            end_gs = nodes.end[i]
        else:
            begin_ltr = reverse_nodes.end[i]
            end_ltr = reverse_nodes.begin[i]
            begin_gs = vhelix.len - reverse_nodes.begin[i] - 1
            end_gs = vhelix.len - reverse_nodes.end[i] - 1
        for j in vhelix.skip[begin_ltr:end_ltr + 1]:
            if j == 1:
                helix_angles_new = np.delete(helix_angles_new, begin_gs - deleted + inserted)
                deleted_this_iteration += 1
        for j in vhelix.loop[begin_ltr:end_ltr + 1]:
            for _ in range(j):
                helix_angles_new = np.insert(helix_angles_new, begin_gs - deleted + inserted, new_angle[i])
                inserted_this_iteration += 1
                
    g = cu.StrandGenerator()
    new_strands = g.generate_or_sq(len(helix_angles_new) + 1, start_pos=start_pos, direction=direction, perp=perp, double=True, rot=rot, angle=helix_angles_new, length_change=length_change, region_begin=new_nodes.begin, region_end=new_nodes.end)
    if use_seq:
        try:
            sequence = [x for x in seqs[vhelix.cad_index]]
        except IndexError:
            base.Logger.die("sequence file contains too few rows compared to the number of virtual helices in the cadnano file, dying")
        if vhelix.num % 2 == 1:
            sequence.reverse()
        if new_strands[0].get_length() != len(sequence):
            base.Logger.log("Cannot change sequence: lengths don't match; virtual helix %s, sequence length %s, virtual helix length %s - are skips/loops accounted for?" % (vhelix.num, len(sequence), new_strands[0].get_length()), base.Logger.WARNING)
        else:
            new_strands[0].set_sequence(sequence)
        sequence2 = [3 - s for s in sequence]
        sequence2.reverse()
        if new_strands[0].get_length() != len(sequence):
            base.Logger.log("Cannot change sequence: lengths don't match; virtual helix %s, sequence length %s, virtual helix length %s - are skips/loops accounted for?" % (vhelix.num, len(sequence), new_strands[0].get_length()), base.Logger.WARNING)
        else:
            new_strands[1].set_sequence(sequence2)

    return new_strands


def add_slice(current_system, vhelix, begin, end, nodes, strands, pos, direction, perp, rot, helix_angles, strand_type, use_seq, seqs):
    # add a slice of the virtual helix to the slice system, taking into account skips and loops
    length_change_begin = 0
    length_change_end = 0
    
    if (vhelix.num % 2 + strand_type) % 2 == 0:  # strand and even num or staple and odd num
        for i in vhelix.skip[:begin]:
            length_change_begin -= int(i)
        for i in vhelix.skip[:end + 1]:
            length_change_end -= int(i)
            
        for i in vhelix.loop[:begin]:
            length_change_begin += int(i)
        for i in vhelix.loop[:end + 1]:
            length_change_end += int(i)
    
        begin_slice = begin + length_change_begin
        end_slice = end + 1 + length_change_end

    else:
        for i in vhelix.skip[end:]:
            length_change_end -= int(i)
        for i in vhelix.skip[begin + 1:]:
            length_change_begin -= int(i)
            
        for i in vhelix.loop[end:]:
            length_change_end += int(i)
        for i in vhelix.loop[begin + 1:]:
            length_change_begin += int(i)

        begin_slice = vhelix.len - begin - 1 + length_change_begin
        end_slice = vhelix.len - end + length_change_end

    new_strands = insert_loop_skip(strands, pos, direction, perp, rot, helix_angles, vhelix, nodes, use_seq, seqs)
    current_system.add_strand(new_strands[strand_type].get_slice(begin_slice, end_slice), check_overlap=False)
    return current_system


def add_slice_nupack(vhelix, strand_number, begin_helix, end_helix, index_lookup, strand_type):
    length_change = 0
    skips = 0
    loops = 0
    if (vhelix.num % 2 + strand_type) % 2 == 0 :  # strand and even num or staple and odd num
        for i in vhelix.skip[begin_helix:end_helix + 1]:
            length_change -= int(i)
            
        for i in vhelix.loop[begin_helix:end_helix + 1]:
            length_change += int(i)
    
    else:
        for i in vhelix.skip[end_helix:begin_helix + 1]:
            length_change -= int(i)
            
        for i in vhelix.loop[end_helix:begin_helix + 1]:
            length_change += int(i)
            
    if (strand_type + vhelix.num % 2) % 2 == 0:
        iter_length = end_helix - begin_helix + 1 + length_change
    else:
        iter_length = begin_helix + 1 - end_helix + length_change

    nucleotide = 0
    while nucleotide < iter_length:
        if (strand_type + vhelix.num % 2) % 2 == 0:
            vhelix_base = nucleotide + begin_helix + skips - loops
        else:
            vhelix_base = begin_helix - nucleotide - skips + loops
        if vhelix.skip[vhelix_base] != 1:
            if (strand_type + vhelix.num % 2) % 2 == 0:
                add_nuc = [nucleotide + x for x in range(vhelix.loop[vhelix_base] + 1)]
            else:
                add_nuc = [nucleotide + x for x in range(vhelix.loop[vhelix_base] + 1)[::-1]]
            if strand_type == 0:
                index_lookup[(vhelix.num, vhelix_base)] = [strand_number, nucleotide]
            elif strand_type == 1:
                index_lookup[(strand_number, nucleotide)] = [vhelix.num, vhelix_base]
            elif strand_type == 2:
                index_lookup.add_scaf(vhelix.num, vhelix_base, strand_number, add_nuc)
            elif strand_type == 3:
                index_lookup.add_stap(vhelix.num, vhelix_base, strand_number, add_nuc)
            nucleotide += 1 + vhelix.loop[vhelix_base]
            loops += vhelix.loop[vhelix_base]
        else:
            if strand_type == 2:
                index_lookup.add_scaf(vhelix.num, vhelix_base, strand_number, [])
            elif strand_type == 3:
                index_lookup.add_stap(vhelix.num, vhelix_base, strand_number, [])
            skips += 1

    return index_lookup


def build_nodes(vh):
    # returns a vh_nodes object which contains the beginning and end square indices of each effective strand. Effective strands
    # on a given vhelix are a combination of information on staple and scaffold strands. Together they tell us which
    # nucleotides need to be held constant when we alter the angles/base-base distances for a section of nucleotides along a final strand
    nodes = vh_nodes()
    if vh.num % 2 == 0:
        direction = 1
    else:
        direction = -1
    for i in range(len(vh.scaf)):
        # need to consider what happens when I add an index to the node list that doesn't fall within the range of square indices in the vhelix
        previd = i - 1 * direction
        nextid = i + 1 * direction
        if previd in range(len(vh.scaf)):
            prev = vh.scaf[previd].type(vh, previd)
            prev_stap = vh.stap[previd].type(vh, previd)
        else:
            prev = False
            prev_stap = False
        if nextid in range(len(vh.scaf)):
            next_ = vh.scaf[nextid].type(vh, nextid)
            next_stap = vh.stap[nextid].type(vh, nextid)
        else:
            next_ = False
            next_stap = False
        # now build the effective strand vh_nodes object
        if not ((prev == prev_stap and (prev == 'begin' or prev == 'end')) or (next_ == next_stap and (next_ == 'begin' or next_ == 'end'))):
            if vh.scaf[i].type(vh, i) == 'empty':
                if vh.stap[i].type(vh, i) == 'begin':
                    nodes.add_end(i)
                elif vh.stap[i].type(vh, i) == 'end':
                    nodes.add_begin(i)
            elif vh.scaf[i].type(vh, i) == 'begin':
                if vh.stap[i].type(vh, i) == 'empty':
                    nodes.add_begin(i)
                elif vh.stap[i].type(vh, i) == 'continue':
                    nodes.add_begin(i)
                    nodes.add_end(i - 1 * direction)
                elif vh.stap[i].type(vh, i) == 'begin':
                    nodes.add_begin(i + 1 * direction)
                    nodes.add_end(i - 1 * direction)
                elif vh.stap[i].type(vh, i) == 'end':
                    nodes.add_begin(i)
            elif vh.scaf[i].type(vh, i) == 'end':
                if vh.stap[i].type(vh, i) == 'empty':
                    nodes.add_end(i)
                elif vh.stap[i].type(vh, i) == 'continue':
                    nodes.add_begin(i + 1 * direction)
                    nodes.add_end(i)
                elif vh.stap[i].type(vh, i) == 'begin':
                    nodes.add_end(i)
                elif vh.stap[i].type(vh, i) == 'end':
                    nodes.add_begin(i + 1 * direction)
                    nodes.add_end(i - 1 * direction)
            elif vh.scaf[i].type(vh, i) == 'continue':
                if vh.stap[i].type(vh, i) == 'begin':
                    nodes.add_begin(i + 1 * direction)
                    nodes.add_end(i)
                elif vh.stap[i].type(vh, i) == 'end':
                    nodes.add_begin(i)
                    nodes.add_end(i - 1 * direction)
                    
    return nodes


def generate_vhelices_origami_sq(vhelix_direction, vhelix_perp, h):
    g = cu.StrandGenerator()
    # generate helix angles
    helix_angles = np.zeros(h.len - 1, dtype=float)
    # hard upper limit on pitch angle seems to be between 54.5 and 55 degrees
    for i in range(len(helix_angles)):
        modi = i % 32
        if modi < 2:
            helix_angles[i] = 28 * np.pi / 180
        elif modi == 2:
            helix_angles[i] = 36 * np.pi / 180
        elif modi == 3:
            helix_angles[i] = 54.375 * np.pi / 180
        elif modi == 4:
            helix_angles[i] = 37 * np.pi / 180
        elif modi in (5, 6):
            helix_angles[i] = 27.6666666666666 * np.pi / 180
        elif modi == 7:
            helix_angles[i] = 30.6666666666666 * np.pi / 180
        elif modi in (8, 9):
            helix_angles[i] = 29.3333333333 * np.pi / 180
        elif modi == 10:
            helix_angles[i] = 34.3333333333 * np.pi / 180
        elif modi == 11:
            helix_angles[i] = 54.5 * np.pi / 180
        elif modi in (12, 13):
            helix_angles[i] = (28.91666666666 * np.pi / 180)
        elif modi in (14, 15, 16, 17):
            helix_angles[i] = 31.16666666666 * np.pi / 180
        elif modi == 18:
            helix_angles[i] = 35.5 * np.pi / 180
        elif modi == 19:
            helix_angles[i] = 52 * np.pi / 180
        elif modi == 20:
            helix_angles[i] = 35.5 * np.pi / 180
        elif modi in (21, 22):
            helix_angles[i] = 27.5 * np.pi / 180
        elif modi == 23:
            helix_angles[i] = 35.5 * np.pi / 180
        elif modi >= 24 and modi < 27:
            helix_angles[i] = 30 * np.pi / 180
        elif modi == 27:
            helix_angles[i] = 52 * np.pi / 180
        elif modi == 28:
            helix_angles[i] = 35.5 * np.pi / 180
        else:
            helix_angles[i] = 30.91666666666 * (np.pi / 180)

    # make sure the helices are periodic in 32 bases
    total_sum = 0
    for i in range(31):
        total_sum += helix_angles[i]

    for i in range(len(helix_angles)):
        if i % 32 == 31:
            helix_angles[i] = 1080 * np.pi / 180 - total_sum
            
    # make the virtual helices
    if h.num % 2 == 0:
        pos = np.array([h.col * DIST_SQUARE, h.row * DIST_SQUARE, 0])
        direction = vhelix_direction
        perp = vhelix_perp
        rot = 0.
        angles = helix_angles
        strands = g.generate_or_sq(h.len, start_pos=pos, direction=direction, perp=perp, double=True, rot=rot, angle=angles)

    else:
        pos = np.array([h.col * DIST_SQUARE, h.row * DIST_SQUARE, (h.len - 1) * base.BASE_BASE])
        direction = -vhelix_direction
        perp = -vhelix_perp
        rot = -np.sum(helix_angles) % (2 * np.pi)
        angles = np.flipud(helix_angles)
        strands = g.generate_or_sq(h.len, start_pos=pos, direction=direction, perp=perp, double=True, rot=rot, angle=angles)

    return (strands[0], strands[1]), helix_angles, pos, rot, direction, perp


def generate_vhelices_origami_he(vhelix_direction, vhelix_perp, h):
    g = cu.StrandGenerator()
    # generate helix angles
    helix_angles = np.zeros(h.len - 1, dtype=float)

    for i in range(len(helix_angles)):
        modi = i % 21
        if modi == 0:
            helix_angles[i] = 32.571 * np.pi / 180
        elif modi == 1:
            helix_angles[i] = 36 * np.pi / 180
        elif modi in (1, 2, 3):
            helix_angles[i] = 42 * np.pi / 180
        elif modi in (5, 6, 7):
            helix_angles[i] = 29.143 * np.pi / 180
        elif modi == 8:
            helix_angles[i] = 32 * np.pi / 180
        elif modi in (9, 10):
            helix_angles[i] = 44 * np.pi / 180
        elif modi in (12, 13, 14):
            helix_angles[i] = 28.571 * np.pi / 180
        elif modi in (16, 17):
            helix_angles[i] = 41.5 * np.pi / 180
        elif modi in (19, 20):
            helix_angles[i] = 28.476 * np.pi / 180
        else:
            helix_angles[i] = 720. / 21 * (np.pi / 180.) 

    # make sure it's periodic
    total_sum = 0
    for i in range(20):
        total_sum += helix_angles[i]

    for i in range(len(helix_angles)):
        if i % 21 == 20:
            helix_angles[i] = 720. * np.pi / 180 - total_sum

    # make the virtual helices
    if h.num % 2 == 0:
        pos = np.array([h.col * np.sqrt(3) * DIST_HEXAGONAL / 2, h.row * 3 * DIST_HEXAGONAL / 2, 0])
        direction = vhelix_direction
        perp = vhelix_perp
        rot = 0.
        strands = g.generate_or_sq(h.len, start_pos=pos, direction=direction, perp=perp, double=True, rot=rot, angle=helix_angles)

    else:
        pos = np.array([h.col * np.sqrt(3) * DIST_HEXAGONAL / 2, h.row * 3 * DIST_HEXAGONAL / 2 + DIST_HEXAGONAL / 2, (h.len - 1) * base.BASE_BASE])
        direction = -vhelix_direction
        perp = -vhelix_perp
        if base.MM_GROOVING:
            rot = -np.sum(helix_angles) % (2 * np.pi) - 0.07
        else:
            rot = -np.sum(helix_angles) % (2 * np.pi)
        angles = np.flipud(helix_angles)
        strands = g.generate_or_sq(h.len, start_pos=pos, direction=direction, perp=perp, double=True, rot=rot, angle=angles)

    return (strands[0], strands[1]), helix_angles, pos, rot, direction, perp


# cadnano object structure
class vstrands (object):

    def __init__(self):
        self.vhelices = []

    def add_vhelix(self, toadd):
        self.vhelices.append(toadd)

    def bbox(self):
        rows = []
        cols = []
        lens = []
        for h in self.vhelices:
            rows.append(h.row)
            cols.append(h.col)
            lens.append(len(h.stap))

        dr = DIST_SQUARE * (max(rows) - min(rows) + 2)
        dc = DIST_SQUARE * (max(cols) - min(cols) + 2)
        dl = 0.34 * (max(lens) + 2)
        
        return 2 * max([dr, dc, dl]) * BOX_FACTOR
    
    def __str__(self):
        a = '{\n"vstrands":[\n'
        if len(self.vhelices) > 0:
            for h in self.vhelices:
                a = a + str(h) + ','
            a = a[0:len(a) - 1]
        a = a + '}\n'
        return a


class vhelix (object):

    def __init__(self):
        self.stapLoop = []
        self.scafLoop = []
        self.skip = []
        self.loop = []
        self.stap_colors = []
        self.row = 0
        self.col = 0
        self.num = 0
        self.stap = []
        self.scaf = []
        self.cad_index = -1
        self.skiploop_bases = 0

    def get_length(self):
        return max (len(self.scaf), len(self.stap))

    len = property (get_length)

    def add_square(self, toadd, which):
        if which == 'stap':
            self.stap.append(toadd)
        elif which == 'scaf':
            self.scaf.append (toadd)
        else:
            base.Logger.log("Cannot add square that is not scaf or stap. Dying now", base.Logger.CRITICAL)
            sys.exit(1)
    
    def __str__(self):
        a = '{\n'

        a = a + '"stapLoop":['
        if len(self.stapLoop) > 0:
            for i in self.stapLoop:
                a = a + str(i) + ','
            a = a[0:len(a) - 1]  # remove last comma
        a = a + '],\n'

        a = a + '"skip":['
        if len(self.skip) > 0:
            for e in self.skip:
                a = a + str(e) + ','
            a = a[0:len(a) - 1]  # remove last comma
        a = a + '],\n'
        
        a = a + '"loop":['
        if len(self.loop) > 0:
            for e in self.loop:
                a = a + str(e) + ','
            a = a[0:len(a) - 1]  # remove last comma
        a = a + '],\n'
        
        a = a + '"stap_colors":['
        if len (self.stap_colors) > 0:
            for e in self.stap_colors:
                a = a + str(e) + ','
            a = a[0:len(a) - 1]  # remove last comma
        a = a + '],\n'

        a = a + '"row":' + str(self.row) + ',\n'
        a = a + '"col":' + str(self.col) + ',\n'
        a = a + '"num":' + str(self.num) + ',\n'
        
        a = a + '"scafLoop":['
        if len(self.scafLoop) > 0:
            for i in self.scafLoop:
                a = a + str(i) + ','
            a = a[0:len(a) - 1]  # remove last comma
        a = a + '],\n'
        
        a = a + '"stap":['
        if len(self.stap) > 0:
            for i in self.stap:
                a = a + str(i) + ','
            a = a[0:len(a) - 1]  # remove last comma
        a = a + '],\n'
        
        a = a + '"scaf":['
        if len(self.scaf) > 0:
            for i in self.scaf:
                a = a + str(i) + ','
            a = a[0:len(a) - 1]  # remove last comma
        a = a + ']\n}'
        return a


class square(object):

    def __init__ (self, V_0=-1, b_0=-1, V_1=-1, b_1=-1):
        """
        V_0, b_0, V_1, b_1 are integer indices correspond to:
        virtual_helix_behind, virtual_base_behind, virtual_helix_ahead, virtual_base_ahead
        """
        self.V_0 = V_0
        self.b_0 = b_0
        self.V_1 = V_1
        self.b_1 = b_1

    def __str__ (self):
        return '[%i,%i,%i,%i]' % (self.V_0, self.b_0, self.V_1, self.b_1)

    def type(self, vhelix, myid):
        # find type of strand (junction) on this square
        # currently direction always equals zero...
        direction = 0
        if self.V_0 == -1 and self.b_0 == -1:
            if self.V_1 == -1 and self.b_1 == -1:
                return 'empty'
            elif self.V_1 == vhelix.num and abs(self.b_1 - myid) == 1:
                if direction == 0:
                    return 'begin'
                else:
                    return 'end'
        elif self.V_0 == vhelix.num and abs(self.b_0 - myid) == 1:
            if self.V_1 == -1:
                if direction == 0:
                    return 'end'
                else:
                    return 'begin'
            elif self.V_1 == vhelix.num and abs(self.b_1 - myid) == 1:
                return 'continue'
            else:
                # join
                if direction == 0:
                    return 'end'
                else:
                    return 'begin'
        else:
            if self.V_1 == vhelix.num and abs(self.b_1 - myid) == 1:
                if direction == 0:
                    return 'begin'
                else:
                    return 'end'

        # shouldn't get to here
        base.Logger.log('unexpected square array', base.Logger.WARNING)

        
def parse_cadnano(path):
    import json
    
    
    try:
        with open(path) as json_data:
            cadnano = json.load(json_data)

            cadsys = vstrands()
            for vstrand in cadnano["vstrands"]:
                vh = vhelix()
                for key, val in list(vstrand.items()):
                    if key == "skip":
                        vh.skip = [abs(int(x)) for x in val]
                    else:
                        setattr(vh, key, val)
                vh.stap = [square(*i) for i in vh.stap]
                vh.scaf = [square(*i) for i in vh.scaf]
                vh.skiploop_bases = len(vh.skip) + sum(vh.loop) - sum(vh.skip)
                cadsys.add_vhelix(vh)
    except IOError:
        print("File '" + path + "' not found, aborting", file=sys.stderr)
        sys.exit(1)
    except ValueError:
        print("Invalid json file '" + path + "', aborting", file=sys.stderr)
        sys.exit(1)
    except:
        print("Caught an error while parsing '" + path + "', aborting", file=sys.stderr)
        sys.exit(1)
        
    return cadsys


def parse_jsonobject(jsonobject):
    
    try: 
        cadsys = vstrands()    
        for vstrand in jsonobject["vstrands"]:
            vh = vhelix()
            for key, val in list(vstrand.items()):
                if key == "skip":
                    vh.skip = [abs(int(x)) for x in val]
                else:
                    setattr(vh, key, val)
            vh.stap = [square(*i) for i in vh.stap]
            vh.scaf = [square(*i) for i in vh.scaf]
            vh.skiploop_bases = len(vh.skip) + sum(vh.loop) - sum(vh.skip)
            cadsys.add_vhelix(vh)
    except ValueError:
        print("Invalid json '" + jsonobject + "', aborting", file=sys.stderr)
        sys.exit(1)
    except:
        print("Caught an error while converting '" + jsonobject + "', aborting", file=sys.stderr)
        sys.exit(1)
        
    return cadsys


def print_usage():
    print("USAGE:", file=sys.stderr)
    print("\t%s cadnano_file lattice_type" % sys.argv[0], file=sys.stderr)
    print("\t[-q\--sequence FILE] [-b\--box VALUE] [-e\--seed VALUE] [-p\--print-virt2nuc] [-o\--print-oxview]", file=sys.stderr)
    exit(1)

def randomSequenceGenerator(cadsys):
    base.Logger.log("No sequence file given, using random sequence", base.Logger.INFO)
    sequences = []
    for vhelix in (cadsys.vhelices):
        seq = []
        for _ in range(vhelix.skiploop_bases):
            seq.append(np.random.randint(0, 4))
        sequences.append(seq)
    
    return sequences

def readingCli(*args):
    if len(sys.argv) < 3:
        print_usage()
        
    shortArgs = 'q:b:e:po'
    longArgs = ['sequence=', 'box=', 'seed=', 'print_lattice_id_map', 'print-virt2nuc', 'print-oxview']
    
    side = False
    sequence_filename = False
    print_lattice_id_map = False
    print_virt2nuc = False
    print_oxview = False
    source_file = sys.argv[1]
    
    np_seed = None

    if sys.argv[2] == "sq":
        lattice_type = "sq"
    elif sys.argv[2] == "he":
        lattice_type = "he"
    else:
        print("Lattice_type should be either 'sq' or 'he'", file=sys.stderr)
        exit(1)
    
    try:
        import getopt
        args, files = getopt.gnu_getopt(sys.argv[3:], shortArgs, longArgs)
        for k in args:
            if k[0] == '-q' or k[0] == "--sequence": 
                sequence_filename = k[1]
            elif k[0] == '-b' or k[0] == "--box": 
                side = float(k[1])
                base.Logger.log("The system will be put in a box of side %s (in oxDNA simulation units)" % str(side), base.Logger.INFO)
            elif k[0] == '-e' or k[0] == "--seed": 
                np_seed = int(k[1])
                np.random.seed(np_seed)
            elif k[0] == '-m' or k[0] == "--print_lattice_id_map":
                print_lattice_id_map == True
            elif k[0] == '-p' or k[0] == "--print-virt2nuc":
                print_virt2nuc = True
            elif k[0] == '-o' or k[0] == "--print-oxview":
                print_oxview = True
            
            
    except Exception:
        print_usage()

    return source_file, lattice_type, sequence_filename, side, np_seed, print_virt2nuc, print_oxview

def parsingCli(source_file, sequence_filename):
    cadsys = parse_cadnano(source_file)
    base.Logger.log("Using json file %s" % source_file, base.Logger.INFO)

    # define sequences by vhelix
    sequence_file = 0
    sequences = []
    if sequence_filename:
        sequence_file = open(sequence_filename, "r")
        base.Logger.log("Using sequence file '%s'" % sequence_filename, base.Logger.INFO)
        # with this we can remove all whitespace and we don't have issues with the different newline sequences (\n vs \r\n)
        pattern = re.compile('\s+')
        lines = sequence_file.readlines()
        for line in lines:
            seq = []
            for x in re.sub(pattern, '', line):
                if x in ["R", "r"]:
                    seq.append(np.random.randint(0, 4))
                else:
                    try:
                        seq.append(base.base_to_number[x])
                    except KeyError:
                        base.Logger.log("KeyError while converting base names to integer; check the sequence file", base.Logger.CRITICAL)
                        sys.exit()
            sequences.append(seq)
    else:
        sequences = randomSequenceGenerator(cadsys)

    return cadsys, sequences

def cadnano_oxdna(output_file, cadsys, lattice_type, input_sequences=None, side=False, print_lattice_id_map=False, print_virt2nuc=False, print_oxview=False):    
    vh_vb2nuc = cu.vhelix_vbase_to_nucleotide()
    vh_vb2nuc_final = cu.vhelix_vbase_to_nucleotide()

    sequences = randomSequenceGenerator(cadsys) if input_sequences is None else input_sequences

    is_1strand_in_multiple_vhelixes_special_case = len(sequences) == 1 and len(cadsys.vhelices) > 1
    # 1strand_in_multiple_vhelixes_special_case: 1 strand system (i.e. NOT double helix) across many vhelices and defined with 1 .sqs line
    if is_1strand_in_multiple_vhelixes_special_case:
        base.Logger.log("One line detected in the sequence file. Since the cadnano file contains more than 1 virtual helix, the sequence found will be used as we were dealing with a single-strand system", base.Logger.INFO)
        single_strand_system = True
    else:
        single_strand_system = False

    vhelix_counter = 0
    if not side:
        side = cadsys.bbox()
        base.Logger.log("Using default box size, a factor %s larger than the size of the cadnano system" % str(BOX_FACTOR), base.Logger.INFO)
    vhelix_direction_initial = np.array([0., 0., 1.])
    vhelix_perp_initial = np.array([1., 0., 0.])
    if lattice_type == "sq":
        vhelix_perp_initial = vhelix_rotation_origami_sq(vhelix_direction_initial, vhelix_perp_initial)
    elif lattice_type == "he":
        vhelix_perp_initial = vhelix_rotation_origami_he(vhelix_direction_initial, vhelix_perp_initial)

    slice_sys = base.System([side, side, side])
    final_sys = base.System([side, side, side])
    strand_number = -1
    partner_list_scaf = []
    partner_list_stap = []
    found_partner = False
    join_list_scaf = []
    join_list_stap = []
    begin_helix = -1
    end_helix = -1
    for h in cadsys.vhelices:
        h.cad_index = vhelix_counter
        if lattice_type == "sq":
            strands, helix_angles, pos, rot, vhelix_direction, vhelix_perp = generate_vhelices_origami_sq(vhelix_direction_initial, vhelix_perp_initial, h)
        elif lattice_type == "he":
            strands, helix_angles, pos, rot, vhelix_direction, vhelix_perp = generate_vhelices_origami_he(vhelix_direction=vhelix_direction_initial, vhelix_perp=vhelix_perp_initial, h=h)
        else:
            print("Unknown lattice type!\n")
            exit(1)

        nodes = build_nodes(h)
        
        # read the scaffold squares and add strands to slice_sys
        i = 0
        for s in h.scaf:
            if s.V_0 == -1 and s.b_0 == -1:
                if s.V_1 == -1 and s.b_0 == -1:
                    pass
                elif s.V_1 == h.num and abs(s.b_1 - i) == 1:
                    if h.num % 2 == 0:
                        strand_number += 1
                    begin_helix = i
                    if h.num % 2 == 1:
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 0, not single_strand_system, sequences)
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 2)
                else:
                    base.Logger.log("unexpected square array", base.Logger.WARNING)
            elif s.V_0 == h.num and abs(s.b_0 - i) == 1:
                if s.V_1 == -1 and s.b_1 == -1:
                    if h.num % 2 == 1:
                        strand_number += 1
                    end_helix = i
                    if h.num % 2 == 0:
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 0, not single_strand_system, sequences)
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 2)
                elif s.V_1 == h.num and abs(s.b_1 - i) == 1:
                    pass
                else:
                    if h.num % 2 == 1:
                        strand_number += 1
                    end_helix = i
                    if h.num % 2 == 0 :
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 0, not single_strand_system, sequences)
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 2)

                    if h.num % 2 == 1:
                        column = i
                    else:
                        column = i
                    for j in range(len(partner_list_scaf)):

                        if [h.num, column] == partner_list_scaf[j]:
                            join_list_scaf[j].insert(0, strand_number)
                            found_partner = True
                    if found_partner == False:
                        join_list_scaf.append([strand_number])
                        partner_list_scaf.append([s.V_1, s.b_1])
                    found_partner = False
            else:
                if s.V_1 == -1 and s.b_1 == -1:
                    base.Logger.log("unexpected square array", base.Logger.WARNING)
                elif s.V_1 == h.num and abs(s.b_1 - i) == 1:
                    if h.num % 2 == 0:
                        strand_number += 1
                    begin_helix = i
                    if h.num % 2 == 1:
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 0, not single_strand_system, sequences)
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 2)

                    for j in range(len(partner_list_scaf)):
                        if h.num % 2 == 1:
                            column = i
                        else:
                            column = i
                        if [h.num, column] == partner_list_scaf[j]:
                            join_list_scaf[j].append(strand_number)
                            found_partner = True
                    if found_partner == False:
                        join_list_scaf.append([strand_number])
                        partner_list_scaf.append([s.V_0, s.b_0])
                    found_partner = False
                else:
                    base.Logger.log("unexpected square array", base.Logger.WARNING)                
            i += 1
            
        if slice_sys.N_strands == 0:
            base.Logger.log("No scaffold strand found in virtual helix n. %d: staples-only virtual helices are not supported" % h.num, base.Logger.WARNING)
            continue

        # read the staple squares and add strands to slice_sys
        i = 0
        for s in h.stap:
            if s.V_0 == -1 and s.b_0 == -1:
                if s.V_1 == -1 and s.b_0 == -1:
                    pass
                elif s.V_1 == h.num and abs(s.b_1 - i) == 1:
                    if h.num % 2 == 1:
                        strand_number += 1
                    begin_helix = i
                    if h.num % 2 == 0:
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 1, not single_strand_system, sequences)
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 3)
                else:
                    base.Logger.log("unexpected square array", base.Logger.WARNING)
            elif s.V_0 == h.num and abs(s.b_0 - i) == 1:
                if s.V_1 == -1 and s.b_1 == -1:
                    if h.num % 2 == 0:
                        strand_number += 1
                    end_helix = i
                    if h.num % 2 == 1:
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 1, not single_strand_system, sequences)
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 3)

                elif s.V_1 == h.num and abs(s.b_1 - i) == 1:
                    pass
                else:
                    if h.num % 2 == 0:
                        strand_number += 1
                    end_helix = i
                    if h.num % 2 == 1:
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 1, not single_strand_system, sequences)
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 3)

                    if h.num % 2 == 0:
                        column = i
                    else:
                        column = i
                    for j in range(len(partner_list_stap)):

                        if [h.num, column] == partner_list_stap[j]:
                            join_list_stap[j].insert(0, strand_number)
                            found_partner = True
                    if found_partner == False:
                        join_list_stap.append([strand_number])
                        partner_list_stap.append([s.V_1, s.b_1])
                    found_partner = False
            else:
                if s.V_1 == -1 and s.b_1 == -1:
                    base.Logger.log("unexpected square array", base.Logger.WARNING)
                elif s.V_1 == h.num and abs(s.b_1 - i) == 1:
                    if h.num % 2 == 1:
                        strand_number += 1
                    begin_helix = i
                    if h.num % 2 == 0:
                        slice_sys = add_slice(slice_sys, h, begin_helix, end_helix, nodes, strands, pos, vhelix_direction, vhelix_perp, rot, helix_angles, 1, not single_strand_system, sequences)
                        vh_vb2nuc = add_slice_nupack(h, strand_number, begin_helix, end_helix, vh_vb2nuc, 3)

                    for j in range(len(partner_list_stap)):
                        if h.num % 2 == 0:
                            column = i
                        else:
                            column = i
                        if [h.num, column] == partner_list_stap[j]:
                            join_list_stap[j].append(strand_number)
                            found_partner = True
                    if found_partner == False:
                        join_list_stap.append([strand_number])
                        partner_list_stap.append([s.V_0, s.b_0])
                    found_partner = False
                else:
                    base.Logger.log("unexpected square array", base.Logger.WARNING)                
            i += 1
        vhelix_counter += 1

    join_lists = [join_list_scaf, join_list_stap]

    # add strands to final_sys that aren't joined
    join_list_unpacked = []
    for a in range(2):
        for i in join_lists[a]:
            join_list_unpacked.extend(i)
    for i in range(len(slice_sys._strands)):
        if i not in join_list_unpacked:
            final_sys.add_strand(slice_sys._strands[i], check_overlap=False)
            vh_vb2nuc_final.add_strand(i, vh_vb2nuc)

    for a in range(2):
        join_list = join_lists[a]
        all_are_joined = False
        restart = False

        # check distance between the backbones we are about to join
        for pair in join_list:
            strand1 = slice_sys._strands[pair[0]]
            strand2 = slice_sys._strands[pair[1]]
            backbone_backbone_dist = strand1._nucleotides[-1].distance(strand2._nucleotides[0], PBC=False)
            absolute_bb_dist = np.sqrt(np.dot(backbone_backbone_dist, backbone_backbone_dist))
            if absolute_bb_dist > 1.0018 or absolute_bb_dist < 0.5525:
                base.Logger.log("the backbone-backbone distance across joints is %f: it will have to be relaxed with preliminary simulations" % absolute_bb_dist, base.Logger.WARNING)

        # match up all the pairs of joins that involve the same strand
        circular = []
        while all_are_joined == False:
            restart = False
            for i in range(len(join_list)):
                if restart == True:
                    break
                for j in range(len(join_list)):
                    if restart == True:
                        break
                    if join_list[i][0] == join_list[j][-1]:
                        if i != j:
                            join_list[j].extend(join_list[i][1:])
                            join_list.pop(i)
                            restart = True
                            break
                        else:
                            if i not in circular:
                                circular.append(i)

            if restart == False:
                all_are_joined = True

        # add joined strands
        for ii, join in enumerate(join_list):
            joined_strand = slice_sys._strands[join[0]]
            if ii in circular:
                for k in range(1, len(join) - 1):
                    joined_strand = joined_strand.append(slice_sys._strands[join[k]])
                joined_strand.make_circular(check_join_len=True)
            else:
                for k in range(1, len(join)):
                    joined_strand = joined_strand.append(slice_sys._strands[join[k]])
                
            final_sys.add_strand(joined_strand, check_overlap=False)

            # This is a bug fix. Ben 12/2/14
            # for a circular strand we need to terminate the strand one element early (so reduce the length
            # of the range by 1), since the final element is just a repeat of the first one.
            if joined_strand._circular:
                joining_range = list(range(len(join) - 2))
            else:
                joining_range = list(range(len(join) - 1))
            # add joined strands to v2n index
            for k in joining_range:
                vh_vb2nuc_final.add_strand(join[k], vh_vb2nuc, continue_join=True)
            vh_vb2nuc_final.add_strand(join[k + 1], vh_vb2nuc, continue_join=False)
                
            if single_strand_system:
                final_sys._strands[0].set_sequence(sequences[0])
    
    if input_sequences is not None and single_strand_system and len(final_sys._strands) > 1:
        base.Logger.log("more than one strand detected - sequence file will not be read", base.Logger.WARNING)
        final_sys._strands[0].set_sequence(np.random.randint(0, 4, len(final_sys._strands[0]._nucleotides)))  # this line does not work

    # Fix to reverse the direction of every strand so that the 3' to 5' direction is the same
    # as in Cadnano. In cadnano the strands point in the 5' to 3' direction, whereas in oxDNA
    # they point in the 3' to 5' direction. Ben 29/11/13
    rev_sys = base.System(final_sys._box)
    for strand in final_sys._strands:
        reverse_nucs = [nuc for nuc in strand._nucleotides]
        reverse_nucs.reverse()
        rev_strand = base.Strand()
        for nuc in reverse_nucs:
            rev_strand.add_nucleotide(base.Nucleotide(nuc.cm_pos, nuc._a1, -nuc._a3, nuc._base, nuc._btype))
        if strand._circular:
            rev_strand.make_circular(check_join_len=True)
        rev_sys.add_strand(rev_strand, check_overlap=False)
    # also reverse the vhelix_vbase_to_nucleotide order so it corresponds to the reversed system
    vh_vb2nuc_rev = cu.vhelix_vbase_to_nucleotide()
    # count the number of nucleotides up to but not including the nucleotides in strand ii
    nnucs_to_here = list(range(rev_sys._N_strands))
    nuc_total = 0
    for strandii, strand in enumerate(rev_sys._strands):
        nnucs_to_here[strandii] = nuc_total
        nuc_total += len(strand._nucleotides)

    id_to_pos = {}
    # fill in the _scaf and _stap dicts for the reverse vhelix_vbase_to_nucleotide object
    for vh, vb in list(vh_vb2nuc_final._scaf.keys()):
        strandii, nuciis = vh_vb2nuc_final._scaf[(vh, vb)]
        rev_nuciis = []
        for nucii in nuciis:
            nuc = len(rev_sys._strands[strandii]._nucleotides) - 1 - (nucii - nnucs_to_here[strandii]) + nnucs_to_here[strandii]
            rev_nuciis.append(nuc)
            id_to_pos[nuc] = (vh, vb, True)
        vh_vb2nuc_rev.add_scaf(vh, vb, strandii, rev_nuciis)
    for vh, vb in list(vh_vb2nuc_final._stap.keys()):
        strandii, nuciis = vh_vb2nuc_final._stap[(vh, vb)]
        rev_nuciis = []
        for nucii in nuciis:
            nuc = len(rev_sys._strands[strandii]._nucleotides) - 1 - (nucii - nnucs_to_here[strandii]) + nnucs_to_here[strandii]
            rev_nuciis.append(nuc)
            id_to_pos[nuc] = (vh, vb, False)

        vh_vb2nuc_rev.add_stap(vh, vb, strandii, rev_nuciis)


    if print_lattice_id_map:
        def search_key(dict, val):
            key_lst = []
            for key in dict:
                if dict[key] == val:
                    key_lst.append(key)
            return key_lst
        
        pos_to_id = {}
        for pos in id_to_pos.values():
            if pos not in pos_to_id.keys():
                pos_to_id[pos] = search_key(id_to_pos, pos)

        with open(output_file+"_cadnano2oxDNA_map.pickle", "wb") as map_dic:
            pickle.dump(pos_to_id, map_dic, protocol=pickle.HIGHEST_PROTOCOL)

        with open(output_file+"_bond_pairs.txt", "w") as pair_file:
            for k in pos_to_id.keys():
                if k[2] == True:
                    if (k[0], k[1], False) in pos_to_id.keys():
                        for i in range(len(pos_to_id[k])):
                            pair_file.write("{} {}\n".format(pos_to_id[k][i], pos_to_id[(k[0],k[1],False)][i]))
    

    # dump the spatial arrangement of the vhelices to a file
    vhelix_pattern = {}
    for i in range(len(cadsys.vhelices)):
        vhelix_pattern[cadsys.vhelices[i].num] = (cadsys.vhelices[i].row,cadsys.vhelices[i].col)
        
    if rev_sys.N == 0:
        base.Logger.log("The generated configuration is empty: this might be due to this conversion module not supporting virtual helices containing no scaffold strands.", base.Logger.CRITICAL)
        exit(1)

    if print_virt2nuc:
        with open("virt2nuc", "wb") as fout:
            pickle.dump((vh_vb2nuc_rev, vhelix_pattern), fout)
            print("## Wrote nucleotides' index conversion data to virt2nuc", file=sys.stderr)

    topology_file = output_file + ".top"
    configuration_file = output_file + ".oxdna"

    if print_oxview:
        # Find colors
        pos_to_color = {}
        for vh in cadsys.vhelices:
            for c in vh.stap_colors:
                vb, color = c
                pos_to_color[(vh.num, vb)] = color

        # Create inverse mapping to find pairs by their positions
        pos_to_id = {}
        for bid, helix in id_to_pos.items():
            pos_to_id.setdefault(helix, []).append(bid)

        # Create mapping to find nucleotide by its index
        id_to_nucleotide = {n.index: n for s in rev_sys._strands for n in s._nucleotides}

        # Find the index of the scaffold strand (there really should be some
        # easier way of doing this)
        strands = [id_to_nucleotide[nucId].strand for nucIds in pos_to_id.values() for nucId in nucIds]
        scaffold_index = max(set(strands), key=strands.count)

        # Make one cluster per domain
        cluster_ids = {}

        # Go through system and add extra information
        # (basepairs, clusters and colors)
        for s in rev_sys._strands:
            for n in s._nucleotides:
                if n.index in id_to_pos:
                    # Find which helix and position the nucleotide had
                    vh, vb = id_to_pos[n.index]
                    paired = pos_to_id[(vh,vb)]

                    # If double-stranded, save basepair
                    if len(paired) > 1:
                        n.pair = id_to_nucleotide[[idx for idx in paired if idx != n.index][0]]

                    # Get the staple strand nucleotide at this position
                    staple_ids = [n for n in pos_to_id[(vh,vb)] if id_to_nucleotide[n].strand != scaffold_index]
                    if len(staple_ids) > 0:
                        staple_nuc = id_to_nucleotide[staple_ids[0]]
                        # One cluster per combination of helix and staple strand
                        n.cluster = cluster_ids.setdefault((vh, staple_nuc.strand), len(cluster_ids)+1)

                        # If this is the staple strand, check if it has a color
                        if n is staple_nuc and (vh, vb) in pos_to_color:
                            n.color = pos_to_color[(vh, vb)]
                    else:
                        n.cluster = -1
        # If there's any colored nucleotide in a strand,
        # use that color for the whole strand
        for s in rev_sys._strands:
            for n in s._nucleotides:
                if n.color is not None:
                    for other in s._nucleotides:
                        other.color = n.color
                    break

        # Print the oxview output
        rev_sys.print_oxview_output(output_file+'.oxview')

    rev_sys.print_lorenzo_output(configuration_file, topology_file)
    
    print("## Wrote data to '%s' / '%s'" % (configuration_file, topology_file), file=sys.stderr)
    print("## DONE", file=sys.stderr)

def main():
    source_file, lattice_type, sequence_filename, side, np_seed, print_lattice_id_map, print_virt2nuc, print_oxview = readingCli()
    cadsys, sequences = parsingCli(source_file, sequence_filename)
    output_file = os.path.abspath(output_file)
    cadnano_oxdna(output_file, cadsys, lattice_type, sequences, side, np_seed, print_lattice_id_map, print_virt2nuc, print_oxview)


if __name__ == '__main__':
    main()