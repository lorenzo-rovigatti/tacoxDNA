#!/usr/bin/env python3

import sys
import re
import os
from .libs.pyquaternion import Quaternion
import numpy as np
from .libs import cadnano_utils as cu
from .libs import base


class Options(object):

	def __init__(self):
		object.__init__(self)
		
		self.seed = None
		self.file_name_in = None
		self.print_oxview = None
		
	def check(self):
		return True


def move_along_vector(point, vector, length):  # rpoly file contains center coordinate of helix, "generate" needs end coordiates of helix:
	move_distance = float(length) * 0.4 / 2.0  # 0.4 is the length of a base pair in oxDNA units, move half the helixlength down
	return [point[0] - move_distance * vector[0] , point[1] - move_distance * vector[1] , point[2] - move_distance * vector[2] ]


def rpoly_to_oxDNA(opts):
	# Read File
	# 'data' stores helix coordinates + rotaion in quaternion
	data = []
	
	rev_helix_connections = []  # staple connection information,
	fwd_helix_connections = []  # scaffold connections
	count = 0
	polyFile = open(opts.file_name_in, 'r')
	
	try:
		for line in polyFile:
			if line.startswith('hb'):
				data.insert(count, line.split(' '))
				count += 1
			elif line.startswith('c'):
				if 'f3' not in line:
					rev_helix_connections.append([int(re.search('c helix_(.+?) ', line).group(1)), int(re.search('\' helix_(.+?) ', line).group(1))])  # Extract connection information
				else:
					fwd_helix_connections.append([int(re.search('c helix_(.+?) ', line).group(1)), int(re.search('\' helix_(.+?) ', line).group(1))])
	except Exception:
		print('Failed to read the file')
	
	generator = cu.StrandGenerator()
	
	staple_fragments = base.System([100, 100, 100])  # temporary system to store staple fragments before later connecting them
	scaffold_fragments = base.System([100,100,100])
	
	# Reads orientation from the "data" and produces rotations from the Quaternian coordinates
	largest_size = 0.0
	for n, i in enumerate(data):
	
		position = [float(i[3]) / 0.84 , float(i[4]) / 0.84 , float(i[5]) / 0.84]  # 0.84 scaling is ad hoc solution to get good looking models
		
		n_bp = int(i[2])
		
		q = Quaternion(w=float(i[9]), x=float(i[6]), y=float(i[7]), z=float(i[8]))  # find the helix roation Info from file
		vec = q.rotate(np.array([0.0, 0.0, 1.0]))  # use it to figure out direction
		vec2 = q.rotate([0.65, -0.76, 0.0])  # ad hoc onversion between rpoly rotation and cadnano utils
		
		new_position = move_along_vector(position , vec , n_bp)  # rpoly helix coordinates are defined in center of helix, cadnano utils have positions in the base of helix.
		
		for j in new_position:  # go through every coordinate to find the largest coordinate to figure out box size
			if j > largest_size:
				largest_size = j
			else:
				pass
		
		# strand 0 is the scaffold and strand 1 is the staple
		new_strands = generator.generate_or_sq(bp=n_bp, start_pos=new_position, direction=vec, perp=vec2)

		#for oxview export, cluster nucleotide by helix
		if opts.print_oxview is not None:
			for i in [0,1]:
				for nucleotide in new_strands[i]._nucleotides:
					nucleotide.cluster = n+1

		# cut strand 1 into two equal lengh staple fragments for later connections
		fragment1, fragment2 = new_strands[1].cut_in_two(copy=False)
		
		# store the fragments in this system for later connections
		staple_fragments.add_strand(fragment1)
		staple_fragments.add_strand(fragment2)

		scaffold_fragments.add_strand(new_strands[0])

	
	output_system = base.System([largest_size * 3.0, largest_size * 3.0, largest_size * 3.0])
	for n in rev_helix_connections:  # iterate through staple strand connections and connect the previously generated fragments
		connect_from = n[0] * 2 - 1
		connect_to = n[1] * 2 - 2
		staple_strand = staple_fragments._strands[connect_from]
		staple_strand = staple_strand.append(staple_fragments._strands[connect_to])
	
		output_system.add_strand(staple_strand)

	scaffold_strand = scaffold_fragments._strands[0]
	for n in fwd_helix_connections[:-1]:
		next_segment_adress = n[1]-1
		next_segment = scaffold_fragments._strands[next_segment_adress]
		scaffold_strand = scaffold_strand.append(next_segment)



	scaffold_strand.make_circular()
	output_system.add_strand(scaffold_strand)
	
	basename = os.path.basename(opts.file_name_in)
	top_file = basename + ".top"
	conf_file = basename + ".oxdna"
	
	output_system.print_lorenzo_output(conf_file, top_file)

	if opts.print_oxview is not None:
		oxview_file = basename + ".oxview"
		output_system.print_oxview_output(oxview_file)


def print_usage():
	print("USAGE:", file=sys.stderr)
	print("\t%s rpoly_file" % sys.argv[0], file=sys.stderr)
	print("\t[-e\--seed=VALUE]", file=sys.stderr)
	print("\t[-o\--print-oxview]", file=sys.stderr)
	exit(1)

	
def parse_options():
	shortArgs = 'e:o'
	longArgs = ['seed=', '--print-oxview']
	
	opts = Options()
	
	try:
		import getopt
		args, files = getopt.gnu_getopt(sys.argv[1:], shortArgs, longArgs)
		for k in args:
			if k[0] == '-e' or k[0] == "--seed": 
				opts.seed = int(k[1])
			if k[0] == '-o' or k[0] == "--print-oxview":
				opts.print_oxview = True
			
		opts.file_name_in = files[0]
	except Exception:
		print_usage()
		
	return opts


def main():
	if len(sys.argv) < 2:
		print_usage()
	
	opts = parse_options()
	if opts.seed is not None:
		np.random.seed(opts.seed)
	rpoly_to_oxDNA(opts)


if __name__ == '__main__':
	main()
