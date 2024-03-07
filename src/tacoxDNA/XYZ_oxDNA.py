#!/usr/bin/env python3
# TODO: how do we choose the simulation box?

import numpy as np
import numpy.linalg as la
import sys
import os
from .libs import topology as top
from .libs import base
from .libs import utils

	
class Options(object):

	def __init__(self):
		object.__init__(self)
		
		self.closed = True
		self.double = True
		self.nicked = False
		self.supercoiling = 0.
		self.writhe = 0.
		self.seed = None
		self.sequence_file = None
		
	def check(self):
		if self.nicked and not self.double:
			print("The --nicked and --ssDNA options are incompatible", file=sys.stderr)
			exit(1)
		
		
def print_usage():
	print("USAGE:", file=sys.stderr)
	print("\t%s centerline_file" % sys.argv[0], file=sys.stderr)
	print("\t[-c\--closed] [-o\--open] [-h\--help] [-d\--dsDNA] [-s\--ssDNA] [-n\--nicked] [-p\--supercoiling VALUE] [-w\--writhe VALUE] [-e\--seed VALUE] [-q\--sequence FILE]", file=sys.stderr)
	exit(1)
		
		
def parse_options():
	shortArgs = 'cohdsnp:w:e:q:'
	longArgs = ['closed', 'open', 'help', 'dsDNA', 'ssDNA', 'nicked', 'supercoiling=', 'writhe=', 'seed=', 'sequence=']
	
	opts = Options()
	
	try:
		import getopt
		args, files = getopt.gnu_getopt(sys.argv[1:], shortArgs, longArgs)
		for k in args:
			if k[0] == '-c' or k[0] == '--closed': opts.closed = True
			if k[0] == '-o' or k[0] == '--open': opts.closed = False
			if k[0] == '-h' or k[0] == '--help': print_usage()
			if k[0] == '-d' or k[0] == '--dsDNA': opts.double = True
			if k[0] == '-s' or k[0] == "--ssDNA": opts.double = False
			if k[0] == '-n' or k[0] == "--nicked": opts.nicked = True
			if k[0] == '-p' or k[0] == "--supercoiling": opts.supercoiling = float(k[1])
			if k[0] == '-w' or k[0] == "--writhe": opts.writhe = float(k[1])
			if k[0] == '-e' or k[0] == "--seed": opts.seed = int(k[1])
			if k[0] == '-q' or k[0] == "--sequence": opts.sequence_file = k[1]
			
		opts.centerline_file = files[0]
	except Exception:
		print_usage()
		
	return opts


# base-base distance along the helical pitch
BASE_BASE = 0.3897628551303122
# distance between the helix centre and the nucleotides' centre of mass
CM_CENTER_DS = 0.6


def main():
	opts = parse_options()
	opts.check()
	
	if opts.seed != None:
		np.random.seed(opts.seed)
	
	# import the coordinates from the user-provided file
	coordxyz = np.loadtxt(opts.centerline_file, float)
	
	# number of base pairs
	nbases = len(coordxyz)
	
	# use the model parameters to scale the distances 
	scaling = BASE_BASE / la.norm(coordxyz[1, :] - coordxyz[0, :]) 
	coordxyz *= scaling

	#initialize vectors
	dist = np.copy(coordxyz)
	dist_norm = np.copy(coordxyz)
	p = np.copy(coordxyz)
	
	ssdna1 = np.copy(coordxyz)
	v_perp_ssdna1 = np.copy(coordxyz)
	
	ssdna2 = np.copy(coordxyz)
	v_perp_ssdna2 = np.copy(coordxyz)
	
	
	# take the bounding box as the simulation box for the output configuration
	boxx = max(coordxyz[:nbases, 0]) - min(coordxyz[:nbases, 0])
	boxy = max(coordxyz[:nbases, 1]) - min(coordxyz[:nbases, 1])
	boxz = max(coordxyz[:nbases, 2]) - min(coordxyz[:nbases, 2])
	boxmax = 1.5 * max(boxx, boxy, boxz)
	
	# centerline base_to_base vectors
	for c in range(0, nbases): 
		ind = c 
		ind1 = (c + 1) % nbases
		#open chain cannot compute dist at c=nbases-1 
		#so we use as reference the previus dist
		if opts.closed or (not opts.closed and c!=nbases-1):	
			dist[ind, :] = coordxyz[ind1, :] - coordxyz[ind, :]
			dist_norm[ind, :] = dist[ind, :] / la.norm(dist[ind, :])
		else:
			dist[ind, :] = dist[ind-1, :]
			dist_norm[ind, :] = dist_norm[ind-1, :]

	# vectors perpendicular between two consecutive centerline vectors (normalized)
	for c in range(0, nbases): 
		ind_1 = (c - 1 + nbases) % nbases 
		ind = c
			
		#opens chain have random p at c=0 and at c=nbases-1 due to absence of neighbours
		if opts.closed or (not opts.closed and c!=0 and c!=nbases-1):
			#check that dist[ind_1, :] and dist[ind, :] are not equals
			if not np.all(np.isclose( dist[ind_1, :] , dist[ind, :]) ):
				p[ind, :] = np.cross(dist[ind_1, :] , dist[ind, :])
				p[ind, :] /= la.norm(p[ind, :])
			#else assign random p
			else:
				rv=np.random.uniform(-1,1,3)
				p[ind, :] = rv - (np.dot(dist_norm[ind, :],rv)) * dist_norm[ind, :]
				p[ind, :] /= la.norm(p[ind, :])
		else:
			rv=np.random.uniform(-1,1,3)
			p[ind, :] = rv - (np.dot(dist_norm[ind, :],rv)) * dist_norm[ind, :]
			p[ind, :] /= la.norm(p[ind, :])

	# chain writhe
	WR = 0.
	if opts.closed:
		WR = top.get_writhe(coordxyz)

	#oxdna equiibrium pitch
	pitch = 10.5

	#global linking number
	LK = round((nbases / pitch) * (opts.supercoiling + 1) + opts.writhe) 
	
	TW = LK - WR 

	#twisting angle between two consecutive bases
	rot_base = TW * 2.0 * np.pi / nbases  

	#recap
	if opts.closed:
		print("Total Linking Number (LK) %f, composed of:" % ((nbases / pitch) * (opts.supercoiling + 1) + opts.writhe), file=sys.stderr)
		print("1) Equilibrium number of DNA turns %f" % (nbases / pitch), file=sys.stderr)
		print("2) Target writhe %f = Topological writhe %f + Turns imposed by supercoiling  %f " % (opts.writhe+(nbases / pitch) *opts.supercoiling,opts.writhe,(nbases / pitch) *opts.supercoiling), file=sys.stderr)
		print("LK has been rounded to %f" % LK, file=sys.stderr)
		print("Initial chain writhe %f" % WR, file=sys.stderr)


	####################################
	# Initialize DNA strand 
	####################################
	
	#First hydrogen-hydrogen bond vector
	v_perp_ssdna1[0, :] = np.cross(dist_norm[0, :], p[0, :]) 
	v_perp_ssdna1[0, :] /= la.norm(v_perp_ssdna1[0, :])
	v_perp_ssdna2[0, :] = -v_perp_ssdna1[0, :]
	
	for c in range(nbases): 
		#dna center of mass positions
		ssdna1[c, :] = coordxyz[c, :] - CM_CENTER_DS * v_perp_ssdna1[c, :] 
		ssdna2[c, :] = coordxyz[c, :] + CM_CENTER_DS * v_perp_ssdna1[c, :] 
	
		#Update v_perp_ssdna1
		
		ind_1 = c 
		ind = (c + 1) % nbases 
		
		alpha = top.py_ang(v_perp_ssdna1[ind_1, :] , p[ind, :] , dist[ind_1, :])
		gamma = rot_base - alpha
		#prevent change when ind=0 on open chain
		if c!=nbases-1:		
			R = utils.get_rotation_matrix(dist[ind, :], gamma)
			v_perp_ssdna1[ind, :] = np.dot(R , p[ind, :])  
			v_perp_ssdna2[ind, :] = -v_perp_ssdna1[ind, :]
		
	# check LK imposed and measured
	if opts.closed:
		TW_measured = top.get_twist(coordxyz, ssdna1)
	box = np.array([boxmax, boxmax, boxmax])
	system = base.System(box)
	
	if opts.sequence_file == None:
		ssdna1_base = np.zeros(nbases, int)
		for c in range(nbases):
			ssdna1_base[c] = np.random.randint(0, 4)
	else:
		try:
			seq_file = open(opts.sequence_file)
		except Exception:
			print("The sequence file '%s' is unreadable" % opts.sequence_file, file=sys.stderr)
			exit(1)
			
		contents = seq_file.read()
		# remove all whitespace from the file's contents
		sequence = ''.join(contents.split())
		if len(sequence) != nbases:
			print("The length of the given sequence (%d) should be equal to the number of coordinates in the centerline file (%d)" % (len(sequence), nbases), file=sys.stderr)
			exit(1)
			
		ssdna1_base = [base.base_to_number[x] for x in sequence]		
			
		seq_file.close()

	strand1 = base.Strand()	
	for c in range(nbases):
		b = ssdna1_base[c]
		strand1.add_nucleotide(base.Nucleotide(ssdna1[c], v_perp_ssdna1[c], dist_norm[c], b, b))
	if opts.closed:
		strand1.make_circular()
	system.add_strand(strand1)

	if opts.double:
		strand2 = base.Strand()
		for c in range(nbases):
			reverse_idx = nbases - 1 - c
			b = 3 - ssdna1_base[reverse_idx]
			strand2.add_nucleotide(base.Nucleotide(ssdna2[reverse_idx], v_perp_ssdna2[reverse_idx], -dist_norm[reverse_idx], b, b))
		if opts.closed and not opts.nicked:
			strand2.make_circular()
		system.add_strand(strand2)
		
	basename = os.path.basename(sys.argv[1])
	topology_file = basename + ".top"
	configuration_file = basename + ".oxdna"
	system.print_lorenzo_output(configuration_file, topology_file)
	
	print("## Wrote data to '%s' / '%s'" % (configuration_file, topology_file), file=sys.stderr)
	print("## DONE", file=sys.stderr)


if __name__ == '__main__':
	main()