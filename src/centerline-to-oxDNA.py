import numpy as np
import math as mt
import numpy.linalg as la
import sys
import os
from libs import topology as top
from libs import base

# angolo definito dando una normale fissata del piano
def py_ang(v1, v2, vplane):
	v1n = v1 / la.norm(v1) 
	v2n = v2 / la.norm(v2) 

	
	# arctan con questi argomenti da angolo tra 0 e pi, e la normale definisce la direzione. Se vuoi fissare la direzione, il segno lo decide quale delle due normali scegli
	return np.arctan2 (la.norm(np.cross(v1n, v2n)) , np.dot (v1n, v2n)) * np.sign(np.dot(np.cross(v1n, v2n), vplane))
	
	
# matrice rotazione asse-angolo (asse non serve normalizzato, angolo in rad)
def get_rotation_matrix(axis, anglest):
	# the argument anglest can be either an angle in radiants
	# (accepted types are float, int or np.float64 or np.float64)
	# or a tuple [angle, units] where angle a number and
	# units is a string. It tells the routine whether to use degrees,
	# radiants (the default) or base pairs turns
	if not isinstance (anglest, (np.float64, np.float32, float, int)):
		if len(anglest) > 1:
			if anglest[1] in ["degrees", "deg", "o"]:
				# angle = np.deg2rad (anglest[0])
				angle = (np.pi / 180.) * (anglest[0])
			elif anglest[1] in ["bp"]:
				angle = int(anglest[0]) * (np.pi / 180.) * (34.81)
			else:
				angle = float(anglest[0])
		else:
			angle = float(anglest[0])
	else:
		angle = float(anglest)  # in degrees, I think
	# print "angolo",angle/np.pi*180
	axis = np.array(axis)
	axis /= np.sqrt(np.dot(axis, axis))

	ct = np.cos(angle)
	st = np.sin(angle)
	olc = 1. - ct
	x, y, z = axis

	# questa e' la matrice che torna
	return np.array([[olc * x * x + ct, olc * x * y - st * z, olc * x * z + st * y],
					[olc * x * y + st * z, olc * y * y + ct, olc * y * z - st * x],
					[olc * x * z - st * y, olc * y * z + st * x, olc * z * z + ct]])
	
	
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
			print >> sys.stderr, "The --nicked and --ssDNA options are incompatible"
			exit(1)
		
		
def print_usage():
        print >> sys.stderr, "USAGE:"
        print >> sys.stderr, "\t%s centerline_file" % sys.argv[0]
        print >> sys.stderr, "\t[-c|--closed] [-o|--open] [-h|--help] [-d\--dsDNA] [-s\--ssDNA] [-n\--nicked] [-p\--supercoiling] [-w\--writhe] [-e\--seed] [-q\--sequence]"
        exit(1)
		
		
def parse_options(argv):
	shortArgs = 'cohdsnp:w:e:q:'
	longArgs = ['closed', 'open', 'help', 'dsDNA', 'ssDNA', 'nicked', 'supercoiling', 'writhe', 'seed', 'sequence']
	
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
			if k[0] == '-p' or k[0] == "--supercoiling=": opts.supercoiling = float(k[1])
			if k[0] == '-w' or k[0] == "--writhe=": opts.writhe = float(k[1])
			if k[0] == '-e' or k[0] == "--seed=": opts.seed = int(k[1])
			if k[0] == '-q' or k[0] == "--sequence=": opts.sequence_file = k[1]
			
		opts.centerline_file = files[0]
	except Exception as e:
		print_usage()
		
	return opts

# parametri
BASE_BASE = 0.3897628551303122  # in direzione normale all'elica, va contata anche la rotazione
CM_CENTER_DS = 0.6

if __name__ == '__main__':
	opts = parse_options(sys.argv)
	opts.check()
	
	if opts.seed != None:
		np.random.seed(opts.seed)
	
	# importa file con coordinate
	coordxyz = np.loadtxt(opts.centerline_file, float)
	
	numrows = len(coordxyz)  # num base pairs
	
	# scaliamo i vettori prendendo 0.3897628551303122 come rif 
	scaling = BASE_BASE / np.sqrt(np.dot(coordxyz[1, :] - coordxyz[0, :], coordxyz[1, :] - coordxyz[0, :]))
	coordxyz *= scaling  # coord che lo strand deve seguire
	
	
	# inizializzazione vettori legati agli strand
	# geometria asse
	dist = np.copy(coordxyz)
	dist_norm = np.copy(coordxyz)
	p = np.copy(coordxyz)
	
	ssdna1 = np.copy(coordxyz)
	v_perp_ssdna1 = np.copy(coordxyz)
	
	ssdna2 = np.copy(coordxyz)
	v_perp_ssdna2 = np.copy(coordxyz)
	
	
	######
	# #Size box sistema## si prende lato massimo
	####
	
	boxx = max(coordxyz[:numrows, 0]) - min(coordxyz[:numrows, 0])
	boxy = max(coordxyz[:numrows, 1]) - min(coordxyz[:numrows, 1])
	boxz = max(coordxyz[:numrows, 2]) - min(coordxyz[:numrows, 2])
	
	boxmax = max(boxx, boxy, boxz)
	
	# scegliamo come boxmax il diametro del cerchio di N particelle a distanza BASE_BASE tra loro
	boxmax = BASE_BASE * numrows / np.pi
	
	# calcoliamo vettori distanza tra elementi successivi asse (non normalizzati)
	for c in range(0, numrows): 
		ind = int(c - mt.floor(c / float(numrows)) * numrows)
		ind1 = int(c + 1 - mt.floor((c + 1) / float(numrows)) * numrows)
	
		dist[ind, :] = coordxyz[ind1, :] - coordxyz[ind, :]
		dist_norm[ind, :] = dist[ind, :] / np.sqrt(np.dot(dist[ind, :], dist[ind, :]))
	
	# calcoliamo vettori perpendicolare a due segmenti successivi (normalizzati)
	for c in range(0, numrows): 
		ind_1 = int(c - 1 - mt.floor((c - 1) / float(numrows)) * numrows)
		ind = int(c - mt.floor(c / float(numrows)) * numrows)
		
			
		p[ind, :] = np.cross(dist[ind_1, :] , dist[ind, :])
		p[ind, :] /= np.sqrt(np.dot(p[ind, :] , p[ind, :]))
	
	# calcoliamo il writhe della configurazione, per capire se stiamo dando una rotazione eccessiva
	WR = 0.
	if opts.closed:
		WR = top.get_writhe(coordxyz)
	
	# il pitch di equilibrio dipende da larghezza sistema, tra 10.50 e 10.55 (N/pitch da il numero di giri di 360 totali, ossia il linking number)
	# per sistemi di 10080 basi ho visto che e' perfettamente 10.50 a 1M
	pitch = 10.5
	LK = round((numrows / pitch) * (opts.supercoiling + 1) + opts.writhe)  # LK deve essere un numero intero? Si perche devi chiudere il giro alla fine
	
	TW = LK - WR  # TW = LK - WR
	rot_base = TW * 2.0 * np.pi / numrows  # angolo in radianti per base 
	
	####################################
	# ##Inizializziamo lo strand di dna
	####################################
	
	# posizione iniziale v_perp_ssdna1 (direzione tra asse e ssdna1), sfruttando il fatto che p(0) e perp sia a dist(-1) che dist(0). Dopo anche se ruoto non cambia nulla. Essendo normalizzati anche v_perp_ssdna1 lo e, ma lo normalizziamo per sicurezza
	v_perp_ssdna1[0, :] = np.cross(dist_norm[0, :], p[0, :]) 
	v_perp_ssdna1[0, :] /= np.sqrt(np.dot(v_perp_ssdna1[0, :], v_perp_ssdna1[0, :])) 
	
	for c in range(numrows): 
		# pos ssdna1
		ssdna1[c, :] = coordxyz[c, :] - CM_CENTER_DS * v_perp_ssdna1[c, :] 
		ssdna2[c, :] = coordxyz[c, :] + CM_CENTER_DS * v_perp_ssdna1[c, :] 
	
		# Aggiorniamo v_perp_ssdna1 in modo che l angolo di twist sia rot_base (vedi Langowski per il calcolo dell angolo)
		
		# riscrivo gli indici per non confondermi
		ind_1 = int(c - mt.floor(c / float(numrows)) * numrows)  # indice del basepair attuale
		ind = int(c + 1 - mt.floor((c + 1) / float(numrows)) * numrows)  # indice del nuovo basepair
		
		alpha = py_ang(v_perp_ssdna1[ind_1, :] , p[ind, :] , dist[ind_1, :])
		gamma = rot_base - alpha
		# gamma = gamma -  int(round(  (gamma)  /(2* np.pi)))  * 2*np.pi # non so se e necessario
	
		# usiamo gamma per ruotare p_i sull asse s_(i+1) in modo da avere l angolo corretto con v_perp_ssdna1 vecchio
		
		R = get_rotation_matrix(dist[ind, :], gamma)
		v_perp_ssdna1[ind, :] = np.dot(R , p[ind, :])  # normalizzato
		v_perp_ssdna2[ind, :] = -v_perp_ssdna1[ind, :]
		
	# check LK imposed and measured
	TW_measured = top.get_twist(coordxyz, ssdna1)
	
	box = np.array([2 * boxmax, 2 * boxmax, 2 * boxmax])
	system = base.System(box)
	
	if opts.sequence_file == None:
		ssdna1_base = np.zeros(numrows, int)
		# assegniamo basi in modo casuale
		for c in range(numrows):
			ssdna1_base[c] = np.random.randint(0, 4)
	else:
		try:
			seq_file = open(opts.sequence_file)
		except Exception:
			print >> sys.stderr, "The sequence file '%s' is unreadable" % opts.sequence_file
			exit(1)
			
		contents = seq_file.read()
		# remove all whitespace from the file's contents
		sequence = ''.join(contents.split())
		if len(sequence) != numrows:
			print >> sys.stderr, "The length of the given sequence (%d) should be equal to the number of coordinates in the centerline file (%d)" % (len(sequence), numrows)
			exit(1)
			
		ssdna1_base = map(lambda x: base.base_to_number[x], sequence)		
			
		seq_file.close()
			

	strand1 = base.Strand()	
	for c in range(numrows):
		b = ssdna1_base[c]
		strand1.add_nucleotide(base.Nucleotide(ssdna1[c], v_perp_ssdna1[c], dist_norm[c], b, b))
	if opts.closed:
		strand1.make_circular()
	system.add_strand(strand1)

	if opts.double:
		strand2 = base.Strand()
		for c in range(numrows):
			reverse_idx = numrows - 1 - c
			b = 3 - ssdna1_base[reverse_idx]
			strand2.add_nucleotide(base.Nucleotide(ssdna2[reverse_idx], v_perp_ssdna2[reverse_idx], -dist_norm[reverse_idx], b, b))
		if opts.closed and not opts.nicked:
			strand2.make_circular()
		system.add_strand(strand2)
		
	basename = os.path.basename(sys.argv[1])
	system.print_lorenzo_output(basename + ".oxdna", basename + ".top")
