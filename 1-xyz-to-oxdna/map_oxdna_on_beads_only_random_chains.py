#primo argomento: file di input # secondo argomento: 0 non-nicked 1 nicked 2 nicked senza una base # terzo argomento pitch

import numpy as np
import math as mt
import numpy.linalg as la
import sys
import function_writhe as fw
import function_twist as ft

#angolo definito dando una normale fissata del piano
def py_ang(v1, v2,vplane):
	v1n = v1 / la.norm( v1 ) 
	v2n = v2 / la.norm( v2 ) 

	
	#arctan con questi argomenti da angolo tra 0 e pi, e la normale definisce la direzione. Se vuoi fissare la direzione, il segno lo decide quale delle due normali scegli
	return np.arctan2 (  la.norm(np.cross(v1n, v2n)) , np.dot (v1n, v2n) ) * np.sign( np.dot( np.cross(v1n,v2n), vplane) )
	
	
#matrice rotazione asse-angolo (asse non serve normalizzato, angolo in rad)
def get_rotation_matrix(axis, anglest):
	# the argument anglest can be either an angle in radiants
	# (accepted types are float, int or np.float64 or np.float64)
	# or a tuple [angle, units] where angle a number and
	# units is a string. It tells the routine whether to use degrees,
	# radiants (the default) or base pairs turns
	if not isinstance (anglest, (np.float64, np.float32, float, int)):
		if len(anglest) > 1:
			if anglest[1] in ["degrees", "deg", "o"]:
				#angle = np.deg2rad (anglest[0])
				angle = (np.pi / 180.) * (anglest[0])
			elif anglest[1] in ["bp"]:
				angle = int(anglest[0]) * (np.pi / 180.) * (34.81)
			else:
				angle = float(anglest[0])
		else:
			angle = float(anglest[0])
	else:
		angle = float(anglest) # in degrees, I think
	#print "angolo",angle/np.pi*180
	axis = np.array(axis)
	axis /= np.sqrt(np.dot(axis, axis))

	ct = np.cos(angle)
	st = np.sin(angle)
	olc = 1. - ct
	x, y, z = axis

	#questa e' la matrice che torna
	return np.array([[olc*x*x+ct, olc*x*y-st*z, olc*x*z+st*y],
					[olc*x*y+st*z, olc*y*y+ct, olc*y*z-st*x],
					[olc*x*z-st*y, olc*y*z+st*x, olc*z*z+ct]])

#parametri
BASE_BASE = 0.3897628551303122  # in direzione normale all'elica, va contata anche la rotazione
CM_CENTER_DS=0.6


#importa file con coordinate
coordxyz=np.loadtxt(sys.argv[1],float)

nicking=int(sys.argv[2])

print >> sys.stderr, nicking

numrows = len(coordxyz)  #num base pairs

#scaliamo i vettori prendendo 0.3897628551303122 come rif 
scaling=BASE_BASE/np.sqrt(np.dot(coordxyz[1,:]-coordxyz[0,:],coordxyz[1,:]-coordxyz[0,:]))
coordxyz*=scaling	#coord che lo strand deve seguire


#inizializzazione vettori legati agli strand
#geometria asse
dist=np.copy(coordxyz)
dist_norm=np.copy(coordxyz)
p=np.copy(coordxyz)

ssdna1=np.copy(coordxyz)
v_perp_ssdna1=np.copy(coordxyz)

ssdna2=np.copy(coordxyz)
v_perp_ssdna2=np.copy(coordxyz)


######
##Size box sistema## si prende lato massimo
####

boxx=max(coordxyz[:numrows,0])-min(coordxyz[:numrows,0])
boxy=max(coordxyz[:numrows,1])-min(coordxyz[:numrows,1])
boxz=max(coordxyz[:numrows,2])-min(coordxyz[:numrows,2])

boxmax=max(boxx,boxy,boxz)

#scegliamo come boxmax il diametro del cerchio di N particelle a distanza BASE_BASE tra loro
boxmax=BASE_BASE*numrows/np.pi

#calcoliamo vettori distanza tra elementi successivi asse (non normalizzati)
for c in range(0,numrows): 
	ind = int( c - mt.floor(c/ float(numrows) ) * numrows		)
	ind1 =int( c+1 - mt.floor((c+1)/ float(numrows) ) * numrows	)

	dist[ind,:]=coordxyz[ind1,:]-coordxyz[ind,:]
	dist_norm[ind,:]=dist[ind,:]/np.sqrt(np.dot(dist[ind,:],dist[ind,:]))
	#print ind,ind1, dist[ind,:],coordxyz[ind,:],coordxyz[ind1,:]

#calcoliamo vettori perpendicolare a due segmenti successivi (normalizzati)
for c in range(0,numrows): 
	ind_1 = int( c-1 - mt.floor((c-1)/ float(numrows) ) * numrows	)
	ind = int( c - mt.floor(c/ float(numrows) ) * numrows		)

	
	
		
	p[ind,:] = np.cross( dist[ind_1,:] , dist[ind,:] )
	p[ind,:]/= np.sqrt(np.dot( p[ind,:] , p[ind,:] ))
	#print c,ind_1,ind,coordxyz[ind_1,:]/scaling,coordxyz[ind,:]/scaling,p[ind,:]
	#print p[ind,0], p[ind,1], p[ind,2]

#calcoliamo il writhe della configurazione, per capire se stiamo dando una rotazione eccessiva
WR=fw.get_writhe(coordxyz)
print >> sys.stderr, "Writhe",WR


#il pitch di equilibrio dipende da larghezza sistema, tra 10.50 e 10.55 (N/pitch da il numero di giri di 360 totali, ossia il linking number)
#per sistemi di 10080 basi ho visto che e' perfettamente 10.50 a 1M
pitch=float(sys.argv[3])
write_nodo=float(sys.argv[4])
sigma_supercoil=float(sys.argv[5])
LK = round((numrows / pitch)*(sigma_supercoil+1) + write_nodo) #LK deve essere un numero intero? Si perche devi chiudere il giro alla fine


TW=LK-WR #TW = LK - WR
rot_base = TW * 2.0 * np.pi / numrows # angolo in radianti per base 

print >> sys.stderr, pitch,LK, TW, rot_base
####################################
###Inizializziamo lo strand di dna
####################################

#posizione iniziale v_perp_ssdna1 (direzione tra asse e ssdna1), sfruttando il fatto che p(0) e perp sia a dist(-1) che dist(0). Dopo anche se ruoto non cambia nulla. Essendo normalizzati anche v_perp_ssdna1 lo e, ma lo normalizziamo per sicurezza
v_perp_ssdna1[0,:] = np.cross(dist_norm[0,:],p[0,:]) 
v_perp_ssdna1[0,:] /= np.sqrt(np.dot(v_perp_ssdna1[0,:],v_perp_ssdna1[0,:])) 

for c in range(numrows): 
	
	
	#pos ssdna1
	ssdna1[c,:] = coordxyz[c,:] - CM_CENTER_DS * v_perp_ssdna1[c,:] 
	ssdna2[c,:] = coordxyz[c,:] + CM_CENTER_DS * v_perp_ssdna1[c,:] 

	#Aggiorniamo v_perp_ssdna1 in modo che l angolo di twist sia rot_base (vedi Langowski per il calcolo dell angolo)
	
	#riscrivo gli indici per non confondermi
	ind_1 = int( c - mt.floor(c/ float(numrows) ) * numrows	)			#indice del basepair attuale
	ind = int( c+1 - mt.floor((c+1)/ float(numrows) ) * numrows		)	#indice del nuovo basepair
	
	alpha = py_ang( v_perp_ssdna1[ind_1,:] , p[ind,:] , dist[ind_1,:])
	gamma = rot_base - alpha
	#gamma = gamma -  int(round(  (gamma)  /(2* np.pi)))  * 2*np.pi # non so se e necessario

	#print alpha,gamma,alpha+gamma #-  int(round(  (alpha+gamma)  /(2* np.pi)))  * 2*np.pi,rot_base #check superimportante , fallo tra -pi e pi


	#usiamo gamma per ruotare p_i sull asse s_(i+1) in modo da avere l angolo corretto con v_perp_ssdna1 vecchio
	
	R = get_rotation_matrix( dist[ind,:], gamma  )
	v_perp_ssdna1[ind,:] = np.dot( R , p[ind,:] ) #normalizzato
	v_perp_ssdna2[ind,:] = -v_perp_ssdna1[ind,:]
	
#check LK imposed and measured
TW_measured = ft.get_twist(coordxyz,ssdna1)
print >>  sys.stderr , TW, TW_measured,TW_measured+WR,LK
	

#stampiamo gli strand #NOTA ssdna2 va stampato inverso, altrimenti non ha elicita giusta nel programma
outfile = open ('generated.conf', 'w')

print >> outfile, "t = 0"
print >> outfile, "b = ", boxmax, boxmax, boxmax
print >> outfile, "E = 0. 0. 0."

for c in range(numrows): 
	print >> outfile,  ssdna1[c,0],ssdna1[c,1],ssdna1[c,2] , v_perp_ssdna1[c,0], v_perp_ssdna1[c,1], v_perp_ssdna1[c,2], dist_norm[c,0] ,dist_norm[c,1],dist_norm[c,2],0,0,0,0,0,0

if(nicking==0 or nicking==1):
	for c in reversed(range(numrows)):
		print >> outfile, ssdna2[c,0],ssdna2[c,1],ssdna2[c,2] , v_perp_ssdna2[c,0], v_perp_ssdna2[c,1], v_perp_ssdna2[c,2], -dist_norm[c,0] ,-dist_norm[c,1],-dist_norm[c,2],0,0,0,0,0,0
if(nicking==2): #togliamo l ultima base
	for c in reversed(range(1,numrows)):
		print >> outfile, ssdna2[c,0],ssdna2[c,1],ssdna2[c,2] , v_perp_ssdna2[c,0], v_perp_ssdna2[c,1], v_perp_ssdna2[c,2], -dist_norm[c,0] ,-dist_norm[c,1],-dist_norm[c,2],0,0,0,0,0,0



################
##stampiamo il file .top
############
number_to_base = {0 : 'A', 1 : 'G', 2 : 'C', 3 : 'T'}
base_to_number = {'A' : 0, 'a' : 0, 'G' : 1, 'g' : 1, 'C' : 2, 'c' : 2, 'T' : 3, 't' : 3}

ssdna1_base = np.zeros(numrows,int)
ssdna2_base = np.zeros(numrows,int)
#assegniamo basi in modo casuale
for c in range(numrows):
	ssdna1_base[c] = np.random.randint(0, 4)
	ssdna2_base[numrows-1-c] = (3-ssdna1_base[c])



out = open ("generated.top", "w")

#header .top
if(nicking==0 or nicking==1):
	print >> out, numrows*2, 2
if(nicking==2):
	print >> out, numrows*2-1, 2

#circolare chiuso
#ssdna1
for c in range(numrows):
	ind_1 = int ( c-1 - mt.floor((c-1)/ float(numrows) ) * numrows )
	ind   = int ( c   - mt.floor((c  )/ float(numrows) ) * numrows )
	ind1  = int ( c+1 - mt.floor((c+1)/ float(numrows) ) * numrows ) 
	print >> out, 1, number_to_base[ssdna1_base[ind]], ind_1 , ind1
			
#ssdna2
for c in range(numrows,2*numrows):
	ind_1 = int ( c-1 - mt.floor((c-1 - numrows )/ float(numrows) ) * numrows )
	ind   = int ( c   - mt.floor((c   - numrows ) / float(numrows) ) * numrows )
	ind1  = int ( c+1 - mt.floor((c+1 - numrows )/ float(numrows) ) * numrows ) 

	if(nicking==0):
		print >> out, 2, number_to_base[ssdna2_base[ind-numrows]], ind_1 , ind1

	if(nicking==1):
		if(c==numrows):
			print >> out, 2, number_to_base[ssdna2_base[ind-numrows]], -1 , ind1
		elif(c==2*numrows-1):
			print >> out, 2, number_to_base[ssdna2_base[ind-numrows]], ind_1 , -1
		else:
			print >> out, 2, number_to_base[ssdna2_base[ind-numrows]], ind_1 , ind1

	if(nicking==2):
		if(c==numrows):
			print >> out, 2, number_to_base[ssdna2_base[ind-numrows]], -1 , ind1
		elif(c==2*numrows-2):
			print >> out, 2, number_to_base[ssdna2_base[ind-numrows]], ind_1 , -1
		elif(c!=2*numrows-1):
			print >> out, 2, number_to_base[ssdna2_base[ind-numrows]], ind_1 , ind1
	
			











