import numpy as np
#import function_writhe as fw
#import function_twist as ft
#import function_linking_number as flk
import sys

#import dsdna as two ssdna, one after the other in order
coordxyz=np.loadtxt(sys.argv[1],float)
numrows = len(coordxyz)
flag=int(sys.argv[2])

#Distinguiamo il caso nick 2
if(int(numrows)%2!=0):
        nicking=2
        nbasi=(numrows+1)/2
else:
        nicking=0 #vale sia per nick 0 che 1
        nbasi=numrows/2

#split the array in to two strands
if(nicking==0):
        ssdna1 = coordxyz[0:nbasi,:]
        ssdna2 = np.flipud (coordxyz[nbasi:nbasi*2,:])
if(nicking==2):
        ssdna1 = coordxyz[1:nbasi,:]    #il nick toglie l ultima base del ssdna2 che si accoppia con la prima di ssdna1
        ssdna2 = np.flipud (coordxyz[nbasi:nbasi*2-1,:])
         
          
         
#print  ssdna1
#ricordiamo che ssdna2 e scritto al contrario
dna_axis= ( ssdna1 + ssdna2)/2
#print dna_axis
    
outfile = open ('centerline.dat', 'w')
for c in range(len(dna_axis)):         
       print >> outfile,  dna_axis[c][0],dna_axis[c][1],dna_axis[c][2] 
