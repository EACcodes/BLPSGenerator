#/usr/bin/python2.7
import sys, time, os, numpy, scipy
from numpy import *
from numpy.fft import *
from scipy import *
from scipy.interpolate import *
from scipy.integrate import *
from scipy.optimize import leastsq

def getxyfromfile1(filename):
    x = []
    lines = os.popen('cat '+filename+' | awk \'{print $1" "$2}\'').read().strip().split('\n')
    for line in lines:
        l = line.split(' ')
        x.append((float(l[0]), float(l[1])))
    return x

#----------------------------------------------------------
# (1) read in file name 
#----------------------------------------------------------
if len(sys.argv) < 2:
	print 'sort the Vq now'
	print 'usage: python this.py <file> nz'
	sys.exit(1)
else:
	file=sys.argv[1]
	z=int(sys.argv[2])
	print 'input file is ',file
	print 'z is ',z

#----------------------------------------------------------
# (2) define some basic parameters 
#----------------------------------------------------------
mev_ang3=4.03231375552813248587
tiny=1.0e-10
qrawmax=20.0

#----------------------------------------------------------
# (3) read in the data from file: data,
# the data contains v(q) from several structures.
# vraw[x][y]
# y=0 --> G norm
# y=1 --> Potential
#----------------------------------------------------------
vraw=getxyfromfile1(file)
vraw.sort()

#----------------------------------------------------------
# set the first value of the following arrays
# see the explanation next for these three.
#----------------------------------------------------------
qraw=[]
vqraw=[]
vnqraw=[]

#----------------------------------------------------------
# put all the value in qraw, vqraw, vnqraw,
# qraw: the G space vector length |q|
# vqraw: vqraw is the local potential
# vnqraw: vnqraw is the vqraw+4pi*Z/G^2,
#         which is the local pseudopotential
#         + z/r potential.
#----------------------------------------------------------
for i in range(0, len(vraw)):
	if vraw[i][0]>tiny and vraw[i][0]<qrawmax and abs(vraw[i][0]-vraw[i-1][0])>tiny:
		qraw.append(vraw[i][0])
		vqraw.append(vraw[i][1])
		vnqraw.append(vraw[i][1]+4*pi*z/vraw[i][0]**2)
print 'reorder, q, full_BLPS, non-Coulombic_BLPS.'


f = open('raw_data.txt', 'w')
f.write('qraw, full_BLPS, non-Coulombic_BLPS\n')
for i in range(0, len(qraw)):
    #f.write('%15.10f %15.10f %15.10f \n' %(qraw[i],vqraw[i]*mev_ang3,vnqraw[i]*mev_ang3))
    f.write('%15.10f %15.10f %15.10f \n' %(qraw[i],vqraw[i],vnqraw[i]))
