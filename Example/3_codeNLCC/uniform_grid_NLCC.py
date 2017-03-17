#/usr/bin/python

import sys, time, os, numpy, scipy
from numpy import *
from numpy.fft import *
from scipy import *
from scipy.interpolate import *
from scipy.integrate import *
import sys, time, os
from numpy import *

maxr=14.0
dr=0.01
nr=int(maxr/dr)

def getfromfile4(filename):
    x = []
    lines = os.popen('cat '+filename+' | awk \'{print $1" "$2" "$3" "$4}\'').read().strip().split('\n')
    for line in lines:
        l = line.split(' ')
        x.append((float(l[0]), float(l[1]), float(l[2]), float(l[3])))
    return x

print 'read in data from NLCC.txt'
raw_data_file='NLCC.txt'
vraw=getfromfile4(raw_data_file)

vr_raw, nlcc0_raw, nlcc1_raw, nlcc2_raw = [0.0], [0.0], [0.0], [0.0]

for i in range(0, len(vraw)):
		vr_raw.append(vraw[i][0])
		nlcc0_raw.append(vraw[i][1])
		nlcc1_raw.append(vraw[i][2])
		nlcc2_raw.append(vraw[i][3])

#ftmp = open('tmp_NLCC.txt', 'w')
#for i in range(1, len(vraw)):
#	ftmp.write('%15.10f %15.10f %15.10f %15.10f \n' %(vraw[i][0],vraw[i][1], vraw[i][2], vraw[i][3]))
#ftmp.close()

print 'setup an uniform real-space grid, the length is ', maxr, ' and dr is ', dr 
r = arange(0.00, maxr, dr)
nlcc0, nlcc1, nlcc2 = [], [], [] 

print 'do interpolation'
f1 = interp1d(vr_raw, nlcc0_raw)
f2 = interp1d(vr_raw, nlcc1_raw)
f3 = interp1d(vr_raw, nlcc2_raw)

print 'print out the NLCC on an uniform grid in ABACUS_NLCC.txt' 
f = open('ABACUS_NLCC.txt', 'w')
f.write('%d %15.10f %15.10f %15.10f %15.10f 0.0000000000000000 0.0000000000000000\n' %(1, 0.0, vraw[0][1], vraw[0][2], vraw[0][3] ))
for i in range(1, len(r)):
	f.write('%d %15.10f %15.10f %15.10f %15.10f 0.0000000000000000 0.0000000000000000\n' %(i+1, r[i],f1(r[i]), f2(r[i]), f3(r[i])))
f.close()


print 'print out the NLCC on an uniform grid in PROFESS_NLCC.txt' 
f = open('PROFESS_NLCC.txt', 'w')
f.write('%d %15.10f\n' % (nr,dr))
f.write('%15.10f %15.10f %15.10f\n' %(r[0], vraw[0][1], vraw[0][2] ))
for i in range(1, len(r)):
	f.write('%15.10f %15.10f %15.10f\n' %(r[i], f1(r[i]), f2(r[i])))
f.close()

print ''
print 'All done.' 
print 'Have a great day!'
