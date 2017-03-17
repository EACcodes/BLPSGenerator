#/usr/bin/python

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
# (1) Read in four parameters: 
# a) The atomic number 
# b) valence electron number 'Z', 
# c) initial V(q) value 'veff0'
# d) a real space radius cutoff 'vr_cut' (usually 4~6 Bohr)
#----------------------------------------------------------
if len(sys.argv) < 5:
	print '--- Now we are going to generate BLPS ---'
	print 'Usage: python this.py <z> <veff0> <vr_cut>'
	print 'The input effective potential should be in "BLPS.txt"'
	sys.exit(1)
else:
	atomic_number=int(sys.argv[1])
	z=int(sys.argv[2])
	veff0=float(sys.argv[3])
	vr_cutoff=float(sys.argv[4])
	print 'The atomic number is ', atomic_number
	print 'Input valence electron number Z is ',z
	print 'Radius cutoff in real space is ', vr_cutoff, ' Bohr'

print '(1) Read in parameters'

#----------------------------------------------------------
# (2) Define some basic parameters 
# mev_ang3 = 13.6058*2*0.529177*0.529177*0.529177
# tiny: smallest possible value of V_eff(r)
# qrawmax : the max value of q allowded.
# qmax: define the max value of uniform q mesh
# dq: define the increasement of q of uniform q mesh
#----------------------------------------------------------
max_r=14.01
dr=0.01
nr=int(max_r/dr)-1
bohr_to_a=0.529177249
mev_ang3=4.03231375552813248587
tiny=1.0e-10

qrawmax=20.0
qmax=30.0
dq=0.003
begin_r=1000

print '(2) Setup parameters'

#----------------------------------------------------------
# (3.1) Read in the potential file and sort the potential 
# according to |q|. The definition of vraw[x][y]:
# y=0 --> |q| 
# y=1 --> The effective potential (BLPS) 
#----------------------------------------------------------
raw_data_file='BLPS.txt'
vraw=getxyfromfile1(raw_data_file)
vraw.sort()

#----------------------------------------------------------
# (3.2) Initialize qraw (|q|), veff_qraw (BLPS), and
# vnc_qraw (non-Coulombic BLPS)
#----------------------------------------------------------
qraw, veff_qraw, vnc_qraw = [0.0], [veff0], [veff0]
for i in range(0, len(vraw)):
	if vraw[i][0]>tiny and vraw[i][0]<qrawmax and abs(vraw[i][0]-vraw[i-1][0])>tiny:
		qraw.append(vraw[i][0])
		veff_qraw.append(vraw[i][1])
		vnc_qraw.append(vraw[i][1]+4*pi*z/vraw[i][0]**2)
print '(3) Read in qraw, veff_qraw(local pp), vnc_qraw(short range).'


#----------------------------------------------------------
# (4) The first imposed constrain for veff(q):
# if qraw is larger than 10.0, then the Veff is zero
# and vnc has a 4piz/|q|^2 tail.
#----------------------------------------------------------
for i in range(0, len(qraw)):
	if qraw[i]>10.0:
		veff_qraw[i]=0.0
		vnc_qraw[i]=4*pi*z/qraw[i]**2
print '(4) Reset veff_qraw, vnc_qraw if q>10.0, the available q number is ',len(qraw)


#----------------------------------------------------------
# (5) Allocate an uniform q array with dq being the interval.
# We seperate the q-domain into two parts to do spline,
# and only the non-Coulombic part is used.
#----------------------------------------------------------
q = arange(0.0, qmax+dq, dq)
veff, vnc = [], []
print '(5) Set up an uniform q mesh with dimension ',len(q)
Sraw1 = UnivariateSpline(qraw,vnc_qraw,s=1)
Sraw2 = UnivariateSpline(qraw,vnc_qraw,s=2)

for i in range(0,len(q)):
    if q[i]<2.0:
        tmp = Sraw1(q[i])
    else:
        tmp = Sraw2(q[i])

    vnc.append(tmp)

    if q[i]==0.0: 
        veff.append(tmp)
    else:
        veff.append(tmp-4*pi*z/q[i]**2)

veff[0] = veff0
vnc[0] = veff0
print '    Get V(q) on an uniform q grid after performing spline'


#f = open('smooth_BLPS.txt', 'w')
#f.write('q non_Coulombic_BLPS full_BLPS \n')
#for i in range(1, len(veff)):
#    f.write('%15.10f %15.10f %15.10f \n' %(q[i],(veff[i]+4*pi*z/q[i]**2),veff[i]))


#----------------------------------------------------------
# (6) The second constrain imposed on v(q):
# The v(q) must always have the same sign.
# After the zerotail, the veff_qraw is just set to zero
#----------------------------------------------------------
zerotail=0.0
counterzero=0
for i in range(begin_r, len(veff)):
    if veff[i]*veff[i-1]<0:
        counterzero=counterzero+1
        if(counterzero>=1):
            zerotail=q[i]
            break
print '(6) Adjust V(q) from', begin_r*dq, ', V(q) does not change the sign unitl:', zerotail

for i in range(0,len(veff)):
    if q[i]>zerotail:
        veff[i]=0
        vnc[i]=4*pi*z/q[i]**2

print '    Set V(q) to 0 at q larger than ', zerotail
#f = open('vnc.txt', 'w')
#for i in range(0, len(veff)):
#    f.write('%15.10f %15.10f %15.10f \n' %(q[i],veff[i], vnc[i]))


#----------------------------------------------------------
# (7) Setup an uniform grid in real space 
# Perform FFT to get v(r) value.
#----------------------------------------------------------
r = arange(0.0, max_r, dr)
r[0]=1e-6
vrp=[]
print '(7) Set up an uniform real-space grid up to ', max_r, ' with dr=', dr

#----------------------------------------------------------
# (8) We obtain smooth non-Coulombic V(q) in q space,
# and then perform 1D FFT to get
# V(r) = 4*pi* (sum_q sin(qr)*q*v(q)/r)
#----------------------------------------------------------
Sq = splrep(q,veff)
globalr = r[0]
def fun(q):
    sev = splev(q,Sq)
    tmp = sev*sin(globalr*q)*q
    return tmp

#----------------------------------------------------------
# Get the radius r[i] and evaluate the integral:
# Veff(r)=FFT[ Veff(q) ]
# the potential is z/r beyond the radius cutoff.
# the prefactor is 4pi/8pi^3=1/2pi^2
#----------------------------------------------------------
for i in range(0,len(r)):
    globalr=r[i]
    # radius cutoff vr_cutoff,
    if r[i]<vr_cutoff:
        #vtmp=scipy.integrate.quad(fun,0,zerotail,limit=150)
        vtmp=scipy.integrate.quad(fun,0,zerotail)
        vtmp=vtmp/r[i]/2.0/pi/pi
        vrp.append(vtmp[0])
    else:
        vtmp=-z/r[i]
        vrp.append(vtmp)
# vrp is long range with tail -z/r
print '(8) Perform an FFT on V(q)' 


f = open('final_vr.txt', 'w')
f.write('r non-Coulombic_BLPS full_BLPS\n')
for i in range(0, len(r)):
    f.write('%15.10f %15.10f %15.10f\n' %(r[i],vrp[i]+z/r[i],vrp[i]))

#if V(r) is not smooth, turn on the smooth:
#Smr = UnivariateSpline(r,vrp,s=2)
#for i in range(50,len(r)):
#    vrp[i] = Smr(r[i])

#----------------------------------------------------------
# (9) After adjusting V(r) in the real space,
# Now we can perfrom FFT on then new V(r) to obtain final V(q).
#----------------------------------------------------------
Sr = splrep(r,vrp)
globalq = q[0]
def fun2(r):
    sev = splev(r,Sr)
    if globalq == 0:
        tmp = (sev + z/r)*r*r
    else:
        tmp = (sev + z/r)*sin(globalq*r)*r/globalq
    return tmp
print '(9) Tune V(r) by forcing it to be -z/r beyond the cutoff'


#----------------------------------------------------------
# (10) Obtain the final v(q) after smoothing
# the potential in the real space
#----------------------------------------------------------
for i in range(0,len(q)):
    globalq = q[i]
    vvtmp = scipy.integrate.quad(fun2,0,16.0)
    vnc[i] = vvtmp[0]
    vnc[i] = vnc[i]*4.0*pi
    if i==0:
        veff[i] = vnc[i]
    else:
        veff[i] = vnc[i] - 4.*pi*z/globalq**2.
    if q[i] > zerotail: break
print '(10) Obtain the final V(q) after smoothing it in real space'

f = open('final_vq.txt', 'w')
f.write('q non-Coulombic_BLPS full_BLPS\n')
for i in range(0, len(veff)):
    f.write('%15.10f %15.10f %15.10f\n' %(q[i],vnc[i], veff[i]))

#----------------------------------------------------------
# output the pseudopotentials in real space,
# which is used by ABINIT
#----------------------------------------------------------
fr = open('psp.rlpot','w')
fr.write('BLPS | realspace | Coulombic Tail @ %2.1f | Vg0 = %3.2f | VgCut @ %3.2f | %s - %s \n' %(vr_cutoff,veff[0],zerotail,time.localtime()[1],time.localtime()[0]))
fr.write('%3.2f   %3.2f    %s    zatom,  zion,  pspd\n' %(atomic_number, z, time.strftime('%Y%m%d', time.localtime())))
fr.write('8    11   0   0   %d   0    pspcode,pspxc,lmax,lloc,mmax,r2well\n' % (nr) )
fr.write('0   -1   0                       rchrg  fchrg   qchrg\n')
fr.write('0    0   0   0   0   nproj\n')
fr.write('0                    extension_switch\n')
fr.write('0\n')
r[0] = 0.0
for i in range(0,len(r)-1):
    fr.write('%i  %10.8f    %15.10f\n' %(i+1, r[i], vrp[i]))
fr.close()


#----------------------------------------------------------
# output the pseudopotentials in reciprocal space,
# which is used by PROFESS
#----------------------------------------------------------
f = open('psp.recpot', 'w')
f.write('START COMMENT \n')
f.write('BLPS: generated on %s \n' %(time.strftime('%Y%m%d', time.localtime())))
f.write('Rho from ABINIT computations\n')
f.write('END COMMENT \n')
f.write('3     5 \n')
f.write('%15.10f \n' %(q[len(q)-1]/bohr_to_a))

for i in range(0, len(veff)):
    f.write('    %15.10f    ' %(veff[i]*mev_ang3))
    if mod(i+1, 3) == 0 and (i+1 != len(veff)-1): f.write('\n')
f.write('\n1000')
f.close()
lines = float(os.popen("cat psp.recpot | wc -l").read().replace('\n', ''))

print ''
print 'All done. To plot the BLPS, please use the files: final_vq.txt and final_vr.txt'
print 'Have a great day!'
