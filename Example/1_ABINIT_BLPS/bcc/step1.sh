#!/bin/bash
nz=11 # nuclear charge
a0=4.1943 # bcc lattice constant 
pseudo="../../../11-Na.GGA.fhi"
step1=1_OutputDensity

test -d $step1 || mkdir $step1

cp ../job.sh ./$step1/

cat > $step1/INPUT << EOF
system.in
system.out
i-system
o-system
tmp-system
$pseudo
EOF

cat > $step1/system.in << EOF
#######################
# ABINIT INPUT FILE
#######################

#----------------------
# Basic Information
#----------------------
ecut    900 eV # energy cutoff
ixc     11     # PBE 
kptopt  1      # type of k-points
ngkpt   20 20 20  # number of kpoints
nstep   50     # maximal number of SCF cycles
tolvrs  1.0d-8 # tolerence for potential
iscf    7      # Pulay mixing of potential
diemac  1000.0

#----------------------
# ATOM INFORMATION
#----------------------
ntypat  1      # number of types of atom
znucl   $nz      # nuclear charge
nband   4      # number of bands
natom   1      # total number of atoms
typat   1*1    # atom number for each type * type index
xangst 0.0 0.0 0.0 # atom positions (length in Angstrom)

#----------------------
# Data Sets
#----------------------
ndtset 1      # number of data sets
acell $a0 $a0 $a0 Angstrom  
rprim          # lattice vectors
-0.5 0.5 0.5
0.5 -0.5 0.5
0.5 0.5 -0.5

#----------------------
# Occupation, Output
#----------------------
occopt  3 # Fermi-Dirac distributions
tsmear  0.1 eV
chkprim 0 # don't need to check whether the cell is primitive cell
prtden  1 # need to output the charge density

# END #
EOF
