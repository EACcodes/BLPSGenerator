#!/bin/bash
step2=2_BLPS_with_SF
step3=3_BLPS_without_SF
code=/home/mohanc/bin/abinip

test -d $step3 || mkdir $step3

cp $step2/INPUT ./$step3/
cp $step2/system.in ./$step3/
cp $step2/res_vion.dat ./$step3/refpot.in
cp ../job_BLPS.sh ./$step3/

cat > ./$step3/param.in << EOF
0       loadDen
1       loadVr
100     Maxcount
0       bKEDFcc
1       bAtomV
1       bAllowsym
1.0e-3  nontmpTol
DIA     str
1.0e-3  stoptol
1.0e-5  penLambda
1       bspin
EOF
