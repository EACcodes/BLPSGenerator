#!/bin/bash
step1=1_OutputDensity
step2=2_BLPS_with_SF

test -d $step2 || mkdir $step2

cp $step1/system.in ./$step2/
cp $step1/refden.in ./$step2/
cp ../job_BLPS.sh ./$step2/

cat > ./$step2/INPUT << EOF
system.in
system.out
i-system
o-system
tmp-system
../../tmp.L0.lps
EOF

cat > ./$step2/param.in << EOF
1       loadDen
0       LoadVr
100     Maxcount
0       bKEDFcc
0       bAtomV
1       bAllowsym
1.0e-6  nontmpTol
BCC     str
1.0e-6  stoptol
5.0e-6  penLambda
1       bspin
EOF
