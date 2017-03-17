#!/bin/bash

# you can put many different vq0 here
for vq0 in 20.0; do
# you can put many rcut here
	for rcut in 5.0; do
      python fitting_BLPS.py 11 1 $vq0 $rcut 
      mv psp.rlpot psp-$vq0-$rcut.rlpot 
	  mv psp.recpot psp-$vq0-$rcut.recpot
    done
done
