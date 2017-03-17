#!/bin/bash

# nz: the number of valence electron number 
nz=1
# data_dir: the directory where you generate the c_vag.out file
data_dir=3_BLPS_without_SF
# fitting_dir: the new directory
fitting_dir=FITTING

# sum_raw_file: the raw data of non-Coulombic BLPS for you to plot
sum_raw_file="total-raw.txt"
test -e $fitting_dir/$sum_raw_file && rm $fitting_dir/$sum_raw_file

# create a new directory
test -d $fitting_dir || mkdir $fitting_dir

#for dir in fcc sc bcc dia; do
for dir in bcc ; do

  echo Now working on: $dir

  # make sure the old BLPS file doesn't exist
  file=$fitting_dir/$dir-BLPS.txt
  test -e $file && rm $file

  # copy the c_vag.out directory to the current directory
  # the c_vag.out file contains effective potential (BLPS) on discretized |q|
  cp $dir/$data_dir/c_vag.out $file

  # get the sorted raw data of BLPS
  cd $fitting_dir
  python ../sort_Vq.py $dir-BLPS.txt $nz
  mv raw_data.txt $dir-raw.txt

  # collect non-Coulombic BLPS from all structures
  # for you to see clearly the shape of the
  # non-Coulombic BLPS
  cat $dir-raw.txt >> $sum_raw_file 
  cd ..
  
  # collect the discretized BLPS from all structures
  # for you to do spline further
  cat $file >> $fitting_dir/BLPS.txt

done
