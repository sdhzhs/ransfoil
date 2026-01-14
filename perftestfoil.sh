#!/bin/sh

if [ -z "$1" ]; then
  echo "Error: Missing parameter!"
  echo "Usage: perftestfoil.sh <profile_directory>"
  echo "Example: perftestfoil.sh results/profile"
  exit 1
fi

runIncomCase() {
 mkdir -p $1
 ransfoil --script test/ransfoil.config.$2 > $1/output.txt
 python comperf.py $1/output.txt $3/$1/output.txt >> $1/output.txt
 if [ $? -eq 0 ]; then
   echo "$1 case: Pass"
 else
   echo "$1 case: Fail"
 fi
}

runMeshCase() {
 mkdir -p $1
 ransfoil --mesh test/ransfoil.config.$2 > $1/output.txt
 python comperf.py $1/output.txt $3/$1/output.txt >> $1/output.txt
 if [ $? -eq 0 ]; then
   echo "$1 case: Pass"
 else
   echo "$1 case: Fail"
 fi
}

runComCase() {
 mkdir -p $1
 ransfoil --script test/ransfoil.config.$2.1 > $1/output.txt
 ransfoil --script test/ransfoil.config.$2.2 >> $1/output.txt
 python comperf.py $1/output.txt $3/$1/output.txt >> $1/output.txt
 if [ $? -eq 0 ]; then
   echo "$1 case: Pass"
 else
   echo "$1 case: Fail"
 fi
}

echo "=========================================="
echo "Performance Regression Test"
echo "=========================================="

runIncomCase naca0012_xyz 0012.xyz "$1"
runIncomCase naca001264_cpt 001264.cpt "$1"
runIncomCase GA_W-1_xyz GA_W-1.xyz "$1"
runIncomCase GA_W-1_cpt GA_W-1.cpt "$1"
runComCase whitcomb_xyz whitcomb.xyz "$1"
runComCase whitcomb_cpt whitcomb.cpt "$1"
runIncomCase goe495_cpt goe495.cpt "$1"
runMeshCase cylinder_grd cylinder.grd "$1"
runIncomCase cylinder_xyz cylinder.xyz "$1"
runIncomCase nasasup5_otype_cpt nasasup5.otype.cpt "$1"
runIncomCase naca0012_mat_xyz 0012.mat.xyz "$1"
runIncomCase naca0012_far_cpt 0012.farboundist.cpt "$1"
runIncomCase naca0012_freebc_cpt 0012.freebc.cpt "$1"
runComCase whitcomb_xyz_super whitcomb.xyz.super "$1"

echo "=========================================="
echo "Performance Regression Test Complete"
echo "=========================================="
