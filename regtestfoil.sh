#!/bin/sh
runIncomCase() {
 mkdir $1
 ransfoil --script test/ransfoil.config.$2 > $1/output.txt
 python comgolden.py $1 $3/$1 >> $1/output.txt
 if [ $? -eq 0 ]; then
   echo "$1 case: Pass"
 else
   echo "$1 case: Fail"
 fi
}

runMeshCase() {
 mkdir $1
 ransfoil --mesh test/ransfoil.config.$2 > $1/output.txt
 python comgolden.py $1 $3/$1 >> $1/output.txt
 if [ $? -eq 0 ]; then
   echo "$1 case: Pass"
 else
   echo "$1 case: Fail"
 fi
}

runComCase() {
 mkdir $1
 ransfoil --script test/ransfoil.config.$2.1 > $1/output.txt
 ransfoil --script test/ransfoil.config.$2.2 >> $1/output.txt
 python comgolden.py $1 $3/$1 >> $1/output.txt
 if [ $? -eq 0 ]; then
   echo "$1 case: Pass"
 else
   echo "$1 case: Fail"
 fi
}

echo "=========================================="
echo "Normal Accuracy Regression Test"
echo "=========================================="

runIncomCase naca0012_xyz 0012.xyz golden
runIncomCase naca001264_cpt 001264.cpt golden
runIncomCase GA_W-1_xyz GA_W-1.xyz golden
runIncomCase GA_W-1_cpt GA_W-1.cpt golden
runComCase whitcomb_xyz whitcomb.xyz golden
runComCase whitcomb_cpt whitcomb.cpt golden
runIncomCase goe495_cpt goe495.cpt golden
runMeshCase cylinder_grd cylinder.grd golden
runIncomCase cylinder_xyz cylinder.xyz golden
runIncomCase nasasup5_otype_cpt nasasup5.otype.cpt golden
runIncomCase naca0012_mat_xyz 0012.mat.xyz golden
runIncomCase naca0012_far_cpt 0012.farboundist.cpt golden
runIncomCase naca0012_freebc_cpt 0012.freebc.cpt golden
runComCase whitcomb_xyz_super whitcomb.xyz.super golden

cd src
gfcompile.sh
./caller.exe > output.txt
if [ $? -eq 0 ]; then
  echo "fortran call lib case: Pass"
else
  echo "fortran call lib case: Fail"
fi
gccompile.sh
./caller.exe >> output.txt
if [ $? -eq 0 ]; then
  echo "c call lib case: Pass"
else
  echo "c call lib case: Fail"
fi
cd ..

echo "=========================================="
echo "Normal Accuracy Regression Test Complete"
echo "=========================================="
