#!/bin/sh
mkdir naca0012_xyz
ransfoil --script test/ransfoil.config.0012.xyz > naca0012_xyz/output.txt
python comgolden.py naca0012_xyz golden/naca0012_xyz >> naca0012_xyz/output.txt
if [ $? -eq 0 ]; then
  echo "naca0012_xyz case: Pass"
else
  echo "naca0012_xyz case: Fail"
fi
mkdir naca001264_cpt
ransfoil --script test/ransfoil.config.001264.cpt > naca001264_cpt/output.txt
python comgolden.py naca001264_cpt golden/naca001264_cpt >> naca001264_cpt/output.txt
if [ $? -eq 0 ]; then
  echo "naca001264_cpt case: Pass"
else
  echo "naca001264_cpt case: Fail"
fi
mkdir GA_W-1_xyz
ransfoil --script test/ransfoil.config.GA_W-1.xyz > GA_W-1_xyz/output.txt
python comgolden.py GA_W-1_xyz golden/GA_W-1_xyz >> GA_W-1_xyz/output.txt
if [ $? -eq 0 ]; then
  echo "GA_W-1_xyz case: Pass"
else
  echo "GA_W-1_xyz case: Fail"
fi
mkdir GA_W-1_cpt
ransfoil --script test/ransfoil.config.GA_W-1.cpt > GA_W-1_cpt/output.txt
python comgolden.py GA_W-1_cpt golden/GA_W-1_cpt >> GA_W-1_cpt/output.txt
if [ $? -eq 0 ]; then
  echo "GA_W-1_cpt case: Pass"
else
  echo "GA_W-1_cpt case: Fail"
fi
mkdir whitcomb_xyz
ransfoil --script test/ransfoil.config.whitcomb.xyz.1 > whitcomb_xyz/output.txt
ransfoil --script test/ransfoil.config.whitcomb.xyz.2 >> whitcomb_xyz/output.txt
python comgolden.py whitcomb_xyz golden/whitcomb_xyz >> whitcomb_xyz/output.txt
if [ $? -eq 0 ]; then
  echo "whitcomb_xyz case: Pass"
else
  echo "whitcomb_xyz case: Fail"
fi
mkdir whitcomb_cpt
ransfoil --script test/ransfoil.config.whitcomb.cpt.1 > whitcomb_cpt/output.txt
ransfoil --script test/ransfoil.config.whitcomb.cpt.2 >> whitcomb_cpt/output.txt
python comgolden.py whitcomb_cpt golden/whitcomb_cpt >> whitcomb_cpt/output.txt
if [ $? -eq 0 ]; then
  echo "whitcomb_cpt case: Pass"
else
  echo "whitcomb_cpt case: Fail"
fi
mkdir goe495_cpt
ransfoil --script test/ransfoil.config.goe495.cpt > goe495_cpt/output.txt
python comgolden.py goe495_cpt golden/goe495_cpt >> goe495_cpt/output.txt
if [ $? -eq 0 ]; then
  echo "goe495_cpt case: Pass"
else
  echo "goe495_cpt case: Fail"
fi
cd scripts
./runwhitcomb.bat
python ../comgolden.py whitcomb1 ../golden/whitcomb_xyz_super > 3.txt
if [ $? -eq 0 ]; then
  echo "whitcomb_xyz_super case: Pass"
else
  echo "whitcomb_xyz_super case: Fail"
fi
cd ../src
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
