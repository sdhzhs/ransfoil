#!/bin/sh
mkdir naca0012_xyz
ransfoil --script test/ransfoil.config.0012.xyz
python comgolden.py naca0012_xyz golden/naca0012_xyz
mkdir naca001264_cpt
ransfoil --script test/ransfoil.config.001264.cpt
python comgolden.py naca001264_cpt golden/naca001264_cpt
mkdir GA_W-1_xyz
ransfoil --script test/ransfoil.config.GA_W-1.xyz
python comgolden.py GA_W-1_xyz golden/GA_W-1_xyz
mkdir GA_W-1_cpt
ransfoil --script test/ransfoil.config.GA_W-1.cpt
python comgolden.py GA_W-1_cpt golden/GA_W-1_cpt
mkdir whitcomb_xyz
ransfoil --script test/ransfoil.config.whitcomb.xyz.1
ransfoil --script test/ransfoil.config.whitcomb.xyz.2
python comgolden.py whitcomb_xyz golden/whitcomb_xyz
mkdir whitcomb_cpt
ransfoil --script test/ransfoil.config.whitcomb.cpt.1
ransfoil --script test/ransfoil.config.whitcomb.cpt.2
python comgolden.py whitcomb_cpt golden/whitcomb_cpt
mkdir goe495_cpt
ransfoil --script test/ransfoil.config.goe495.cpt
python comgolden.py goe495_cpt golden/goe495_cpt
cd scripts
./runwhitcomb.bat
python ../comgolden.py whitcomb1 ../golden/whitcomb_xyz_super
cd ../src
gfcompile.sh
./caller.exe
gccompile.sh
./caller.exe
cd ..
