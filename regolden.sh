#!/bin/sh
# Usage: [SRC_DIR=<results dir>] ./regolden.sh <case1> [case2 ...]
#        [SRC_DIR=<results dir>] ./regolden.sh all
# Special case name "libcall" regoldens src/output.txt.
# SRC_DIR (default: current dir) is the directory holding the new case outputs.

SRC_DIR=${SRC_DIR:-.}

ALL_CASES="naca0012_xyz naca001264_cpt GA_W-1_xyz GA_W-1_cpt whitcomb_xyz whitcomb_cpt
goe495_cpt cylinder_grd cylinder_xyz nasasup5_otype_cpt naca0012_mat_xyz naca0012_far_cpt
naca0012_freebc_cpt whitcomb_xyz_super naca0012_pseudo_xyz naca001264_pseudo_local_cpt
libcall"

if [ $# -eq 0 ]; then
    echo "Usage: $0 <case1> [case2 ...] | all" >&2
    exit 1
fi

if [ $# -eq 1 ] && [ "$1" = "all" ]; then
    set -- $ALL_CASES
fi

status=0
for c in "$@"; do
    if [ "$c" = "libcall" ]; then
        cp "$SRC_DIR/src/output.txt" golden/ || status=1
        continue
    fi
    if [ ! -d "$SRC_DIR/$c" ]; then
        echo "Skipping $c: $SRC_DIR/$c not found" >&2
        status=1
        continue
    fi
    mkdir -p "golden/$c"
    cp "$SRC_DIR/$c/Aeroreport.dat" "$SRC_DIR/$c/Wallsol.dat" "golden/$c/" || status=1
done

echo "Regolden completed."
exit $status