#!/usr/bin/env bash
A1=6
TH=${OMP_NUM_THREADS:-2}
BV=30
ITER=100
ANG=4
SCH=0

#[SN-KERNEL] Version 1.0 Linear tetrahedral mesh generation...........
#[SN-KERNEL] Usage:   makemesh.x   X    Y    Z  BCmodel...............
#[SN-KERNEL] Example: makemesh.x  10   20   10  0      ...............
#[SN-KERNEL] X, Y, Z specifies the number of X-Y-Z meshes (E=X*Y*Z*12)
#[SN-KERNEL] BC 0-6  specifies how many box surfaces have vacuum b.c.s
echo "[BASH] makemesh.x $A1 $A1 $A1 0"
makemesh.x $A1 $A1 $A1 0 | tee makemesh_2.out
echo "[BASH] cp grid_tet_mesh.ascii inputmesh.ascii"
cp grid_tet_mesh.ascii inputmesh.ascii

#[SN-KERNEL] Version 1.0 mesh processing utility.....................
#[SN-KERNEL] Usage:   processmesh.x  Partition Threads...............
#[SN-KERNEL] Example: processmesh.x  1         32     ...............
#[SN-KERNEL] Scheme   Use (1) METIS, (other) GREEDY..................
#[SN-KERNEL] Threads  Maximum threads to setup for in the mini-app...
echo "[BASH] processmesh.x 1 $TH"
processmesh.x 1 $TH | tee processmesh_2.out

# This tests out the fortran version
echo "[BASH] SNaCFE.x $SCH $ITER $BV $ANG $TH"
SNaCFE.x $SCH $ITER $BV $ANG $TH | tee SNaCFE_2.out

# This tests out the c version
echo "[BASH] cSNaCFE.x $SCH $ITER $BV $ANG $TH"
cSNaCFE.x $SCH $ITER $BV $ANG $TH | tee cSNaCFE_2.out
