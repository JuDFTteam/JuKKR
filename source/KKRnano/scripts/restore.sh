#!/bin/sh
# Copies files needed to restart from directory given as 1st argument to current dir
cp $1/bin.vpotnew .
cp $1/bin.meshes .
cp $1/bin.meshes.idx .
cp $1/bin.vpotnew.idx .
#cp $1/bin.mesh.* .
#cp $1/new_meshes .
#cp $1/new_meshes.idx .
#cp $1/new_mesh.* .
cp $1/bin.pot.* .
cp $1/bin.energy_mesh .
cp $1/nonco_angle_out.dat .
