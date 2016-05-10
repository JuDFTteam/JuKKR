#!/bin/sh
# Copies files needed to restart from directory given as 1st argument to current dir
cp $1/vpotnew .
cp $1/meshes .
cp $1/meshes.idx .
cp $1/vpotnew.idx .
cp $1/mesh.* .
cp $1/pot.* .
cp $1/energy_mesh .
