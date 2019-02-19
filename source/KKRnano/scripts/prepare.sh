#!/bin/sh
echo Using calculated potential as new starting potential...
cp bin.vpotnew bin.vpotnew.0
cp bin.vpotnew.idx bin.vpotnew.0.idx
cp bin.meshes bin.meshes.0
cp bin.meshes.idx bin.meshes.0.idx
#cp new_meshes new_meshes.0
#cp new_meshes.idx new_meshes.0.idx
cp bin.energy_mesh bin.energy_mesh.0
cp nonco_angle_out.dat nonco_angle.dat
