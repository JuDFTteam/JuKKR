#!/bin/sh

num_atoms_in=8       # number of atoms in the potential to be inflated
inflation_factor=2   # repetition of unit cell in x-,y- and z-direction
sf_header_lines=3    # number of header lines in 'shapefun'

fac_x=${inflation_factor}
fac_y=${inflation_factor}
fac_z=${inflation_factor}

prod_facs=$((${fac_x} * ${fac_y} * ${fac_z})) 
num_atoms_out=$((${prod_facs} * ${num_atoms_in}))

echo "Create calculation data for " ${fac_x} "x" ${fac_y} "x" ${fac_z} " supercell"
python2.7 repeatcell_xyz.py ${num_atoms_out} ${fac_x} ${fac_y} ${fac_z}
./extend_shapefun.sh ${prod_facs} ${sf_header_lines}
mv shapefun_extended newcell/shapefun


echo "Calculation data created for " ${num_atoms_out} " atoms" 

