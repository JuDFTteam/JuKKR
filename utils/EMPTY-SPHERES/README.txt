Empty sphere finder
===================

This program finds optimal positions for empty spheres in 3D-periodic lattices
by means of a Monte-Carlo optimization.

On input it requires the Bravais vectors, number and positions of basis sites 
of the lattice (NFIX, RFIX), and wished number of empty spheres (NVAR) that it
positions on output at RVAR.

Created by Phivos Mavropoulos, 2013-2014.


Compilation
-----------
For compilation settings see `makefile`.
Compile with `make` in this folder


Example Usage
-------------
For example input see inputcard in this folder


Developer comments
------------------

01.05.2015
Fixed a bug in subr. rationalbasis that created problems when the fixed atoms
were outside the primitive cell defined by the three Bravais vectors.

Changed default of NSUPER from 2 to 1.
