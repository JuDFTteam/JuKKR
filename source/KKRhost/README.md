# Jülich KKR code for bulk and interfaces

## Description

Some short description ...

## Installation

### Dependencies
- Fortran compiler
- cmake

### Compiling the code

```
mkdir build
cd build
FC=<compiler-you-want-to-use> cmake -D<options> ..
```

where with `FC` you coose the compiler, e.g. `FC=gfortran`, `FC=mpif90` (gfortran with MPI), `FC=ifort`, `FC=mpiifort`, ...

List of default values for `-D<options>` (used if not specified):
```
 -DENABLE_MPI=ON
 -DENABLE_OMP=OFF
 -DENABLE_COV=OFF
 -DCMAKE_BUILD_TYPE=Release
 -DENABLE_BdG=OFF
```

## Further reading

- The code's [wiki page](https://iffwiki.fz-juelich.de/kkr/doku.php)

- The [source code documentation](https://kkr.iffgit.fz-juelich.de/kkrjm/)

