# The Jülich KKR codes

## Description

The Korringa-Kohn-Rostoker (KKR) Greens function method is a highly accurate all-electron method to perform density functional theory calculations. The most important features of the Jülich KKR codes include the possibility to perform relativistic calculations, predict scattering effects, and treat finite-sized clusters or very large systems.

## Installation

### Dependencies
- a Fortran compiler (tested with `gfortran` and `ifort` but `ifort` is recommended)
- [cmake](https://cmake.org)

### Compiling the code

The easiest way to set up the code is to execute the `install.py` script.

Alternatively you can configure the build manually as shown below. 

```
mkdir build
cd build
FC=<compiler-you-want-to-use> cmake -D<options> ..
```

where with `FC` you coose the compiler, e.g. `FC=gfortran`, `FC=mpif90` (gfortran with MPI), `FC=ifort`, `FC=mpiifort`, ...

In case cmake does not find the compiler you have chosen try wrapping the compiler like this (here for `mpiifort` as an example): `FC=$(which mpiifort)`

List of default values for `-D<options>` (used if not specified):
```
 -DENABLE_MPI=ON # use MPI parallelization by default
 -DENABLE_OMP=OFF # OpenMP level of parallelization turned off
 -DENABLE_COV=OFF # do not write coverage reports while running (slows code down)
 -DCMAKE_BUILD_TYPE=Release # use release version (alternative: Debug)
 -DENABLE_BdG=OFF # disable Bogoliubov-de-Gennes formalism
```

## Further reading

- The code's [wiki page](https://iffwiki.fz-juelich.de/kkr/doku.php)

- The [source code documentation](https://kkr.iffgit.fz-juelich.de/jukkr/)


## Found a bug?

If you find any bugs, please file a new issue on the [gitlab page](https://iffgit.fz-juelich.de/kkr/jukkr/issues)


