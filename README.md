# The Jülich KKR codes

## Description

The Korringa-Kohn-Rostoker (KKR) Greens function method is a highly accurate all-electron method to perform density functional theory calculations. The most important features of the Jülich KKR codes include the possibility to perform relativistic calculations, predict scattering effects, and treat finite-sized clusters or very large systems.

## Installation

### Dependencies
- a Fortran compiler (tested with `gfortran` and `ifort` but `ifort` is recommended)
- [cmake](https://cmake.org)
- an installation of [LAPACK](http://www.netlib.org/lapack/)
- a compiler supporting MPI (optional but strongly recommended)

### Compiling the code

The easiest way to set up the code is to execute the `install.py` script which will guide you through the installation. Afterwards you shoud go to the `build` directory and execute `make` which will start the compilation of the code. The compiled executable will then be placed in the `build` directory.

## Further reading

- The code's [wiki page](https://iffgit.fz-juelich.de/kkr/jukkr/wikis/home)

- The [source code documentation](https://kkr.iffgit.fz-juelich.de/jukkr/)

## Development

Code development is done on the main [gitlab server hosted at the FZ Jülich](https://iffgit.fz-juelich.de/kkr/jukkr) (and not on the [mirrored repository on github](https://github.com/JuDFTteam/JuKKR)).

Creating / editing issues, branches etc. is only supported on the [gitlab version hosted at the FZ Jülich](https://iffgit.fz-juelich.de/kkr/jukkr) which requires sign in (possible with free github account).

## Found a bug?

If you find any bugs, please file a new issue on the [gitlab issues page](https://iffgit.fz-juelich.de/kkr/jukkr/issues).


