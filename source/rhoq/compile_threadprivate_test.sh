ifort -openmp -r8 -CB -check all -check uninit -ftrapuv -gen-interfaces -warn all -warn notruncated_source -fpe0 -debug extended -traceback -g test_threadprivate.f90
