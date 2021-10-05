#!/bin/bash
#inputgenerator.py Dimensions dimensions.txt > Dimensions_mod.f90
cp InputParams_mod.F90 InputParams_mod.F90.bak
./inputgenerator.py InputParams InputParamsNew.txt > InputParams_mod.F90
