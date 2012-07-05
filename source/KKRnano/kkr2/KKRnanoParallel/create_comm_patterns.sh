#!/bin/sh

# Create files for different datatypes

sed -e 's/_TYPE/Z/' comm_patterns_TYPE_mod.F90.template > comm_patternsZ_mod.F90
sed -e 's/_TYPE/D/' comm_patterns_TYPE_mod.F90.template > comm_patternsD_mod.F90
sed -e 's/_TYPE/I/' comm_patterns_TYPE_mod.F90.template > comm_patternsI_mod.F90
sed -e 's/_TYPE/C/' comm_patterns_TYPE_mod.F90.template > comm_patternsC_mod.F90
