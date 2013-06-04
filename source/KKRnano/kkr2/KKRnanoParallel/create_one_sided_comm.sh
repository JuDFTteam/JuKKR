#!/bin/sh

# Create files for different datatypes

sed -e 's/_TYPE/Z/' one_sided_comm_TYPE_mod.F90 > one_sided_commZ_mod.F90
sed -e 's/_TYPE/I/' one_sided_comm_TYPE_mod.F90 > one_sided_commI_mod.F90
#sed -e 's/_TYPE/D/' one_sided_comm_TYPE_mod.F90 > one_sided_commD_mod.F90
#sed -e 's/_TYPE/C/' one_sided_comm_TYPE_mod.F90 > one_sided_commC_mod.F90
