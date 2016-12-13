#!/bin/sh

# Create files for different datatypes

sed -e 's/_TYPE/Z/' two_sided_comm_TYPE_mod.F95 > two_sided_commZ_mod.F90
sed -e 's/_TYPE/I/' two_sided_comm_TYPE_mod.F95 > two_sided_commI_mod.F90
sed -e 's/_TYPE/D/' two_sided_comm_TYPE_mod.F95 > two_sided_commD_mod.F90
# sed -e 's/_TYPE/C/' two_sided_comm_TYPE_mod.F95 > two_sided_commC_mod.F90
