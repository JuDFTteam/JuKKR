#!/bin/sh

# Create files for different datatypes

sed -e 's/_TYPE/Z/' one_sided_comm_TYPE_mod.F90 > one_sided_commZ_mod.F90
sed -e 's/_TYPE/I/' one_sided_comm_TYPE_mod.F90 > one_sided_commI_mod.F90
sed -e 's/_TYPE/D/' one_sided_comm_TYPE_mod.F90 > one_sided_commD_mod.F90
#sed -e 's/_TYPE/C/' one_sided_comm_TYPE_mod.F90 > one_sided_commC_mod.F90

#sed -e 's/_TYPE/Z/' TruncationZone_TYPE_com_mod.F90 > TruncationZoneZ_com_mod.F90
#sed -e 's/_TYPE/I/' TruncationZone_TYPE_com_mod.F90 > TruncationZoneI_com_mod.F90
#sed -e 's/_TYPE/D/' TruncationZone_TYPE_com_mod.F90 > TruncationZoneD_com_mod.F90
#sed -e 's/_TYPE/C/' TruncationZone_TYPE_com_mod.F90 > TruncationZoneC_com_mod.F90
