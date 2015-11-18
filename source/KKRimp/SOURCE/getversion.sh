echo "      character(len=*), parameter :: version='" > tmpver
git describe >> tmpver
echo "'" >> tmpver
tr -d '\n' < tmpver > version.f90
echo "" >> version.f90
rm -f tmpver
