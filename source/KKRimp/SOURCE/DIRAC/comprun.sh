clear
cd /work/iff_th1/pkordt/kkrflex/SOURCE/
echo "compiling..."
# make clean
make
cp kkrflex.exe ../1Fe_0Cu_sratest2/
cd ../1Fe_0Cu_sratest2/
ulimit -s unlimited
ulimit -v unlimited
echo "executing program kkrflex.exe..."
time kkrflex.exe
cd ../SOURCE/DIRAC
./plotR.py
exit 0