rm a.out
echo 'compile'
ifort -O3 -mkl -openmp test2.f90

export OMP_STACKSIZE=1g
export OMP_NESTED=false

echo
echo '1'
export OMP_NUM_THREADS=1;  ./a.out >  out
echo '2'
export OMP_NUM_THREADS=2;  ./a.out >> out
#echo '3'
#export OMP_NUM_THREADS=3;  ./a.out >> out
echo '4'
export OMP_NUM_THREADS=4;  ./a.out >> out
echo '5'
export OMP_NUM_THREADS=5;  ./a.out >> out
#echo '6'
#export OMP_NUM_THREADS=6;  ./a.out >> out
#echo '7'
#export OMP_NUM_THREADS=7;  ./a.out >> out
#echo '8'
#export OMP_NUM_THREADS=8;  ./a.out >> out
#echo '9'
#export OMP_NUM_THREADS=9;  ./a.out >> out
echo '10'
export OMP_NUM_THREADS=10; ./a.out >> out
#echo '11'
#export OMP_NUM_THREADS=11; ./a.out >> out
#echo '12'
#export OMP_NUM_THREADS=12; ./a.out >> out
#echo '13'
#export OMP_NUM_THREADS=13; ./a.out >> out
#echo '14'
#export OMP_NUM_THREADS=14; ./a.out >> out
echo '15'
export OMP_NUM_THREADS=15; ./a.out >> out
#echo '16'
#export OMP_NUM_THREADS=16; ./a.out >> out
#echo '17'
#export OMP_NUM_THREADS=17; ./a.out >> out
#echo '18'
#export OMP_NUM_THREADS=18; ./a.out >> out
#echo '19'
#export OMP_NUM_THREADS=19; ./a.out >> out
echo '20'
export OMP_NUM_THREADS=20; ./a.out >> out


cat out
