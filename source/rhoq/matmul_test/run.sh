rm a.out
ifort -v -O3 -mkl -xhost -vec-report1 test.f90 ; ./a.out | tee out
