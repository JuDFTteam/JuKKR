##########################################################################################
# Run auto tests manually, to do verify step go to 'tests' firectory and run 'pytest -v'
##########################################################################################


##########################################################################################
# compile code
cd SOURCE
#make clean; make hybrid; cp kkrflex.exe kkrflex.exe_hybrid
#make clean; make mpi; cp kkrflex.exe kkrflex.exe_mpi
cd ..

##########################################################################################
# then run checks with old solver to check parallelism

if [[ ! -d tests/test_case_kkrflex_host_in_host/imp/test_run_mpi_1 ]]; then
  echo "test 1"
  cd tests/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template/ test_run_mpi_1
  cd test_run_mpi_1/
  ln -s ../../../../SOURCE/kkrflex.exe_mpi kkrflex.exe
  ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -d tests/test_case_kkrflex_host_in_host/imp/test_run_mpi_8 ]]; then
  echo "test 2"
  cd tests/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template/ test_run_mpi_8
  cd test_run_mpi_8/
  ln -s ../../../../SOURCE/kkrflex.exe_mpi kkrflex.exe
  mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -d tests/test_case_kkrflex_host_in_host/imp/test_run_hybrid_1_1 ]]; then
  echo "test 3"
  cd tests/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template/ test_run_hybrid_1_1
  cd test_run_hybrid_1_1/
  ln -s ../../../../SOURCE/kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 1 ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -d tests/test_case_kkrflex_host_in_host/imp/test_run_hybrid_1_8 ]]; then
  echo "test 4"
  cd tests/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template/ test_run_hybrid_1_8
  cd test_run_hybrid_1_8/
  ln -s ../../../../SOURCE/kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -d tests/test_case_kkrflex_host_in_host/imp/test_run_hybrid_8_1 ]]; then
  echo "test 5"
  cd tests/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template/ test_run_hybrid_8_1
  cd test_run_hybrid_8_1/
  ln -s ../../../../SOURCE/kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=8; mpirun -np 1 ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -d tests/test_case_kkrflex_host_in_host/imp/test_run_hybrid_2_4 ]]; then
  echo "test 6"
  cd tests/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template/ test_run_hybrid_2_4
  cd test_run_hybrid_2_4/
  ln -s ../../../../SOURCE/kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=2; mpirun -np 4 ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -d tests/test_case_kkrflex_host_in_host/imp/test_run_hybrid_4_2 ]]; then
  echo "test 7"
  cd tests/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template/ test_run_hybrid_4_2
  cd test_run_hybrid_4_2/
  ln -s ../../../../SOURCE/kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=4; mpirun -np 2 ./kkrflex.exe | tee out
  cd ../../../../
fi

##########################################################################################
# now run tests for features of the code (SOC+mpi, SOC+hybrid, Jijs, wf-saving, sratrick)

if [[ ! -d tests/test_case_kkrflex_host_in_host/imp/test_run_tmatnew_hybrid_2_4 ]]; then
  echo "test 8"
  cd tests/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template_tmatnew/ test_run_tmatnew_hybrid_2_4
  cd test_run_tmatnew_hybrid_2_4/
  ln -s ../../../../SOURCE/kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=2; mpirun -np 4 ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -d tests/test_case_kkrflex_host_in_host/imp/test_run_tmatnew_mpi_8 ]]; then
  echo "test 9"
  cd tests/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template_tmatnew/ test_run_tmatnew_mpi_8
  cd test_run_tmatnew_mpi_8/
  ln -s ../../../../SOURCE/kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -f tests/test_case_kkrflex_host_in_host/imp/test_run_nosavewf/out ]]; then
  echo "test 10"
  cd tests/test_case_kkrflex_host_in_host/imp/test_run_nosavewf/
  ln -s ../../../../SOURCE/kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi
  
if [[ ! -f tests/test_case_kkrflex_host_in_host/imp/test_run_nosratrick/out ]]; then
  echo "test 11"
  cd tests/test_case_kkrflex_host_in_host/imp/test_run_nosratrick/
  ln -s ../../../../SOURCE/kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi
  
if [[ ! -f tests/test_case_kkrflex_host_in_host/imp/test_run_Jij/out ]]; then
  echo "test 12"
  cd tests/test_case_kkrflex_host_in_host/imp/test_run_Jij/
  ln -s ../../../../SOURCE/kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi
  
if [[ ! -f tests/test_case_kkrflex_host_in_host/imp/test_run_Jij_savewf/out ]]; then
  echo "test 13"
  cd tests/test_case_kkrflex_host_in_host/imp/test_run_Jij_savewf/
  ln -s ../../../../SOURCE/kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi
  
if [[ ! -f tests/test_case_kkrflex_host_in_host/imp/test_run_Jij_nosratrick/out ]]; then
  echo "test 14"
  cd tests/test_case_kkrflex_host_in_host/imp/test_run_Jij_nosratrick/
  ln -s ../../../../SOURCE/kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi

##########################################################################################
# finally check outcome of calculations
#
# do this manually: (requires installation of pytest)
#
# cd tests
# pytest -v

