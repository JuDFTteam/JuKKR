##########################################################################################
# Run auto tests manually, to do verify step go to 'tests' firectory and run 'pytest -v'
##########################################################################################

##########################################################################################
# then run checks with old solver to check parallelism

if [[ ! -d KKRimp/test_case_kkrflex_host_in_host/imp/test_run_mpi_1 ]]; then
  echo "test 1"
  cd KKRimp/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template/ test_run_mpi_1
  cd test_run_mpi_1/
  ln -s ../../../../kkrflex.exe_mpi kkrflex.exe
  ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -d KKRimp/test_case_kkrflex_host_in_host/imp/test_run_mpi_8 ]]; then
  echo "test 2"
  cd KKRimp/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template/ test_run_mpi_8
  cd test_run_mpi_8/
  ln -s ../../../../kkrflex.exe_mpi kkrflex.exe
  mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -d KKRimp/test_case_kkrflex_host_in_host/imp/test_run_hybrid_1_1 ]]; then
  echo "test 3"
  cd KKRimp/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template/ test_run_hybrid_1_1
  cd test_run_hybrid_1_1/
  ln -s ../../../../kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 1 ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -d KKRimp/test_case_kkrflex_host_in_host/imp/test_run_hybrid_1_8 ]]; then
  echo "test 4"
  cd KKRimp/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template/ test_run_hybrid_1_8
  cd test_run_hybrid_1_8/
  ln -s ../../../../kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -d KKRimp/test_case_kkrflex_host_in_host/imp/test_run_hybrid_8_1 ]]; then
  echo "test 5"
  cd KKRimp/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template/ test_run_hybrid_8_1
  cd test_run_hybrid_8_1/
  ln -s ../../../../kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=8; mpirun -np 1 ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -d KKRimp/test_case_kkrflex_host_in_host/imp/test_run_hybrid_2_4 ]]; then
  echo "test 6"
  cd KKRimp/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template/ test_run_hybrid_2_4
  cd test_run_hybrid_2_4/
  ln -s ../../../../kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=2; mpirun -np 4 ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -d KKRimp/test_case_kkrflex_host_in_host/imp/test_run_hybrid_4_2 ]]; then
  echo "test 7"
  cd KKRimp/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template/ test_run_hybrid_4_2
  cd test_run_hybrid_4_2/
  ln -s ../../../../kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=4; mpirun -np 2 ./kkrflex.exe | tee out
  cd ../../../../
fi

##########################################################################################
# now run tests for features of the code (SOC+mpi, SOC+hybrid, Jijs, wf-saving, sratrick)

if [[ ! -d KKRimp/test_case_kkrflex_host_in_host/imp/test_run_tmatnew_hybrid_2_4 ]]; then
  echo "test 8"
  cd KKRimp/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template_tmatnew/ test_run_tmatnew_hybrid_2_4
  cd test_run_tmatnew_hybrid_2_4/
  ln -s ../../../../kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=2; mpirun -np 4 ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -d KKRimp/test_case_kkrflex_host_in_host/imp/test_run_tmatnew_mpi_8 ]]; then
  echo "test 9"
  cd KKRimp/test_case_kkrflex_host_in_host/imp/
  cp -r test_run_template_tmatnew/ test_run_tmatnew_mpi_8
  cd test_run_tmatnew_mpi_8/
  ln -s ../../../../kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi

if [[ ! -f KKRimp/test_case_kkrflex_host_in_host/imp/test_run_nosavewf/out ]]; then
  echo "test 10"
  cd KKRimp/test_case_kkrflex_host_in_host/imp/test_run_nosavewf/
  ln -s ../../../../kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi
  
if [[ ! -f KKRimp/test_case_kkrflex_host_in_host/imp/test_run_nosratrick/out ]]; then
  echo "test 11"
  cd KKRimp/test_case_kkrflex_host_in_host/imp/test_run_nosratrick/
  ln -s ../../../../kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi
  
if [[ ! -f KKRimp/test_case_kkrflex_host_in_host/imp/test_run_Jij/out ]]; then
  echo "test 12"
  cd KKRimp/test_case_kkrflex_host_in_host/imp/test_run_Jij/
  ln -s ../../../../kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi
  
if [[ ! -f KKRimp/test_case_kkrflex_host_in_host/imp/test_run_Jij_savewf/out ]]; then
  echo "test 13"
  cd KKRimp/test_case_kkrflex_host_in_host/imp/test_run_Jij_savewf/
  ln -s ../../../../kkrflex.exe_hybrid kkrflex.exe
  export OMP_NUM_THREADS=1; mpirun -np 8 ./kkrflex.exe | tee out
  cd ../../../../
fi
  
if [[ ! -f KKRimp/test_case_kkrflex_host_in_host/imp/test_run_Jij_nosratrick/out ]]; then
  echo "test 14"
  cd KKRimp/test_case_kkrflex_host_in_host/imp/test_run_Jij_nosratrick/
  ln -s ../../../../kkrflex.exe_hybrid kkrflex.exe
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

