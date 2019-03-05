#!/bin/bash



### README
#
# This is a slurm script (e.g. on CLAIX2018) that submit the tests which are perfomed automatically on a cluster.
# Note that the path to srcdir has to be changed accordingly to run the tests.
# The tests create these directories: build build_1 build_2 build_3 venv_pytest
# The output of the tests can be found in tests_OUTPUT
#
### Now script starts:




### Job name
#SBATCH --job-name=tests_JuKKR

### File for the output
#SBATCH --output=tests_OUTPUT

### Time your job needs to execute, e. g. 50 min
#SBATCH --time=04:00:00

### Use more than one node for parallel jobs on distributed-memory systems, e. g. 2
#SBATCH --nodes=1

### Number of CPUS per task (for distributed-memory parallelisation, use 1)
#SBATCH --cpus-per-task=1

### Disable hyperthreading by setting the tasks per core to 1
#SBATCH --ntasks-per-core=1

### Number of processes per node, e. g. 6 (6 processes on 2 nodes = 12 processes in total)
#SBATCH --ntasks-per-node=8 

### The last part consists of regular shell commands:
### Set the number of threads in your cluster environment to 1, as specified above
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

### Change to working directory
srcdir=/Users/ruess/sourcecodes/JuKKR/
cd $srcdir

# needed for ifflsurm:
export FI_PROVIDER=tcp

### compile KKRhost code and prepare test directory
./install.py --program=kkrhost --compiler=mpiifort --parallelization=hybrid \
  && cd build \
  && make -j 8 \
  && cp -r ../tests/KKRhost/ . \
  && cd KKRhost \
  && ./tools/run_all_manual.sh

### Now compile and run tests for KKRimp 
cd $srcdir
./install.py --program=kkrimp --compiler=mpiifort --parallelization=hybrid \
  && cd build \
  && make -j 8 \
  && ln -s kkrflex.exe kkrflex.exe_mpi && ln -s kkrflex.exe kkrflex.exe_hybrid \
  && cp -r ../tests/KKRimp/ . \
  && cd KKRimp \
  && ./kkrimp_run_tests_manual.sh

### Now compile and run tests for voronoi 
cd $srcdir
./install.py --program=voronoi --compiler=mpiifort --parallelization=hybrid \
  && cd build \
  && make -j 8 \
  && cp -r ../tests/voronoi/ . \
  && cd voronoi \
  && ./run_all_manual.sh

### Now compile and run tests for PKKprime 
cd $srcdir
./install.py --program=pkkprime --compiler=mpiifort --parallelization=hybrid \
  && cd build \
  && make -j 8


## Finally run pytest with output
cd $srcdir
if [[ ! -d venv_pytest ]]; then
  mkdir venv_pytest
  virtualenv-2.7 venv_pytest
fi
source venv_pytest/bin/activate
pip install numpy
pip install pytest

cd $srcdir/build_1/KKRhost
pytest -v --ignore=tools/aiida_simple_test.py 

cd $srcdir/build_2/KKRimp
pytest -v
