#download aiida-kkr to have tests directory
git clone --depth=1 https://github.com/JuDFTteam/aiida-kkr.git

# fake code install (reuse already compiled code)
mkdir aiida-kkr/aiida_kkr/tests/jukkr
ln -s ../../kkr.x aiida-kkr/aiida_kkr/tests/jukkr/kkr.x
ln -s ../../voronoi.exe aiida-kkr/aiida_kkr/tests/jukkr/voronoi.exe
ln -s ../../../ElementDataBase aiida-kkr/aiida_kkr/tests/jukkr/ElementDataBase


# now run aiida-kkr tests
cd aiida-kkr/aiida_kkr/tests/

# change dbsetup to match slurm settings
sed -i "s/kkr_codename = 'kkrhost'/kkr_codename = 'KKRhost'/g" dbsetup.py 
#sed -i "s/computername = 'localhost'/computername = 'slurmcontrol'/g" dbsetup.py 
#sed -i "s/kkr_codename = 'kkrhost'/kkr_codename = 'KKRcode'/g" dbsetup.py 

# first run tests without doing actual kkr or kkrimp calculations (only voronoi included)
./run_all.sh

# now do kkr calculations
pytest --cov-report=term-missing --cov-append --cov=aiida_kkr --ignore=jukkr -k Test_dos_workflow
pytest --cov-report=term-missing --cov-append --cov=aiida_kkr --ignore=jukkr -k Test_gf_writeout_workflow
pytest --cov-report=term-missing --cov-append --cov=aiida_kkr --ignore=jukkr -k Test_scf_workflow

# and finally kkrimp calculations
#pytest --cov-report=term-missing --cov-append --cov=aiida_kkr --ignore=jukkr -k Test_kkrimp_scf_workflow
#pytest --cov-report=term-missing --cov-append --cov=aiida_kkr --ignore=jukkr -k Test_kkrimp_full_workflow


# then upload test results to codecov
#pip install codecov
#codecov -t ad573476-72fe-49e3-910e-314f40cf6f6e
