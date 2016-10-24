echo 'compile'
cd ../
#make clean
#make serial
#make explicit1
#make explicit2

cd -
echo 'run tests'

../rhoq.x; cp out_timing.000.txt time_imp
../rhoq.x_expl1; cp out_timing.000.txt time_exp1
../rhoq.x_expl2; cp out_timing.000.txt time_exp2
