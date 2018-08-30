cd tests

echo "###########################################"
echo "run:intel:serial_1:"
echo ""
./run_serial.py 1

echo "###########################################"
echo "run:intel:serial_2:"
echo ""
./run_serial.py 2

echo "###########################################"
echo "run:intel:serial_3:"
echo ""
./run_serial.py 3

echo "###########################################"
echo "run:intel:parallel_2:"
echo ""
./run_parallel.py 2

echo "###########################################"
echo "run:intel:MPIatom_7:"
echo ""
./run_parallel.py -7

echo "###########################################"
echo "run:intel:MPIenerg_8:"
echo ""
./run_parallel.py -8


echo "###########################################"
echo "run:intel:Jijs_4:"
echo ""
./run_serial.py 4
./run_parallel.py -4

echo "###########################################"
echo "run:intel:kkrflex_5:"
echo ""
./run_serial.py 5
./run_parallel.py -5

echo "###########################################"
echo "run:intel:FERMIOUT_6:"
echo ""
./run_serial.py 6
./run_parallel.py -6

echo "###########################################"
echo "run:intel:qdos_12:"
echo ""
./run_parallel.py -12

cd ../

echo "###########################################"
echo "run:intel:OPERATOR_10:"
echo ""
tests/tools/run_test_10.sh

echo "###########################################"
echo "run:intel:DTM_GMATLL_11:"
echo ""
tests/tools/run_test_11.sh


#not working at the moment:
#echo "run:intel:rhoq_13:"
#tests/tools/run_test_13.sh
