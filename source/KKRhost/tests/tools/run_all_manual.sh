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


echo "###########################################"
echo "run:intel:OPERATOR_10:"
echo ""
cd ../
#tests/tools/run_test_10.sh
cd tests

echo "###########################################"
echo "run:intel:DTM_GMATLL_11:"
echo ""
cd ../
#tests/tools/run_test_11.sh
cd tests

echo "###########################################"
echo "run:intel:rhoq_13:"
echo ""
cd ../
#tests/tools/run_test_13.sh
cd tests

echo "###########################################"
echo "run:intel:ASA_14:"
echo ""
./run_parallel.py -14

echo "###########################################"
echo "run:intel:CPA_15:"
echo ""
./run_parallel.py -15

echo "###########################################"
echo "run:intel:Dirac_16:"
echo ""
./run_parallel.py -16

echo "###########################################"
echo "run:intel:lambda_xc_17:"
echo ""
./run_parallel.py -17

echo "###########################################"
echo "run:intel:noco_18:"
echo ""
./run_parallel.py -18


echo "###########################################"
echo "run:intel:decimate_19:"
echo ""
cd ../
tests/tools/run_test_19.sh
cd tests


#SOC tests

echo "###########################################"
echo "run:intel:Au_bulk_SOC_1_1:"
echo ""
./run_parallel.py -1001

echo "###########################################"
echo "run:intel:Fe_slab_SOC_2_1:"
echo ""
./run_parallel.py -1002

echo "###########################################"
echo "run:intel:Si_LLY_SOC_3_1:"
echo ""
./run_parallel.py -1003

echo "###########################################"
echo "run:intel:Jijs_SOC_4_1:"
echo ""
./run_parallel.py -1004

echo "###########################################"
echo "run:intel:kkrflex_SOC_5_1:"
echo ""
./run_parallel.py -1005

echo "###########################################"
echo "run:intel:FERMIOUT_SOC_6_1:"
echo ""
./run_parallel.py -1006

echo "###########################################"
echo "run:intel:qdos_SOC_12_1:"
./run_parallel.py -1012
echo ""

echo "###########################################"
echo "run:intel:ASA_SOC_14_1:"
./run_parallel.py -1014
echo ""

echo "###########################################"
echo "run:intel:CPA_SOC_15_1:"
./run_parallel.py -1015
echo ""


echo "###########################################"
echo "run:intel:decimate_SOC_19_1:"
echo ""
cd ../
tests/tools/run_test_19_1.sh
cd tests

