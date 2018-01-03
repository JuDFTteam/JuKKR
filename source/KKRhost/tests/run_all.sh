mkdir test_run1 test_run2 test_run3

cd test_run1/
ln -s ../test_inputs/test_1_*/* ../../kkr.x .; ./kkr.x | tee out_kkr
rm -f gmat tmat gref *for* inputcard_generated.txt

cd ../test_run2/
ln -s ../test_inputs/test_2_*/* ../../kkr.x .; ./kkr.x | tee out_kkr
rm -f gmat tmat gref *for* inputcard_generated.txt

cd ../test_run3
ln -s ../test_inputs/test_3_*/* ../../kkr.x .; ./kkr.x | tee out_kkr
rm -f gmat tmat gref *for* inputcard_generated.txt
