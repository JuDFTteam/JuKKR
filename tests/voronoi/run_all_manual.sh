
#run_voronoi:intel:oldstyle:
cp -r test_inputs/test01_oldstyle test_run01
cd test_run01 && ln -s ../../../ElementDataBase . && ../../voronoi.exe | tee out_voronoi
cd ../

#run_voronoi:intel:newstyle:
cp -r test_inputs/test02_newstyle test_run02
cd test_run02 && ln -s ../../../ElementDataBase . && ../../voronoi.exe | tee out_voronoi
cd ../
