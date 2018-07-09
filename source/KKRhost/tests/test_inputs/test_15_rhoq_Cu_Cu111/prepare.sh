#!/usr/bin/env bash

# compile KKRcode in hybrid version
cp inc.p ../../
cd ../../
make clean && make hybrid -j1
cd -
cp ../../kkr.x .
# compile rhoq code
cd rhoq_module/
make
cd ..
