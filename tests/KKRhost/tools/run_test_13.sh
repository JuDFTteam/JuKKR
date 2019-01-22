#!/usr/bin/env bash

mkdir test_run13
cd test_run13
ln -s ../test_inputs/test_13_*/* .

./prepare.sh
./run.sh

cd ../
