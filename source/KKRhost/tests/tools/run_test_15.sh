#!/usr/bin/env bash

mkdir tests/test_run15
cd tests/test_run15
ln -s ../test_inputs/test_15_*/* .

./prepare.sh
./run.sh

