#!/usr/bin/env bash

mkdir tests/test_run13
cd tests/test_run13
ln -s ../test_inputs/test_13_*/* .

./prepare.sh
./run.sh

