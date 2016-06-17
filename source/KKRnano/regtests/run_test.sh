#!/bin/sh
day=`date "+%Y%m%d"` ## today's date
nohup ./tests.py < /dev/null > ${day}_tests.txt 
./clearfiles.sh
head -5 ${day}_tests.txt
