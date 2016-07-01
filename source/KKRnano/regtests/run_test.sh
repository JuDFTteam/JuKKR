#!/bin/sh
# day=`date "+%Y%m%d"` ## today's date
day=`date "+%Y%m%d%H%M"` ## today's date, hour, minute
nohup ./tests.py < /dev/null > ${day}_tests.txt 
./clearfiles.sh
head -5 ${day}_tests.txt
