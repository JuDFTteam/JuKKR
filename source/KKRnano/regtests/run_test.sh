#!/bin/sh
ARCHIVE=loggs

## today's date, hour, minute
day=`date "+%Y%m%d%H%M"` 

## remove old soft link
rm -f ./last

## run the test scripts
nohup python2.7 ./tests.py < /dev/null > ${day}_tests.txt

## clean the working directory
./clearfiles.sh 

## make sure that the archive directory for loggs exists
mkdir -p ${ARCHIVE}

## move the last log file there
mv ${day}_tests.txt ${ARCHIVE}

## set a soft link so it is easier to identify the log file of the last test run 
ln -s ${ARCHIVE}/${day}_tests.txt ./last

## show the first 5 lines, enough to distinguish between success and failure
head -5 last
