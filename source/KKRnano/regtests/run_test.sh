#!/bin/sh
ARCHIVE=loggs

# day=`date "+%Y%m%d"` ## today's date
day=`date "+%Y%m%d%H%M"` ## today's date, hour, minute


rm -f ./last ## remove old soft link

nohup ./tests.py < /dev/null > ${day}_tests.txt

./clearfiles.sh
mv ${day}_tests.txt ${ARCHIVE}

ln -s ${ARCHIVE}/${day}_tests.txt ./last

head -5 last
