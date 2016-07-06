#!/bin/sh
ARCHIVE=loggs

# day=`date "+%Y%m%d"` ## today's date
day=`date "+%Y%m%d%H%M"` ## today's date, hour, minute

ln -s ${ARCHIVE}/${day}_tests.txt ./current
nohup ./tests.py < /dev/null > ${day}_tests.txt
./clearfiles.sh
mv ${day}_tests.txt ${ARCHIVE}
rm -f ./current ## remove old soft link

rm -f ./last ## remove old soft link
ln -s ${ARCHIVE}/${day}_tests.txt ./last

head -5 last
