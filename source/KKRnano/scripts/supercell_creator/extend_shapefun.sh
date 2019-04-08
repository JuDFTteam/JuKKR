#!/bin/bash

IN=shapefun
OUT=shapefun_extended

if [[ $1 -gt 0 ]] ; then TIMES=${1} ; else echo "enter integer number (number of shapefun copies) as first parameter... exit!" ; exit ; fi
if [[ $2 -gt 0 ]] ; then HEAD=${2} ; else echo "enter integer number (number of header lines) as second parameter... exit!" ; exit ; fi

nshin=`head -n 1 $IN`
nshout=`echo "${nshin}*$TIMES" | bc`
nline=`echo "($nshout-1)/4+1" | bc`

echo $nshin $nshout $nline

printf "%5d\n" $nshout > $OUT

for (( i=1; i<=$nline; i++ )) ; do
  echo "  0.100000000000E+01  0.100000000000E+01  0.100000000000E+01  0.100000000000E+01" >> $OUT
done

for (( i=1; i<=$TIMES; i++ )) ; do
    tail -n +$(( $HEAD + 1 )) $IN >> $OUT    
done

exit

#for (( i=1; i<=$TIMES; i++ )) ; do
#    echo "${IN}"
#done | xargs cat > "${OUT}"


