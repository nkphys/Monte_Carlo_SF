#!/bin/bash -x

#(ky_i*Parameters_.lx) + kx_i

L=6
file_in="Akw.txt"
file_out="Akw_arranged.txt"
rm $file_out

total_lines=$(awk 'END{print FNR}' ${file_in})
k_plot=0

##----[\GAMMA TO X]-------
ky_i=0
for kx_i in {0..3}
do
let k_i=${ky_i}*${L}+${kx_i}

COUNTER=1
while [  $COUNTER -le $total_lines ]
 do
      val=$(awk -v row=$COUNTER 'NR==row{print $3}' ${file_in})
#	echo "$val"
     # if [ ${val} -eq ${k_i} ]
#	then
#	line=$(awk -v row="$COUNTER" 'NR==row' $file_in)
#	echo "$k_plot   $line" >> $file_out
#	fi
      let COUNTER=COUNTER+1
      echo $COUNTER
done

#awk -v k=$k_i '$3 == k' $file_in >> $file_out
echo '' >> $file_out
let k_plot=k_plot+1
done
##----------------------

##----(X TO M]-------

##----------------------
