#!/bin/bash
#   makes animated gif to demonstrate growth in instability of Lax-Wendroff 
#    problem with poor choice of Courant condition, a = 1.1
script='gifplot.gps'

rm -f $script

echo 'set terminal gif animate delay 50' >> $script
echo "set output 'instab.gif'" >> $script
echo "set xrange [0:1]" >> $script
for i in *.dat
do
    echo "plot '$i' using 1:2 with lines" >> $script
done

gnuplot $script
rm $script
