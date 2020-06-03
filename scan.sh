#!/bin/bash
MATLAB=/Applications/MATLAB_R2018b.app/bin/matlab

for hh in {10..25..5}
do
for ww in {0.5e17,0.7925e17,1.256e17,1.9907e17,3.155e17,5.0e17}
do 
$MATLAB -noFigureWindows -batch "coupled_rf_transp($hh,$ww)" &
#echo "$hh $ww"
done
done

