#!/bin/bash
MATLAB=/Applications/MATLAB_R2018b.app/bin/matlab

for hh in {10..20..5}
do 
$MATLAB -noFigureWindows -batch "coupled_rf_transp($hh)" &
#echo "$hh"
done

