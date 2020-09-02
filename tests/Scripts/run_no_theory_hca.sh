#!/bin/bash

shopt -s nullglob

resultFold='resultData/'
#rm -r 'resultData/'
mkdir -p 'resultData/'

cat $2> ${resultFold}/'hcasettings.txt'
ls "$1"/*.{tif,txt}  >${resultFold}/'tifstorun.txt'

#matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath('$3')); hca_no_theory('resultData/tifstorun.txt', 'resultData/hcasettings.txt');quit;";


