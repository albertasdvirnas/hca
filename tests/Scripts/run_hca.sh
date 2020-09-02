#!/bin/bash

shopt -s nullglob

resultFold='resultData/'
#rm -r 'resultData/'
#mkdir 'resultData/'

cat $2> ${resultFold}/'hcasettings.txt'
ls $1/*.{tif,txt}  >${resultFold}/'tifstorun.txt'
ls $3/*.mat >${resultFold}/'theoryfiles.txt'


matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath('$1')); addpath(genpath(pwd)); hca_script('resultData/tifstorun.txt', 'resultData/hcasettings.txt','resultData/theoryfiles.txt');quit;";


