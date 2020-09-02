#!/bin/bash

shopt -s nullglob

resultFold='resultData/'
#rm -r 'resultData/'i
mkdir -p 'resultData/'

cat $2> ${resultFold}/'runsettings.txt'
ls "$1"/*.{fasta,fna,fa}  >${resultFold}/'fastatorun.txt'

matlab -nodisplay -nosplash -nodesktop -r "addpath(genpath('$1')); addpath(genpath(pwd)); hca_theory_script('resultData/fastatorun.txt', 'resultData/runsettings.txt');quit;";



