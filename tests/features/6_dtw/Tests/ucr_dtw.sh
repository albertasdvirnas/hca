#!/bin/bash


query1=$1
query2=$2

data=$3
M=$4
R=$5
resultFold=$6
result=$7
ucrcode=$8


mkdir -p resultFold

g=$(sed 's/dataSim/resultData/g' <<< "$f")
$8 $query1 $data $M $R > $result
$8 $query2 $data $M $R >> $result

