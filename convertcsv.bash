#!/bin/bash

file=$1
tsvfile=${file%.*}.tsv
echo Creating tab delimited file: $tsvfile
cat $file | sed 's/,/\t/g' > $tsvfile 
cat $tsvfile | awk '{print $1"\t"$2"\t"$3}' > tmp.tsv
tail -n +2 tmp.tsv > $tsvfile
rm tmp.tsv
cat $tsvfile
