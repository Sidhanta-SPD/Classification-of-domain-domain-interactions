#!/bin/bash
for i in $(less uniprots)
do
wget  -O- 'https://www.uniprot.org/uniprot/'$i'' >> database.fsa
#echo $i
done
