#!/bin/bash


#gunzip *.gz
#rm Kd_ligand_table.csv


for file in ????.pdb
do
Kd=$(./get_ligand_Kd.pl ${file%.pdb})
#echo $Kd >> Kd_ligand_table.csv
ligand=$(echo $Kd | cut -d ',' -f2)
echo $ligand
cat $file | grep ^ATOM > p_$file
cat $file | grep ^HETATM | grep $ligand > l_$file
done
