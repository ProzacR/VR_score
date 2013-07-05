#!/bin/bash


gunzip *.gz
#rm Kd_ligand_table.csv


for file in ????.pdb
do
Kd=$(./get_ligand_Kd.pl ${file%.pdb})
echo $Kd >> Kd_ligand_table.csv
ligand=$(echo $Kd | cut -d ',' -f2)
echo $ligand
cat $file | grep ^ATOM > p_$file
cat $file | grep ^HETATM | grep $ligand > l_$file


echo l_$file >> lig_convert_log.txt
babel -h l_$file l_${file%.pdb}.mol2 &>> lig_convert_log.txt
babel --partialcharge mmff94 l_${file%.pdb}.mol2 l_${file%.pdb}_c.mol2 &>> lig_convert_log.txt


echo p_$file >> prot_convert_log.txt
babel -h p_$file p_${file%.pdb}.mol2 &>> prot_convert_log.txt
babel --partialcharge mmff94 p_${file%.pdb}.mol2 p_${file%.pdb}_c.mol2 &>> prot_convert_log.txt


checkl=$(ls -lh l_${file%.pdb}_c.mol2 | awk '{print $5}')
checkp=$(ls -lh p_${file%.pdb}_c.mol2 | awk '{print $5}')
echo "$file, $checkl, $checkp" >> convert_log.txt

done
