#!/bin/bash
#
# shell script used to create data set from id list
# id list is from pdb.org condenced view saved as text file
#
# VR


pdbid=($(awk '{print $1}' id_list.txt))
echo "" > Kd_table.csv
for id in "${pdbid[@]}"
do
echo $id
Kd=$(./get_ligand_Kd.pl $id)
echo $Kd
if [[ $Kd ]]
then
 echo $Kd >> Kd_table.csv
 wget 'http://pdb.org/pdb/files/'$id.pdb.gz
fi
done
