#/bin/bash
#
# train VR-score
#
# VR

while read line; do 
    #echo $line
    id=$(echo $line | awk -F ";" '{print $1}')
    echo -n "$id"
    perl -I/home/vytautas/bin/VR_score/ ~/bin/VR_score/VR_score.pl ~/training_set/p_${id}_c.mol2 ~/training_set/l_${id}_c.mol2 ::0
    Kd=$(echo $line | awk -F ";" '{print $4}')
    echo "Kd $Kd"
done < Final_table3.csv
