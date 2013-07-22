#/bin/bash
#
# train VR-score
#
# VR

while read line; do 
    #echo $line
    id=$(echo $line | awk -F "," '{print $1}')
    if [ -s ~/training_set_hCAII/p_${id}_c.mol2 ] && [ -s ~/training_set_hCAII/l_${id}_c.mol2 ]
    then
        echo -n "$id"
        if perl -I/home/vytautas/bin/VR_score/ ~/bin/VR_score/VR_score.pl ~/training_set_hCAII/p_${id}_c.mol2 ~/training_set_hCAII/l_${id}_c.mol2 ::0
        then
            Kd=$(echo $line | awk -F "," '{print $4}')
            echo "$Kd"
        fi
    fi
done < Final_table3.csv
