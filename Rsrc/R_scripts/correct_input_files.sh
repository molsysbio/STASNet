#!/bin/bash
#-*- coding:utf-8 -*-

# This script corrects linebreaks mistakes and aliases the names of the stimulators and inhibitors

if [ $# -lt 1 ]
then
    echo "You need to give at least one file that has to be corrected"
    exit -1
fi

while [[ -e $1 ]]
do
    new_file=$1.corrected
    touch $new_file
    sed s/\\r/\\n/g $1 > $new_file

    # Stimulator/Inhibitor alias
    sed -i s/\".*BSA\"/BSA/g $new_file
    sed -i s/IGF\ 1/IGF/g $new_file
    sed -i s/TGFa/TGFA/g $new_file
    sed -i s/U0126/MEK/g $new_file
    sed -i s/AZD\ 6244/MEK/g $new_file
    sed -i s/LY294002/PI3K/g $new_file
    sed -i s/None//g $new_file
    sed -i s/SB216763/GSK3AB/g $new_file

    # Protein alias, standardizes the notations
    sed -i s/[Ii][Gg][Ff]-*IR/IGFR/g $new_file
    sed -i s/[Ii][Rr][Ss]-*1/IRS1/g $new_file
    sed -i s/[Gg][Ss][Kk]3-*[aA]\\/*[Bb]/GSK3AB/g $new_file
    sed -i s/[Pp]70-*[Ss]6[Kk]/P70S6K/g $new_file

    sed -ie -f ~/bin/sed_upper $new_file # Uppercase
    # For the interpretation in the R script
    sed -i s/MEDIUM/Medium/g $new_file
    sed -i s/BLANK/blank/g $new_file
    sed -i s/\\tT\\t/\\tt\\t/g $new_file
    sed -i s/\\tC\\t/\\tc\\t/g $new_file
    sed -i s/\\tM\\t/\\tm\\t/g $new_file
    sed -i s/INHIBITOR/inhibitor/g $new_file
    sed -i s/STIMULATOR/stimulator/g $new_file
    sed -i s/TYPE/type/g $new_file

    if [ $false ]
    then
    while read line
    do
        insert=true
        for word in line
        do
            if [[ $word == "m" ]]
            then
                insert=false
            fi
        done
        if [ $insert ]
        then
            $line > $new_file
        fi
    done < $new_file
    fi

    shift
done




