#!/bin/bash
#-*- coding:utf8 -*-

# Create subfolders for the cell lines data in the argument folder
if [[ $# == 0 ]]
then
    echo "This command needs a folder where to look"
    exit 1
elif [[ $# > 1 ]]
then
    echo "This command only takes one argument, the others will be ignored"
fi

file=$1
if [[ ! $1 =~ "/$" ]]
then
    file=$1/
fi

for cell_line in `ls $file*.data`
do
    base=`echo $cell_line | sed -r s_/.+/\([^/]+\)\\.data_\\\\1_`
    if [[ ! -e $base ]]
    then
        mkdir $base
    fi
    cp $file$base.data $base
    cp $file$base.var $base
done

