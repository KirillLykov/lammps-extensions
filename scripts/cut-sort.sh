#! /bin/bash
if [[ $# == 1 ]] && [[ $1 == "--help" ]] ; then 
echo "Used for sorting part of atoms files. Input parameters are filename, from which line to which you need to sort"
exit
fi

if [ $# -ne 3 ]; then
echo "Usage: cut-sort.sh <file_name> <start line> <end line>"
exit
fi

until=$(($2 - 1))
after=$(($3 + 1))
head -n $until $1 > OUT
sed -n "$2,$3p" $1 | sort -k1 >> OUT
tail -n +$after $1 >> OUT

