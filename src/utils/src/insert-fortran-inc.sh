#!/bin/sh

usage="usage: $0 file1 file2\n";
usage="$usage\t\t Include file2 into file1 just before END INTERFACE.\n";
if [ $# = 2 ]
then
    LINE=$(awk '/END INTERFACE/ {print NR}' $1);
    LENGHT=$(wc -l $1 | awk '{print $1}');

    head -n $(($LINE - 1)) $1 >> tmp;
    cat $2 >> tmp;
    tail -n $(($LENGHT - $LINE + 1)) $1 >> tmp;
    mv tmp $1;
else
    echo -e $usage;
fi;