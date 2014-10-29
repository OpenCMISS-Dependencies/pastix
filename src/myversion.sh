#!/bin/sh

if [ -e ".git" ]; then
    rev=`git rev-parse --short HEAD`;
    echo $rev;
    exit;
else
    if [ -x "`which svnversion`" ]; then
	rev=`LC_ALL=C svnversion 2>/dev/null | sed -e 's/ //g' -e 's/,//g'`
	if [ "$rev" != 'exporté' ]
	then
	    echo $rev;
	    exit;
	fi;
    fi
fi

if [ -f Revision ]
then
    cat ./Revision
else
    if [ -f ../../Revision ]
    then
	cat ../../Revision
    else
	echo "NoRevisionFound"
    fi
fi;
