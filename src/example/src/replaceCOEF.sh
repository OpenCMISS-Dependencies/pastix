#/bin/sh
source=$1
destination=$2
shift 2
echo $* | grep DTYPE_COMPLEX >/dev/null
complex=$?
echo $* | grep DDISTRIBUTED >/dev/null
distributed=$?

if [ $complex -eq 0 ]
    then
    if  [ $distributed -eq 0 ]
	then
	sed -e 's/__IS_COMPLEX__/.true./g' -e 's/__COEF__/COMPLEX(KIND=MURGE_COEF_KIND)/g' -e 's/ifdef\ DISTRIBUTED/if\ 1/g' $source > $destination;
    else
	sed -e 's/__IS_COMPLEX__/.true./g' -e 's/__COEF__/COMPLEX(KIND=MURGE_COEF_KIND)/g' -e 's/ifdef\ DISTRIBUTED/if\ 0/g' $source > $destination;
    fi
else
    if  [ $distributed -eq 0 ]
	then
	sed -e 's/__IS_COMPLEX__/.false./g' -e 's/__COEF__/REAL(KIND=MURGE_COEF_KIND)/g' -e 's/ifdef\ DISTRIBUTED/if\ 1/g' $source > $destination;
    else
	sed -e 's/__IS_COMPLEX__/.false./g' -e 's/__COEF__/REAL(KIND=MURGE_COEF_KIND)/g' -e 's/ifdef\ DISTRIBUTED/if\ 0/g' $source > $destination;
    fi
fi

