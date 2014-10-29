#!/bin/bash

###OAR -p "cluster='borderline'"
#OAR -d _PATH_
#OAR -l /nodes=_NBNODE_,walltime=_TIME_
#OAR -n _NAME_
#OAR -E _FILE_ERR_
#OAR -O _FILE_OUT_
####OAR --notify "mail:_USER_"

i=0
machine=`uname -n | sed "s/\([a-zA-Z]*\)[-.][-0-9a-zA-Z.]*/\1/"`;

if [ $machine == "bordereau" ]
then
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/softs/acml/gnu64/lib/;
    export PATH=/home/bordeaux/mfaverge/bin:/home/bordeaux/mfaverge/bin:/usr/local/bin:/usr/bin:/bin:/usr/bin/X11:/usr/games;
elif [ $machine == "borderline" ]
then
    module load acml/gcc_int64/64/3.6.0;
    module unload mvapich2-tm/intel/64/1.0.1;
    source $HOME/env-mpichmad.sh
fi

# On récupère le nombre de noeuds présent dans le fichier
nbnode=`cat $OAR_NODEFILE | sort | uniq | wc -l | awk '{ print $1 }'`;

_EXECCMD_ mpirun _PROTOCOL_ -np $nbnode ./pastix _DRIVER_


