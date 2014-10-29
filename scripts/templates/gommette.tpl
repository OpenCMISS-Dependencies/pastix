#!/bin/bash

###OAR -p "cluster='borderline'"
#OAR -d _PATH_
#OAR -l /nodes=_NBNODE_,walltime=_TIME_
#OAR -n _NAME_
#OAR -E _FILE_ERR_
#OAR -O _FILE_OUT_
####OAR --notify "mail:_USER_"

source ~/opt/scotch/int64/env.sh;
#module load scotch/current/int64;
module load acml/gcc-int64/64/3.6.0;
export MPI_NMAD_STRATEGY=_NMADSTRAT_;
export MPI_NMAD_PROTOCOL=_PROTOCOL_;
export MV2_RNDV_PROTOCOL=RGET;

function LoadMPD
{

    # Le fichier $OAR_NODE_FILE contient une ligne par core, donc on reduit Å‡ une ligne par noeud
    echo "-- Creation du fichier host --"
    cat $OAR_NODE_FILE | uniq > /tmp/hosts.$OAR_JOBID;
    
    # On rÅÈcupÅËre le nombre de noeuds prÅÈsent dans le fichier
    nbnode=`wc -l /tmp/hosts.$OAR_JOBID | awk '{ print $1 }'`;
    
    # On passe le noeud courant en derniÅËre ligne (on cherche pas Å‡ savoir pq :-))
    sed "1h;1d;${nbnode}G" /tmp/hosts.$OAR_JOBID > /tmp/hosts2.$OAR_JOBID
    
    # for i in `cat /tmp/hosts2.$OAR_JOBID`
    # do
    #     echo $i;
    #     oarsh $i `(mkdir -p /tmp/mfaverge/$OAR_JOBID && cp -r _PATH_ /tmp/mfaverge/$OAR_JOBID)`
    # done

    echo "-- Lancement de mpd --"
    
    # On boucle sur le lancement du dÅÈmon mpd en cas de pb
    i=0;
    ret=1;
    while [ $ret -ne 0 ]
      do
      i=$(($i + 1))
      echo "Tentative $i\n"
      `sleep 2`
      mpdcleanup -f /tmp/hosts2.$OAR_JOBID
      mpdboot -n $nbnode -f /tmp/hosts2.$OAR_JOBID -r oarsh
      ret=$?;
    done
}

export PM2_FLAVOR=_PM2_FLAVOR_
if [ -z $PM2_FLAVOR ] || [ $PM2_FLAVOR == "pastix-wom" ]
then 
    module load mvapich2-tm/gcc/64;
    LoadMPD;
else
    source $HOME/env/env-mpichmad.sh
fi

_EXECCMD_ mpirun -np _NBPROC_ _EXECUTABLE_ _DRIVER_

# for i in `cat /tmp/hosts2.$OAR_JOBID`
# do
#     oarsh $i cp /tmp/mfaverge/$OAR_JOBID/trace* _PATH_
# done
mpdallexit