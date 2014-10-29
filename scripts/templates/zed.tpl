#!/bin/bash

###OAR -p "cluster='chinqchint'"
#OAR -d _PATH_
#OAR -l /nodes=_NBNODE_,walltime=_TIME_
#OAR -n _NAME_
#OAR -E _FILE_ERR_
#OAR -O _FILE_OUT_
####OAR --notify "mail:_USER_"

source $HOME/env/env-scotch.sh
source $HOME/env/env-gotoblas.sh
source $HOME/env/env-mpich2.sh

export PM2_FLAVOR=_PM2_FLAVOR_

# Le fichier $OAR_NODE_FILE contient une ligne par core, donc on reduit à une ligne par noeud
echo "-- Creation du fichier host --"
cat $OAR_NODE_FILE | sort | uniq | sed 's/\..*//' > /tmp/hosts.$OAR_JOBID;

# On récupère le nombre de noeuds présent dans le fichier
nbnode=`wc -l /tmp/hosts.$OAR_JOBID | awk '{ print $1 }'`;

# On passe le noeud courant en dernière ligne (on cherche pas à savoir pq :-))
sed "1h;1d;${nbnode}G" /tmp/hosts.$OAR_JOBID > /tmp/hosts2.$OAR_JOBID

echo "-- Lancement de mpd --"
    
# On boucle sur le lancement du démon mpd en cas de pb
i=0;
ret=1;
while [ $ret -ne 0 ]
  do
  i=$(($i + 1))
  echo "Tentative $i\n"
  `sleep 2`;
  mpdboot -n $nbnode -f /tmp/hosts2.$OAR_JOBID -r oarsh
  ret=$?;
done
    

_EXECCMD_ mpirun -np _NBPROC_ ./pastix _DRIVER_ 

