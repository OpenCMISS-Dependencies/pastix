#PBS -S /bin/bash
#PBS -N Pastix__NBNODE_P__NBTHREAD_T
#PBS -e _FILE_ERR_
#PBS -o _FILE_OUT_
#PBS -l select=_NBNODE_:ncpus=_NBTHREAD_:mpiprocs=_TASKPERNODE_
#PBS -l walltime=_TIME_

cd _PATH_

## Redefinition des modules
. /etc/profile.d/modules.sh
## Initialisation environnement des modules
module purge
## Chargement du module intel (compilateur)
module load intel
## Chargement du module intelmpi (librairie // MPI de intel)
module load intelmpi
## recuperation de la liste des nœuds de calcul utilises
## et creation du fichier mpd.hosts necessaire a intelmpi
cat $PBS_NODEFILE | uniq > mpd.hosts
## Demarrage du daemon mpd de intelmpi
mpdboot --rsh=ssh -v -n `cat mpd.hosts|wc -l` -f mpd.hosts
## lancement executable
mpiexec -np _NBNODE_ ./_EXECUTABLE_
##Arret des daemons mpd
mpdallexit
