#PBS -N PaStiX
#PBS -o _FILE_OUT_
#PBS -e _FILE_ERR_
#PBS -j oe
#PBS -l walltime=_TIME_ 
#PBS -l nodes=_NBNODE_
## Positionne tous les noeud sur un même switch, pour plus de 16 le job
## restera bloque indefiniment
#PBS -W x=\"NODESET:ONEOF:FEATURE:ib001:ib002:ib005:ib006\"

source /etc/profile.d/modules.sh

cd _PATH_

# Chargement des modules :
module purge
#module load compiler/intel-11.1-069 mpi/openmpi-1.4.1
module load compiler/intel-11.1-069 mpi/mvapich2-1.4.1 mpiexec-0.83 

# Positionnement dans le bon repertoire qui va bien :
# cat $PBS_NODEFILE
cat $PBS_NODEFILE | sort | uniq > hostfile
#_EXECCMD_ -n _NBNODE_ -hostfile hostfile _EXECUTABLE_  _DRIVER_
_EXECCMD_ -n _NBNODE_  _EXECUTABLE_  _DRIVER_

