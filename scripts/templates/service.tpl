#PBS -S /bin/bash
#PBS -N Pastix__NBNODE_P__NBTHREAD_T
#PBS -o _FILE_OUT_
#PBS -e _FILE_ERR_
#PBS -j oe
#PBS -m ae -M _USER_
#PBS -l select=_NBNODE_:ncpus=_NBTHREAD_:mpiprocs=_TASKPERNODE_
#PBS -l walltime=_TIME_

cd _PATH_
/usr/pbs/bin/mpiexec ./_EXECUTABLE_

