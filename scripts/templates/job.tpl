#!/bin/ksh 
#@ job_type = parallel
#  job_type = serial
#@ job_name = _NAME_

#@ total_tasks = _NBPROC_
#  node = _NBNODE_
#  tasks_per_node = _TASKPERNODE_
#  blocking = 1 
#@ multi_threads = _MULTITHREAD_
#  large_page = Y
#@ threads_per_task = _NBTHREAD_
#  requirements = ( machine == "mc12" )
#@ node_usage = not_shared

## Repertoire ou se trouvent les fichiers I/O, working_dir de votre JOB
#@ initialdir = _PATH_
#@ error  = _FILE_ERR_
#@ output = _FILE_OUT_
#@ notification = complete
#@ notify_user  = _USER_
#@ comment = pastix

#@ wall_clock_limit = _TIME_
#@ resources = mbytes(_MEM_)
#@ network.MPI = sn_all,not_shared,us

#@ environment = COPY_ALL; \
AIXTHREAD_MNRATIO=1:1 ; \
AIXTHREAD_SCOPE=S ; \
EXTSHM=ON ; \
MEMORY_AFFINITY=MCM@SHMEM=RR ; \
MP_SHARED_MEMORY=YES ; \
MP_INFOLEVEL=2 ; MP_MULTIPLE_THREAD=YES ; MP_PRINTENV=YES ;
# MEMORY_AFFINITY=MCM@SHMEM=RR ; \
# rset = rset_mcm_affinity
# mcm_affinity_option = mcm_mem_pref, mcm_distribute

## Si je veux lancer depuis REGATTA
#@ class = cluster_mc
#@ queue

#export MALLOCTYPE=watson
#export MALLOCTYPE=debug
#export XLSMPOPTS=stack=8589934592
#export MEMORY_AFFINITY=MCM
#export MP_SHARED_MEMORY=yes
#export MP_WAIT_MODE=poll
#export MP_EAGER_LIMIT=0
#export MP_EAGER_LIMIT=65536
#export MP_EUILIB=us
#export MP_BUFFER_MEM=64000000
#export MP_SINGLE_THREAD=no
#export AIXTHREAD_SCOPE=S
#export AIXTHREAD_RWLOCK_DEBUG=OFF
#export AIXTHREAD_MUTEX_DEBUG=OFF
#export AIXTHREAD_COND_DEBUG=OFF
#export OMP_NUM_THREADS=1

_EXECCMD_  _PATH_/_EXECUTABLE_ _DRIVER_
