#!/bin/csh
# @ class = PC
# @ initialdir = _PATH_
# @ input = /dev/null
# @ output = _FILE_OUT_
# @ error  = _FILE_ERR_
# @ job_type=parallel
# @ job_name = _NAME_
# @ node = _NBNODE_
# @ tasks_per_node = _TASKPERNODE_
# @ notify_user = _USER_
# @ notification = complete
# @ network.MPI=csss,shared,US,AVERAGE
# @ resources = ConsumableCpus(1) ConsumableMemory(900mb)
# @ wall_clock_limit = _TIME_
#  requirements = (Feature == "F17")
# @ node_usage = not_shared
# @ environment = COPY_ALL;\
# LANG=En_US;\
# MP_SHARED_MEMORY=YES;\
# MP_LABELIO=yes;\
# MP_EUILIB=us;\
# MP_HOSTFILE=NULL;\
# MP_RMPOOL=0;\
# OMP_NUM_THREADS=1;\
# XLSMPOPTS=parthds=1;\
# MP_WAIT_MODE=poll;\
# MP_CSS_INTERRUPT=yes;
# @ restart = no
# @ queue
_PATH_/pastix

