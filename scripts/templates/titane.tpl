#!/bin/bash

#MSUB -r _NAME_
#MSUB -n _NBPROCTOTAL_
#BSUB -W 10
#MSUB -T _TIME_
#MSUB -e _FILE_ERR_
#MSUB -o _FILE_OUT_
#BSUB -u _USER_
#BSUB -N

export DAT_OVERRIDE=$HOME/dat.conf 

SRCPATH=_PATH1_/_PATH2_;

#****************************************************
#  SCRIPT SUBMISSION
#****************************************************

cd ${SRCPATH}
# export MV2_SRQ_SIZE=1000

export VIADEV_CLUSTER_SIZE=LARGE
source /applications/intel/impi/3.2.1.009/bin64/mpivars.sh
export I_MPI_WAIT_MODE=enable
export I_MPI_SPIN_COUNT=100000
/applications/intel/impi/3.2.1.009/bin64/mpirun -r ssh -perhost 1 -np _NBPROC_ ./_EXECUTABLE_ _DRIVER_
