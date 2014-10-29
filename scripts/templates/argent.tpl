#!/bin/bash

#BSUB -J _NAME_
#BSUB -n _NBPROCTOTAL_
#BSUB -W _TIME_
#BSUB -e _FILE_ERR_
#BSUB -o _FILE_OUT_
#BSUB -u _USER_
#BSUB -N

set -vx

SRCPATH=_RESPATH_/_RESPATH2_;
DSTPATH=${SCRATCHDIR}/_RESPATH2_;

#****************************************************
#  CREATION OF THE WORK DIRECTORY
#****************************************************
# create the work directory (in $SCRATCHDIR)
mkdir -p ${DSTPATH};

# recovery of the input data file
cp ${SRCPATH}/iparm.txt ${DSTPATH}/
cp ${SRCPATH}/dparm.txt ${DSTPATH}/
cp ${SRCPATH}/pastix    ${DSTPATH}/
ln -sf _MATRICE_ ${DSTPATH}/_LINK_

#****************************************************
#  SCRIPT SUBMISSION
#****************************************************
cd ${DSTPATH}
# export MV2_SRQ_SIZE=1000
export VIADEV_CLUSTER_SIZE=LARGE
srun -n _NBPROC_ -c _NBTHREAD_ -t _TIME_ -K ./pastix _DRIVER_
