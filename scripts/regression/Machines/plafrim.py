# -*- coding: utf-8 -*-

import os
from remoteMachine import RemoteMachine

machines = {}


remote = {}
users=[u"ramet",u"lacoste"]
machinename=u"devel09.plafrim.cluster-%s"
for user in users:
    remote[machinename % user] = RemoteMachine(u"mygale.bordeaux.inria.fr",
                                               user,
                                               [u"ssh"],
                                               [u"scp"],
                                               # not in tmp as it is not shared 
                                               # between machines
                                               u'/lustre/%s/results' % user,
                                               u'qsub',
                                               qstat = u'qstat -f | grep __JOBID__ | wc -l',
                                               sleeptime = 5,
                                               modules=["compiler/gcc",
                                                        "mpi/openmpi",
                                                        "lib/mkl",
                                                        "hwloc",
                                                        "scotch/int32/5.1.12b"])


script = u"""
#PBS -N __JOBID__
#PBS -j oe
#PBS -l walltime=00:05:00 
#PBS -l nodes=__PROCNBR__:ppn=__THRDNBR__
### Positionne tous les noeud sur un meme switch, pour plus de 16 le job
### restera bloque indefiniment
### P_B_S -W x=\"NODESET:ONEOF:FEATURE:ib001:ib002:ib005:ib006\"


# Chargement des modules :
module purge
__MODULES__
# Positionnement dans le bon repertoire qui va bien :
# cat $PBS_NODEFILE
# cat $PBS_NODEFILE | sort | uniq > hostfile
cd __DIR__
module list  > __LOG__ 2>&1
echo 'mpirun -n __PROCNBR__ __CMD__' >> __LOG__ 2>&1
mpirun -npernode 1 -n __PROCNBR__ __CMD__ >> __LOG__ 2>&1

iter=0
while [ ! -s __LOG__  -a $iter -le 10 ]
  do
  sleep 1
  iter=`expr $iter + 1`
  done

# On attend encore 2 secondes de plus au cas ou.
sleep 2
if [ ! -s __LOG__ ]; then echo "fichier de log absent..."; fi
"""

name   = u"plafrim-gnu-openmpi-mkl-%s"
libmkl = u'-L${BLAS_DIR} ${BLAS_LIB}'
libgfortran = u"-lgfortran"

for user in users :
    machines[name % user] = remote[machinename % user].createConfig(name % user, 
                                                                    machinename % user, 
                                                                    32, 
                                                                    script,
                                                                    scotch_home = u"${SCOTCH_DIR}/..",
                                                                    hwloc_home  = u"$(shell pkg-config --libs-only-L hwloc|sed -e 's/-L//')",
                                                                    fcflags     = u'', 
                                                                    ccflags     = u'', 
                                                                    ccopt       = u'-O3', 
                                                                    ccdeb       = u'-g3',
                                                                    libblas     = libmkl, 
                                                                    libfortran  = libgfortran,
                                                                    libmaths    = u'-lm', 
                                                                    librt       = u'-lrt',
                                                                    ar          = u'ar', 
                                                                    arflags     = u'-ruv', 
                                                                    make        = u'make', 
                                                                    make_j      = u'make -j 8')
