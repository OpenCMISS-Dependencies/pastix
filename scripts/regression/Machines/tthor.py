import os
from remoteMachine import RemoteMachine

machines = {}
machinename = u'Ubuntu-VirtualBox-ramet'
address = u'localhost'
remote = {
    machinename : RemoteMachine(address,
                                u"ramet",
                                [u"ssh"],
                                [u"scp"],
                                u'/home/ramet/Work/pastix/results',
                                u'bash')
    }

script = u"""
#!/bin/bash

#source /opt/intel/Compiler/11.0/083/mkl/tools/environment/mklvarsem64t.sh
#source /opt/hwloc/env.sh

cd __DIR__
mpirun -np __PROCNBR__ __CMD__ > __LOG__ 2>&1 
"""
libmkl      = u"-lblas"
libgfortran = u"-lgfortran"

name = u"tthor"
machines[name] = remote[machinename].createConfig(name, 
                                                  machinename,
                                                  4, 
                                                  script,
                                                  scotch_int32 = u"/home/ramet/Work/scotch",
                                                  cc          = u"gcc", 
                                                  fc          = u"gfortran", 
                                                  mpicc       = u'mpicc', 
                                                  mpifc       = u'mpif90',
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
