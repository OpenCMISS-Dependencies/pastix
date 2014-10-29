import os
from remoteMachine import RemoteMachine

machines = {}
users = [u"mfaverge"]

machinename = u'remus.eecs.utk.edu-%s'
address = u'remus.eecs.utk.edu'
remote = {}
for user in users:
    remote[machinename % user] = RemoteMachine(address,
                                               user,
                                               [u"ssh"],
                                               [u"scp"],
                                               # not in tmp as it is not shared 
                                               # between machines
                                               u'/home/%s/results-python' % user,
                                               u'bash')
    
script = u"""
#!/bin/bash

source /home/mfaverge/trace/opt/trace_env.sh
source /home/mfaverge/mklvars.sh intel64
export LD_LIBRARY_PATH=/home/mfaverge/opt/x86_64/lib:$LD_LIBRARY_PATH
export PATH=/home/mfaverge/opt/x86_64/bin:$PATH
export PKG_CONFIG_PATH=/home/mfaverge/opt/x86_64/lib/pkgconfig/:${PKG_CONFIG_PATH}

cd __DIR__
echo remus > hostfile
#echo romulus >> hostfile
mpirun -f hostfile -np __PROCNBR__ __CMD__ > __LOG__ 2>&1
"""

libmkl      = u"-L/mnt/scratch/sw/intel/11.1.069/mkl/lib/em64t/ -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
libgfortran = u"-lgfortran"

name = u"remus-gnu-mpich2-mkl-%s"

for user in users :
    machines[name % user] = remote[machinename % user].createConfig(name % user, 
                                                                    machinename % user,
                                                                    48, 
                                                                    script,
                                                                    scotch_int64 = u"/home/mfaverge/INRIA/scotch_5.1.11/",
                                                                    #metis_int   = u"",
                                                                    #metis_int32 = u"",
                                                                    hwloc_home  = u"/mnt/scratch/plasma/sw/hwloc-1.0.1/",
                                                                    cc          = u"gcc", 
                                                                    fc          = u"gfortran", 
                                                                    mpicc       = u'/home/mfaverge/opt/x86_64/bin/mpicc',
                                                                    mpifc       = u'/home/mfaverge/opt/x86_64/bin/mpif90',
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


script = u"""
#!/bin/bash

source /mnt/scratch/sw/intel/11.1.069/mkl/tools/environment/mklvarsem64t.sh 
source /mnt/scratch/sw/intel/11.1.069/bin/iccvars.sh intel64
source /mnt/scratch/sw/intel/11.1.069/bin/ifortvars.sh intel64

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/scratch/plasma/sw/hwloc-1.0.1/lib
export PATH=$PATH:/mnt/scratch/plasma/sw/hwloc-1.0.1/bin
export PKG_CONFIG_PATH=${PKG_CONFIG_PATH}:/mnt/scratch/plasma/sw/hwloc-1.0.1/lib/pkgconfig/


cd __DIR__
/opt/mpich2/bin/mpirun -np __PROCNBR__ __CMD__ > __LOG__ 2>&1 
"""
 
name = u'remus-intel-mpich2-mkl-%s'
libifcore = u'-L/opt/intel/Compiler/11.0.069/lib/intel64/ -lifcore'
cc = u"/mnt/scratch/sw/intel/11.1.069/bin/intel64/icc", 
fc = u"/mnt/scratch/sw/intel/11.1.069/bin/intel64/ifort", 
for user in users :
    machines[name % user] = remote[machinename % user].createConfig(name, 
                                                                    machinename,
                                                                    48, 
                                                                    script,
                                                                    scotch_home = u"/home/mfaverge/INRIA/scotch_5.1.11/",
                                                                    #metis_int   = u"/opt/metis/int",
                                                                    #metis_int32 = u"/opt/metis/int32",
                                                                    hwloc_home  = u"/mnt/scratch/plasma/sw/hwloc-1.0.1",
                                                                    cc          = cc, 
                                                                    fc          = fc, 
                                                                    mpicc       = u'/home/mfaverge/opt/x86_64/bin/mpicc -cc=%s'   % (cc), 
                                                                    mpifc       = u'/home/mfaverge/opt/x86_64/bin/mpif90 -fc=%s' % (fc),
                                                                    fcflags     = u'', 
                                                                    ccflags     = u'', 
                                                                    ccopt       = u'-O3', 
                                                                    ccdeb       = u'-g3',
                                                                    libblas     = libmkl, 
                                                                    libfortran  = libifcore,
                                                                    libmaths    = u'-lm', 
                                                                    librt       = u'-lrt',
                                                                    ar          = u'ar', 
                                                                    arflags     = u'-ruv', 
                                                                    make        = u'make', 
                                                                    make_j      = u'make -j 8')

