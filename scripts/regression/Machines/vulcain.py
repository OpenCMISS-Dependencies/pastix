import os
from remoteMachine import RemoteMachine

machines = {}
users = [u"pastix-user", u"lacoste"]

machinename = u'vulcain.bordeaux.inria.fr-%s'
address = u'vulcain.bordeaux.inria.fr'
remote = {}
for user in users:
    remote[machinename % user] = RemoteMachine(address,
                                               user,
                                               [u"ssh"],
                                               [u"scp"],
                                               # not in tmp as it is not shared 
                                               # between machines
                                               u'/media/vulcain/%s/results-python' % user,
                                               u'bash')
    
script = u"""
#!/bin/bash

source /opt/intel/composer_xe_2013/mkl/bin/mklvars.sh intel64
source /opt/hwloc/env.sh
source /opt/mpich2/env.sh

cleanup()
{
    trap - ALRM               #reset handler to default
    kill -ALRM $a 2>/dev/null #stop timer subshell if running
    kill $! 2>/dev/null &&    #kill last job
      exit 124                #exit with 124 if it was running
}

watchit()
{
    trap "cleanup" ALRM
    sleep $1& wait
    kill -ALRM $$
}


cd __DIR__
echo vulcain > hostfile
echo hephaistos >> hostfile
watchit 10 & a=$!
mpirun -f hostfile -np __PROCNBR__ timeout 10 __CMD__ > __LOG__ 2>&1 & wait $!;RET=$?
kill -ALRM $a
wait $a
exit $RET
"""
libmkl      = u"-L/opt/intel/composer_xe_2013/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
libgfortran = u"-lgfortran"

name = u"vulcain-gnu-mpich2-mkl-%s"

for user in users :
    machines[name % user] = remote[machinename % user].createConfig(name % user, 
                                                                    machinename % user,
                                                                    32, 
                                                                    script,
                                                                    scotch_home = u"/opt/scotch/",
                                                                    metis_int   = u"/opt/metis/",
                                                                    metis_int32 = u"/opt/metis/",
                                                                    hwloc_home  = u"/opt/hwloc",
                                                                    cc          = u"gcc", 
                                                                    fc          = u"gfortran", 
                                                                    mpicc       = u'/opt/mpich2/bin/mpicc', 
                                                                    mpifc       = u'/opt/mpich2/bin/mpif90',
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

source /opt/intel/composer_xe_2013/bin/compilervars.sh intel64
source /opt/intel/composer_xe_2013/mkl/bin/mklvars.sh intel64
source /opt/hwloc/env.sh

cd __DIR__
/opt/mpich2/bin/mpirun -np __PROCNBR__ __CMD__ > __LOG__ 2>&1 
"""
 
name = u'vulcain-intel-mpich2-mkl-%s'
libifcore = u'-L/opt/intel/Compiler/11.0/083/lib/intel64/ -lifcore'
for user in users :
    machines[name % user] = remote[machinename % user].createConfig(name, 
                                                                    machinename,
                                                                    32, 
                                                                    script,
                                                                    scotch_home = u"/opt/scotch/",
                                                                    metis_int   = u"/opt/metis/int",
                                                                    metis_int32 = u"/opt/metis/int32",
                                                                    hwloc_home  = u"/opt/hwloc",
                                                                    cc          = u"/opt/intel/Compiler/11.0/083/bin/intel64/icc", 
                                                                    fc          = u"/opt/intel/Compiler/11.0/083/bin/intel64/ifort", 
                                                                    mpicc       = u'/opt/intel/mpich2/bin/mpicc -cc=/opt/intel/Compiler/11.0/083/bin/intel64/icc', 
                                                                    mpifc       = u'/opt/intel/mpich2/bin/mpif90 -f90=/opt/intel/Compiler/11.0/083/bin/intel64/ifort',
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

