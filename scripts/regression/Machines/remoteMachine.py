from copy import deepcopy
import subprocess
import time
import os

class RemoteMachine():
    """ Defines procedures to do remote actions on remote machines """
    def __init__(self, address, user, ssh, scp, 
                 resultdir, submit, qstat = None, sleeptime = 30,
                 mpdboot = None, mpdclean = None, modules=None):
        self.address    = address
        self.username   = user
        self.ssh        = ssh
        self.scp        = scp
        self.resultdir  = resultdir
        self.sleeptime  = sleeptime
        self.submit     = submit
        self.qstat      = qstat
        self.mpdboot    = mpdboot
        self.mpdclean   = mpdclean
        self.modules    = modules
        self.module_loads = u""
        if not modules == None:
            for module in modules:
                self.module_loads += u"module load %s\n" % module

        
        
    def remoteCallCommand(self, cmd):
        """
        Return the subprocess calling the given command on the given
        machine
        """ 
        mycmd = deepcopy(self.ssh)
        mycmd.append(self.username+u'@'+self.address)
        mycmd.append(u'-o')
        mycmd.append(u'ForwardAgent=yes')
        mycmd.append(" ".join(cmd))
        if __debug__:
            print u'Running "%s" on %s' %(str(cmd), self.address)
        return subprocess.Popen(mycmd, 
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)

    def fileCopyFrom(self, source, destination):
        """
        Return the subprocess copying the given file from the given
        machine onto the destination
        """ 
        mycmd = deepcopy(self.scp)
        mycmd.append(self.username+u'@'+self.address+u':'+source)
        mycmd.append(destination)
        return  subprocess.Popen(mycmd,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)


    def fileCopyTo(self, source, destination):
        """
        Return the subprocess copying the given file from the given
        machine onto the destination
        """ 
        mycmd = deepcopy(self.scp)
        mycmd.append(source)
        mycmd.append(self.username+u'@'+self.address+u':'+destination)
        return  subprocess.Popen(mycmd,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)

    def createConfig(self, name, machinename, nproc, script, 
                     arch=u"i686_pc_linux", 
                     scotch_home=None,
                     scotch_int32=None, scotch_int64=None,
                     scotch_int=None, scotch_long=None,
                     metis_home=None,
                     metis_int32=None, metis_int64=None,
                     metis_int=None, metis_long=None,
                     hwloc_home=None,
                     cc = u"gcc", fc =u"gfortran",
                     mpicc= u'mpicc', mpifc = u'mpif90',
                     fcflags = u'', ccflags = u'',
                     ccopt=u'-O3', ccdeb=u'-g3',
                     libblas = u'-lblas', libfortran=u'-lgfortran',
                     libmaths=u'-lm', librt=u'-lrt',
                     ar = u'ar', arflags = u'-ruv',
                     make= u'make', make_j=u'make -j 8'):
        config = {}
        config[u"NAME"        ] = name
        config[u"machinename" ] = machinename
        config[u"username"    ] = self.username
        config[u"nproc"       ] = nproc
        config[u"ROOT"        ] = None # Will be set later.
        config[u"ARCH"        ] = arch
        defined = False
        if not scotch_home == None:
            defined = True
            config[u"SCOTCH_INT32"] = os.path.join(scotch_home, u"int32")
            config[u"SCOTCH_INT64"] = os.path.join(scotch_home, u"int64")
            config[u"SCOTCH_LONG" ] = os.path.join(scotch_home, u"long")
            config[u"SCOTCH_INT"  ] = os.path.join(scotch_home, u"int")
        if not scotch_int32 == None:
            defined = True
            config[u"SCOTCH_INT32"] = scotch_int32
        if not scotch_int64 == None:
            defined = True
            config[u"SCOTCH_INT64"] = scotch_int64
        if not scotch_int == None:
            defined = True
            config[u"SCOTCH_INT"]   = scotch_int
        if not scotch_long == None:
            defined = True
            config[u"SCOTCH_LONG"] = scotch_long
        if defined == False:
            print "warning: no scotch defined for machine %s" % name

        if not metis_home == None:
            defined = True
            config[u"METIS_INT32"] = os.path.join(metis_home, u"int32")
            config[u"METIS_INT64"] = os.path.join(metis_home, u"int64")
            config[u"METIS_LONG" ] = os.path.join(metis_home, u"long")
            config[u"METIS_INT"  ] = os.path.join(metis_home, u"int")
        if not metis_int32 == None:
            defined = True
            config[u"METIS_INT32"] = metis_int32
        if not metis_int64 == None:
            defined = True
            config[u"METIS_INT64"] = metis_int64
        if not metis_int == None:
            defined = True
            config[u"METIS_INT"]   = metis_int
        if not metis_long == None:
            defined = True
            config[u"METIS_LONG"] = metis_long
        if defined == False:
            print "warning: no metis defined for machine %s" % name

        config[u"CC"          ] =  cc
        config[u"FC"          ] =  fc
        config[u"MPCC"        ] =  mpicc
        config[u"MPFC"        ] =  mpifc
        config[u"FCFLAGS"     ] =  fcflags
        config[u"CCFLAGS"     ] =  ccflags
        config[u"COPT"        ] =  ccopt
        config[u"CDEB"        ] =  ccdeb
        config[u"LIBBLAS"     ] =  libblas
        config[u"LIBFORTRAN"  ] =  libfortran
        config[u"LIBMATH"     ] =  libmaths
        config[u"LIBREALTIME" ] =  librt
        config[u"AR"          ] =  ar
        config[u"ARFLAGS"     ] =  arflags
        config[u"MAKE"        ] =  make
        config[u"MAKE_J"      ] =  make_j
        config[u"result_dir"  ] =  self.resultdir
        config[u"script"      ] =  script.replace(u"__MODULES__",
                                                  self.module_loads)
        config[u"soumission"  ] =  self.submit
        config[u"MODULES"     ] =  self.modules
        config[u"MODULE_LOADS"] =  self.module_loads
        if not hwloc_home == None:
            config[u'HWLOC_HOME'  ] =  hwloc_home
        
        return config

    def waitJobs(self, ident):
        flag = True
        if self.qstat == None: flag = False
        while (flag):
            qstat = self.qstat.replace(u'__JOBID__', ident)
            process = subprocess.Popen(self.qstat, shell=True, 
                                       stdout=subprocess.PIPE, 
                                       stderr=subprocess.STDOUT)
            output  = process.communicate()[0]
            if (int(output) == 0): break
            time.sleep(self.sleeptime)
        # Because some logs seems not to be synced
        time.sleep(60)
