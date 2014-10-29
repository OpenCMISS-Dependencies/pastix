#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement # This isn't required in Python 2.6
import os
import subprocess
import sys
import hashlib
import re
from project import Project

import xml.etree.cElementTree as xml

from copy import deepcopy


class PaStiX(Project):
    """Class defining PaStiX project."""
    _name = u'PaStiX'
    ## Compilation option
    _liste_default_yes = [u"NUMA_ALLOC",
                          u"MEMORY_USAGE"]

    _liste_default_not = [u"FORCE_CONSO",
                          u"NOSMP_RAFF",
                          u"STATS_SOPALIN",
                          u"ONLY_LOAD_SYMBMTX",
                          u"COMPACT_SMX",
                          u"DOUBLE_DOWN_TIME",
                          u"FORCE_DIAGDOM",
                          u"TEST_IRECV",
                          u"NO_MPI_TYPES",
                          u"THREAD_COMM",
                          u"COMM_REORDER",
                          u"PASTIX_FUNNELED",
                          u"MULT_SMX",
                          u"PASTIX_ESC",
                          u"PASTIX_DYNSCHED",
                          u"TRACE_SOPALIN",
                          u"TEST_REPART"]
    _optionsComp = {}
    _binaries = {}


        ## CONFIGURATION FILE

    _config  = u"""
HOSTARCH    = __HOSTARCH__
VERSIONBIT  = __VERSIONBIT__
EXEEXT      =
OBJEXT      = .o
LIBEXT      = .a
CCPROG      = __CC__ __CFLAGS__
CFPROG      = __FC__ __FFLAGS__
CF90PROG    = gfortran -ffree-form
MCFPROG     = __MPFC__  __FFLAGS__
CF90CCPOPT  = -ffree-form -x f95-cpp-input
# Compilation options for optimization (make expor)
CCFOPT      = __COPTFLAGS__

LKFOPT      =
MKPROG      = __MAKE__
MPCCPROG    = __MPCC__ __CFLAGS__
CPP         = cpp
ARFLAGS     = __ARFLAGS__
ARPROG      = __AR__
EXTRALIB    = __LIBFORT__ __LIBMATH__ __LIBREALTIME__

VERSIONMPI  = __VERSIONMPI__
VERSIONSMP  = __VERSIONSMP__
VERSIONSCH  = _sc
VERSIONINT  = __VERSIONINT__
VERSIONPRC  = __VERSIONPRC__
VERSIONFLT  = __VERSIONFLT__
VERSIONORD  = __VERSIONORD__

###################################################################
#                          INTEGER TYPE                           #
###################################################################
# uncomment the following lines for integer type support (Only 1)

CCTYPES     = __INTSIZE__

###################################################################
#                           FLOAT TYPE                            #
###################################################################
CCTYPESFLT  = __PRECISION__ __COMPLEX__


###################################################################
#                          MPI/THREADS                            #
###################################################################

# uncomment the following lines for sequential (NOMPI) version
CCTYPES    := $(CCTYPES) __MPI__

# uncomment the following lines for non-threaded (NOSMP) version
CCTYPES    := $(CCTYPES) __SMP__

# Uncomment the following line to enable a progression thread
CCPASTIX   := $(CCPASTIX) __OPT__

# Uncomment the following lines for Out-of-core
CCPASTIX   := $(CCPASTIX) __OOC__

###################################################################
#                      GRAPH PARTITIONING                         #
###################################################################

# uncomment the following lines for using metis ordering
#VERSIONORD  = _metis
METIS_HOME  = __METIS_DIR__

# Scotch always needed to compile
SCOTCH_HOME = __SCOTCH_DIR__

CCPASTIX   := $(CCPASTIX) __PARTITIONNER__
EXTRALIB   := $(EXTRALIB) __LIBPART__

###################################################################
#                Portable Hardware Locality                       #
###################################################################
# If HwLoc library is available, uncomment the following lines to bind correctly threads on cpus
HWLOC_HOME ?= __HWLOC_HOME__
HWLOC_INC  ?= $(HWLOC_HOME)/include
HWLOC_LIB  ?= $(HWLOC_HOME)/lib
CCPASTIX   := $(CCPASTIX) -I$(HWLOC_INC) -DWITH_HWLOC
EXTRALIB   := $(EXTRALIB) -L$(HWLOC_LIB) -lhwloc

###################################################################
#                             MARCEL                              #
###################################################################

# Uncomment following lines for marcel thread support
#VERSIONSMP := $(VERSIONSMP)_marcel
#CCPASTIX   := $(CCPASTIX) `pm2-config --cflags` -I${PM2_ROOT}/marcel/include/pthread
#EXTRALIB   := $(EXTRALIB) `pm2-config --libs`
# ---- Thread Posix ------
EXTRALIB   := $(EXTRALIB) -lpthread

# Uncomment following line for bubblesched framework support (need marcel support)
#VERSIONSCH  = _dyn
#CCPASTIX   := $(CCPASTIX) -DPASTIX_BUBBLESCHED

###################################################################
#                              BLAS                               #
###################################################################

# Choose Blas library (Only 1)
# Do not forget to set BLAS_HOME if it is not in your environnement
# BLAS_HOME=/path/to/blas
#----  Blas    ----
BLASLIB  = __LBLAS__

###################################################################
#                         PYTHON WRAPPER                          #
###################################################################
#MPI4PY_DIR    =
#SWIG_INC      = $(MPI4PY_DIR)/src/include/mpi4py/
#MPI4PY_LIBDIR = $(MPI4PY_DIR)/build/lib.linux-x86_64-2.6/
#CCTYPES       = -fPIC

###################################################################
#                          DO NOT TOUCH                           #
###################################################################

FOPT      := $(CCFOPT)
FDEB      := $(CCFDEB)
CCHEAD    := $(CCPROG) $(CCTYPES) $(CCFOPT)
CCFOPT    := $(CCFOPT) $(CCTYPES) $(CCPASTIX)
CCFDEB    := $(CCFDEB) $(CCTYPES) $(CCPASTIX)


###################################################################
#                        MURGE COMPATIBILITY                      #
###################################################################

MAKE     = $(MKPROG)
CC       = $(MPCCPROG)
CFLAGS   = $(CCFOPT) $(CCTYPESFLT)
FC       = $(MCFPROG)
FFLAGS   = $(CCFOPT)
LDFLAGS  = $(EXTRALIB) $(BLASLIB)
"""
    _config_name = u"config.in"
    _optionsComp = {
        u'INTEGER' : {
            u'database_id' : u'INTEGER',
            u'default' : {
                u'replace'      : { u'__VERSIONINT__': u'int_' ,
                                    u'__INTSIZE__'   : u''},
                u'database_val' : u'int'
                },
            u'int32'   : {
                u'replace'      : { u'__VERSIONINT__': u'int32_' ,
                                    u'__INTSIZE__'   : u'-DINTSIZE32' },
                u'database_val' : u'int32',
                u'searchInName' : re.compile(u'_int32')
                },
            u'int64'   : {
                u'replace'      : { u'__VERSIONINT__': u'int64_',
                                    u'__INTSIZE__'   : u'-DINTSIZE64'  },
                u'database_val' : u'int64',
                u'searchInName' : re.compile(u'_int64')
                },
            u'long'    : {
                u'replace'      : { u'__VERSIONINT__': u'long_' ,
                                    u'__INTSIZE__'   : u'-DLONG' },
                u'database_val' : u'long',
                u'searchInName' : re.compile(u'_long')
                },
            u'int'     : {
                u'replace'      : { u'__VERSIONINT__': u'int_',
                                    u'__INTSIZE__'   : u''  },
                u'database_val' : u'int'
                }
            },
        
        u'NOMPI' : {
            u'database_id' : u'MPI',
            u'default' : {
                u'define'      : u'',
                u'replace'     : {u'__VERSIONMPI__'    : u'mpi_',
                                  u'__MPI__'           : u''},
                u'database_val': u'MPI'
                },
            True       : {
                u'replace'      : {
                    u'__VERSIONMPI__'    : u'nompi_',
                    u'__MPI__'           : u'-DFORCE_NOMPI\n'+
                                           u'MPCCPROG    = $(CCPROG)\n'+
                                           u'MCFPROG     = $(CFPROG)'},
                u'database_val' : u'NOMPI',
                u'searchInName' : re.compile(u'_nompi')
                }
            },

        u'NOSMP' : {
            u'database_id' : u'SMP',
            u'default'     : {
                u'replace'      : {u'__VERSIONSMP__'    : u'smp_',
                                   u'__SMP__'           : u''},
                u'database_val' : u'SMP'
                },
            True           : {
                u'replace'      : {u'__VERSIONSMP__'    : u'nosmp_',
                                   u'__SMP__'           : u'-DFORCE_NOSMP'},
                u'database_val' : u'NOSMP',
                u'searchInName' : re.compile(u'_nosmp')
                },
            },

        u'PRECISION' : {
            u'database_id' : u'PRECISION',
            u'default'     : {
                u'replace'      : { u'__VERSIONPRC__' : u'simple_',
                                    u'__PRECISION__'  : u''},
                u'database_val' : u'simple'},
            u'double'      : {
                u'define'       : u'',
                u'replace'      : { u'__VERSIONPRC__' : u'double_',
                                    u'__PRECISION__'  : u'-DPREC_DOUBLE'},
                u'database_val' : u'double',
                u'searchInName' : re.compile(u'_double')}
            },

        u'COMPLEX' : {
            u'database_id'  : u'COMPLEX',
            u'default'      : {
                u'define'       : u'',
                u'replace'      : { u'__VERSIONFLT__' : u"real_",
                                    u'__COMPLEX__'    : u''},
                u'database_val' : u'real'},
            True            : {
                u'define'       : u'',
                u'replace'      : { u'__VERSIONFLT__' : u"complex_",
                                    u'__COMPLEX__'    : u'-DTYPE_COMPLEX'},
                u'database_val' : u'complex',
                u'searchInName' : re.compile(u'_complex')}
            },

        'OOC':  {
            u'database_id' : u'OOC',
            u'default'     : {
                u'replace'        : {u'__OOC__' : u''},
                u'database_val'   : u'UNDEFINED'},
            True           : {
                u'replace'        : {u'__OOC__' : u'-DOOC_FTGT'+
                                                  u' -DOOC_PERCENT_COEFTAB'+
                                                  u' -DOOC_CLOCK'+
                                                  u' -DOOC_DETECT_DEADLOCKS'},
                u'database_val'   : u'DEFINED',
                u'searchInName'   : re.compile(u'_ooc')}
            },

        u'PARTITIONER' : {
            u'database_id' : u'ORDERING',
            u'default' : {
                u'replace' : {
                    u'__PARTITIONNER__' : u'-DWITH_SCOTCH'+
                                          u' -I${SCOTCH_HOME}/include',
                    u'__LIBPART__'      : u'-L${SCOTCH_HOME}/lib'+
                                          u' -lscotch -lscotcherrexit',
                    u'__VERSIONORD__'   : u'_scotch'},
                u'database_val' : u'scotch'
                },
            u'scotch' : {
                u'replace' : {
                    u'__PARTITIONNER__' : u'-DWITH_SCOTCH'+
                                          u' -I${SCOTCH_HOME}/include',
                    u'__LIBPART__'      : u'-L${SCOTCH_HOME}/lib'+
                                          u' -lscotch -lscotcherrexit',
                    u'__VERSIONORD__'   : u'_scotch'},
                u'database_val' : u'scotch',
                u'searchInName' : re.compile(u'_scotch')
                },
            u'ptscotch' : {
                u'replace' : {
                    u'__PARTITIONNER__' : u'-DDISTRIBUTED -DWITH_SCOTCH'+
                                          u' -I${SCOTCH_HOME}/include',
                    u'__LIBPART__'      : u'-L${SCOTCH_HOME}/lib'+
                                          u' -lptscotch -lptscotcherrexit',
                    u'__VERSIONORD__'   : u'_ptscotch'},
                u'database_val' : u'ptscotch',
                u'searchInName' : re.compile(u'_ptscotch')
                },
            u'metis' : {
                u'replace' : {
                    u'__PARTITIONNER__' : u'-DMETIS -I${METIS_HOME}/include',
                    u'__LIBPART__'      : u'-L${METIS_HOME}/lib -lmetis',
                    u'__VERSIONORD__'   : u'_metis'},
                u'database_val' : u'metis',
                u'searchInName' : re.compile(u'_metis')
                }
            }
        }

    _libraries = {
        u"libpastix.a" : [u" __MAKE__ clean",
                          u" __MAKE_J__ all",
                          u" __MAKE_J__ drivers"]}
    _no_type_binaries = [u'simple', u'simple_dist', u'simple_param'] 

    _s_binaries = []
    for binary in _no_type_binaries:
        _s_binaries.append(u's%s' % binary)

    _d_binaries = []
    for binary in _no_type_binaries:
        _d_binaries.append(u'd%s' % binary)

    _c_binaries = []
    for binary in _no_type_binaries:
        _c_binaries.append(u'c%s' % binary)

    _z_binaries = []
    for binary in _no_type_binaries:
        _z_binaries.append(u'z%s' % binary)
    
    _all_binaries = []
    _all_binaries.extend(_no_type_binaries)
    _all_binaries.extend(_s_binaries)
    _all_binaries.extend(_d_binaries)
    _all_binaries.extend(_c_binaries)
    _all_binaries.extend(_z_binaries)

    for binary in _all_binaries:
        _binaries[binary] = {
            u'filename'  : binary,
            u'build_dir' : u"example/src",
            u'make_cmd'  : u"../bin/%s" % binary,
            u'binary'    : u'example/bin/%s' %  binary}
    

    _parse_util = {
        u"ORD_TIME" : {
            u'parse_key' : re.compile(u"Time to compute ordering [ ]*([-.e+0-9]*) s"),
            u'type'      : u"FLOAT"
            },
        u"ANAL_TIME" : {
            u'parse_key' : re.compile(u"Time to analyze  [ ]*([-.e+0-9]*) s"),
            u'type'      : u"FLOAT"
            },
        u"PRED_FACT_TIME" : {
            u"parse_key" : re.compile(u'Prediction Time to factorize \(IBM PWR5 ESSL\)[ ]*([-.e+0-9]*) s'),
            u"type"      : u"FLOAT"
            },
        u"FACT_TIME" : {
            u"parse_key" : re.compile(u"Time to factorize[ ]*([-.e+0-9]*) s"),
            u"type"      : u"FLOAT"
            },
        u"SOLVE_TIME" : {
            u"parse_key" : re.compile(u"Time to solve[ ]*([-.e+0-9]*) s"),
            u"type"      : u"FLOAT"
            },
        u"REF_TIME" :  {
            u"parse_key" : re.compile(u"Time for refinement[ ]*([-.e+0-9]*) s"),
            u"type"      : u"FLOAT"
            },
        u"MAX_MEM" : {
            u"parse_key" : re.compile(u"Max memory used after clean[ ]*([.e+0-9]*) ([KMG]*o)"),
            u"type"      : u"memory"
            },
        u"MEM_END" : {
            u"parse_key" : re.compile(u"Memory used after clean[ ]*([.e+0-9]*) ([KMG]*o)"),
            u"type"      : u"memory"
            },
        u"OVERHEAD" : {
            u"parse_key" : re.compile(u"Total.*overhead : ([.e+0-9]*)"),
            u"type"      : u"FLOAT"
            },
        u"MAX_OVERHEAD" : {
            u"parse_key" : re.compile(u"Maximum.*overhead : ([.e+0-9]*)"),
            u"type"      : u"FLOAT"
            },
        u"LOCAL_OVERHEAD" : {
            u"parse_key" : re.compile(u"Local.*overhead : ([.e+0-9]*)"),
            u"type"      : u"double"
            },
        u"FILLRATE2" : {
            u'parse_key' : re.compile(u"fillrate2 ([.e+0-9]*)"),
            u'type'      : u'FLOAT'
            },
        u"NB_ITER" : {
            u'parse_key' : re.compile(u"Refinement[ ]*([0-9]+) iterations"),
            u'type'      : u'SMALLINT'
            },
        u"PRECISION" : {
            u'parse_key2' : re.compile(u", norm=([0-9e\.-]+)"), # unused ?
            u'parse_key'  : re.compile(u"Precision : .* = ([0-9e\.-]+)"),
            u'type'       : u'FLOAT'
            },
        u"VERSION" : {
            u'parse_key' : re.compile(u"Version[ ]*:[ ](.*)\n"),
            u'type'      : "VARCHAR(30)"
            },
        u"MATRIX_SIZE" : {
            u'parse_key' : re.compile(u"  Matrix size  [ ]*([0-9]*) x [0-9]*"),
            u'type'      : u"INT"
            },
        u"NNZA" : {
            u'parse_key' : re.compile(u'[0-9]* : Number of nonzeros \(local block structure\)[ ]*([0-9]*)'),
            u'type'      : u"BIGINT"
            },
        u"NNZL" : {
            u'parse_key' : re.compile(u'Number of nonzeros in factorized matrice[ ]*([0-9]*)'),
            u'type'      : u"BIGINT"
            },

        u"NNZL_BLOCK" : {
            u'parse_key' : re.compile(u'Number of nonzeros (block structure)[ ]*([0-9]*)'),
            u'type'      : u"BIGINT"
            },

        u"OPC" : {
            u'parse_key' : re.compile(u'Number of operations \(L[ULt]*\) [ ]*([-.e+0-9]*)'),
            u'type'      : u"FLOAT"
            },
        u"FILL_IN" : {
            u'parse_key' : re.compile(u"Fill-in [ ]*([-.e+0-9]*)"),
            u'type'      : u"FLOAT"
            },
        u"FILL_IN_BLOCK" : {
            u'parse_key' : re.compile(u"Fill-in (block structure)[ ]*([-.e+0-9]*)"),
            u'type'      : u"FLOAT"
            },

        u"SOLVMTX_SIZE" : {
            u'parse_key' : re.compile(u'SolverMatrix size \(without coefficient\)[ ]*([-.e+0-9]*) ([KMG]*o)'), 
            u'type'      : u"memory"
            },

        u"MAX_COEFTAB_SIZE" : {
            u'parse_key' : re.compile(u'Maximum coeftab size \(cefficients\)[ ]*([-.e+0-9]*) ([KMG]*o)'),
            u'type'      : u"memory"
            },

        u"NB_TASK_ESP" : {
            u'parse_key' : re.compile(u"Number of tasks added by esp[ ]*([-0-9]*)"),
            u'type'      : u"INT"
            },
        u"NB_TASK" : {
            u'parse_key' : re.compile(u"Number of tasks  [ ]*([-0-9]*)"),
            u'type'      : u"INT"
            }
        }
    
        ### PARSE LOG COUNTERS

    _parse_count = {
        u"ERRORS" : {
            u'parse_key' : re.compile(u"ERROR:", flags=re.IGNORECASE),
            u'count'     : 0},

        u"WARNINGS" : {
            u'parse_key' : re.compile(u"WARNING:", flags=re.IGNORECASE),
            u'count'     : 0}
        }

        ## ConfigFlags
    _configFlags = {
        u'__LBLAS__' : {
            u'machine_key' : u'LIBBLAS',
            u'default'     : u'-lblas'
            },
        u'__LIBFORT__' :  {
            u'machine_key' : u'LIBFORTRAN',
            u'default'     : u'-lgfortran'
            },
        u'__LIBMATH__' : {
            u'machine_key' : u'LIBMATH',
            u'default'     : u'-lm'
            },
        u'__LIBREALTIME__' : {
            u'machine_key' : u'LIBREALTIME',
            u'default'     : u'-lrt'
            },
        u'__CFLAGS__' : {
            u'machine_key' : u'CCFLAGS',
            u'default'     : u''
            },
        u'__FFLAGS__' : {
            u'machine_key' : u'FCFLAGS',
            u'default'     : u''
            },
        u'__CC__' : {
            u'machine_key' : u'CC',
            u'default'     : u'gcc'
            },
        u'__FC__' : {
            u'machine_key' : u'FC',
            u'default'     : u'gfortran'
            },
        u'__MPCC__' : {
            u'machine_key' : u'MPCC',
            u'default'     : u'mpicc'
            },
        u'__MPFC__' : {
            u'machine_key' : u'MPFC',
            u'default'     : u'mpif90'
            },
        u'__COPTFLAGS__' : {
            u'machine_key' : u'COPT',
            u'default'     : u''
            },
        u'__AR__' : {
            u'machine_key' : u'AR',
            u'default'     : u'ar'
            },
        u'__ARFLAGS__' : {
            u'machine_key' : u'ARFLAGS',
            u'default'     : u'-crs'
            },
        u'__HOSTARCH__' : {
            u'machine_key' : u'ARCH',
            u'default'     : u'i686_pc_linux'
            },
        u'__VERSIONBIT__' : {
            u'machine_key' : u'VERSIONBIT',
            u'default'     : u'_32bit'
            },
        u'__MAKE__' : {
            u'machine_key' : u'MAKE',
            u'default'     : u'make'
            }
        }


    def Configure(self,  machine, options, compilations, ident):
        """Configure the project on the given machine"""
        config_xml = Project.Configure(self,
                                       machine,
                                       options,
                                       compilations,
                                       ident)


    # Definition de la flavor a utiliser
    #$param_exec[_PM2_FLAVOR_] = "";
    #if ( $cas =~ /.*_FL([-a-z0-9]*)FL.*/ )[
    #	$param_exec[_PM2_FLAVOR_] = "$1";
    #   ]

    # Utilisation des thread Marcels a la place des threads POSIX
    #  if ( $cas =~ /.*_Marcel.*/ || $cas =~ /.*_dynSched.*/ )[
    #	if ( $param_exec[_PM2_FLAVOR_] eq "")[
    #      print "ERROR : Flavor non definie";
    #      exit;
    #	]
    #	$param_exec[_LDFLAGS_] .= " `pm2-config --libs`";
    #       $param_exec[_CCFLAGS_] .= " `pm2-config --cflags` -I$ENV[ PM2_ROOT]/marcel/include/pthread";
    #  ] else [
    #     $param_exec[_LDFLAGS_] .= " -lpthread";
    #]

    # Utilisation des bulles avec les ordo de Marcel
    #if ( $cas =~ /.*_dynSched.*/ )[
    #    $param_exec[_DYN_]      = "VERSIONSCH  = _dyn"; 
    #    $param_exec[_CCFLAGS_] .= " -DPASTIX_DYNSCHED -DPASTIX_BUBBLESCHED";
    #]
    #else
    #[
    #    $param_exec[_DYN_] = "VERSIONSCH = _static"; 
    #]

    def Run(self, machines, run_dict, ident):
        Project.Run(self,  machines, run_dict, ident)

    def _runOne(self, binary, case, parameters, machine,
                options, revision, ident):
        launch = u""
        filename = re.compile(".*/([^/]*)")
        binary = self._binaries[binary][u'filename']
        iparms = parameters[u"iparms"]
        amag_key = u"IPARM_AMALGAMATION_LEVEL"
        lof_key  = u"IPARM_LEVEL_OF_FILL"
        ooc_key  = u'IPARM_OOC_LIMIT'
        esp_thr_key = u'IPARM_ESP_THRESHOLD'
        dist_lvl_key = u'IPARM_DISTRIBUTION_LEVEL'

        re_dynsched = re.compile(u"PASTIX_DYNSCHED")
            
        for nthread in parameters[u'THREAD']:
            for nproc in parameters[u'PROC']:
                if (nproc*nthread <= machine[u"nproc"]):
                    run_xml = xml.Element(u"Run")
                    arguments = [u'%s/%s' % (self.install, binary)]
                    if ( re.search(u"\.mm", case) or
                         re.search(u"\.mtx", case) or
                         re.search(u'\.ijv', case)):
                        arguments.append(u'-mm')
                    elif ( re.search(u"\.rua", case) or
                           re.search(u"\.rsa", case)):
                        arguments.append(u'-rsa')

                    arguments.append(case)
                    arguments.append(u'-t')
                    arguments.append(unicode(str(nthread)))

                    if iparms.has_key(u"IPARM_VERBOSE"):
                        arguments.append(u'-v')
                        arguments.append(unicode(str(iparms[u"IPARM_VERBOSE"])))
                        self._addTextChild(run_xml, u'option',
                                           str(iparms[u"IPARM_VERBOSE"]),
                                          {u'id'   : u'IPARM_VERBOSE',
                                           u'type' : u'INT'})

                    
                    if iparms.has_key(ooc_key):
                        arguments.append(u"-ooc")
                        arguments.append(unicode(str(iparms[ooc_key])))
                        self._addTextChild(run_xml, u'option',
                                           str(iparms[ooc_key]),
                                           {u'id'   : ooc_key,
                                            u'type' : u'INT'})

                    if (iparms.has_key(esp_thr_key)):
                        arguments.append(u"-iparm 43 1 -iparm 56 ")
                        arguments.append(unicode(str(iparms[esp_thr_key])))

                    if (iparms.has_key(dist_lvl_key)):
                        arguments.append(u"-iparm 35")
                        arguments.append(unicode(str(iparms[dist_lvl_key])))
                    if (re_dynsched.search(options[u"NAME"])):
                        arguments.append(u"-iparm IPARM_THREAD_COMM_MODE API_THREAD_COMM_ONE")
                    if (iparms.has_key(amag_key) and 
                        iparms.has_key(lof_key)):
                        arguments.append(u"-incomp")
                        arguments.append(unicode(str(iparms[amag_key])))
                        arguments.append(unicode(str(iparms[lof_key])))

                        am_lvl = str(iparms[amag_key])
                        self._addTextChild(run_xml, u'option',
                                           am_lvl,
                                           {u'id'   : amag_key,
                                            u'type' : u'INT'})
                        lof = str(iparms[lof_key])
                        self._addTextChild(run_xml, u'option',
                                           lof,
                                           {u'id'   : lof_key,
                                            u'type' : u'INT'})

                    self._addTextChild(run_xml, u"BINARY_NAME",  binary)
                    date = self._idToDate(ident)
                    self._addTextChild(run_xml, u"PROJECT",      self._name)
                    self._addTextChild(run_xml, u"DATE",         date)
                    self._addTextChild(run_xml, u"MD5CONF",      self.md5conf)
                    self._addTextChild(run_xml, u"MD5DIFF",      self.md5diff)
                    self._addTextChild(run_xml, u"OPT_SET_NAME",
                                       parameters[u"NAME"])
                    mycase = filename.search(case)
                    if mycase: mycase = mycase.group(1)
                    else     : mycase = case
                    self._addTextChild(run_xml, u"CASE_NAME",    mycase)
                    self._addTextChild(run_xml, u"MPI_PROC_NBR",
                                       unicode(str(nproc)))
                    self._addTextChild(run_xml, u"THREAD_NBR",
                                       unicode(str(nthread)))
                    self._addTextChild(run_xml, u"USER",         self.user)
                    self._addTextChild(run_xml, u"MACHINE",
                                       machine["NAME"])
                    self._addTextChild(run_xml, u"BUILD",
                                       options["NAME"])
                    self._addTextChild(run_xml, u"REVISION",     revision)

                    name = u"-".join([parameters["NAME"],
                                      unicode(os.path.basename(case)),
                                      binary,
                                      unicode(str(nproc)),
                                      unicode(str(nthread))])
                    script = machine["script"].replace(u"__PROCNBR__",
                                                       str(nproc))
                    script = script.replace(u"__THRDNBR__", str(nproc))
                    script = script.replace(u'__JOBID__', ident)
                    script = script.replace(u"__CMD__", u" ".join(arguments))
                    script = script.replace(u"__DIR__", self.install)
                    script = script.replace(u"__LOG__", name+u".log")

                    self._addTextChild(run_xml, u"logfile",      
                                       u'%s/%s.log' % (self.install, name))
                    self._addTextChild(run_xml, u"SCRIPT",
                                       script)
                    self.doc_xml.append(run_xml)

                    f=open(u'%s/%s.sh' % (self.install, name) , u'w')
                    f.write(script)
                    f.close()
                    launch += u'%s %s/%s.sh\n' % (machine[u"soumission"],
                                                  self.install,
                                                  name)
        return launch

    def getRevision(self):
        revision = Project.getRevision(self)
        if revision == -1:
            if os.path.isfile(u'%s/Revision' %(self.branchdir)):
                revision = int(open(u'%s/Revision' %(self.branchdir),
                                    u'r').read())
        return revision

    def export(self, archive, rootdir, branches, ident):
        export_path = u'/tmp/regression-local-%s/ricar' % ident
        rev_search = re.compile('r([0-9]+)')
        subprocess.Popen([u"rm", u"-rf", export_path], env=self._env).wait()

        if __debug__:
            print u"Creating archive from %s %s" % (rootdir, str(branches))

        for branchname in branches:
            branchname = re.sub("^/", "", branchname)
            branchname = re.sub("/$", "", branchname)
            subprocess.Popen([u"mkdir", u"-p",
                              os.path.dirname(os.path.join(export_path,
                                                           branchname))],
                             env=self._env, stdout=open('/dev/null',
                                                        'w')).wait()
            diffFile    = os.path.join( os.path.join(export_path, branchname),
                                        u'__DIFF__' )

            cmd = [u'rm', u'-rf',
                   os.path.join(export_path, branchname)]
            if __debug__:
                print str(cmd)
            subprocess.Popen(cmd, env=self._env).wait()
            if (re.search(u'git\+ssh', rootdir)):
                if __debug__:
                    print u"Downloading last revision"
                directory = os.path.join(export_path,branchname)
                cmd = [u'git', u'clone', rootdir,
                       os.path.join(export_path,branchname)]
                process = subprocess.Popen(cmd, env=self._env,
                                           stdout=open('/dev/null', 'w')).wait()
                if branchname == u"master":
                    cmd = [u'git', u'checkout', branchname]
                else:
                    cmd = [u'git', u'checkout', '-b', branchname,
                           u'origin/%s' %(branchname)]
                process = subprocess.Popen(cmd, env=self._env,
                                           stdout=open('/dev/null', 'w'),
                                           cwd=directory).wait()
                cmd = [u'git', u'rev-parse', 'HEAD']
                p = subprocess.Popen(cmd, env=self._env,
                                     stdout=subprocess.PIPE,
                                     cwd=directory)
                
                revision = p.communicate()[0].strip()

                cmd = [u'touch', 'config.in']
                p = subprocess.Popen(cmd, env=self._env,
                                     stdout=subprocess.PIPE,
                                     cwd=directory)
                p.wait()

                cmd = [u'make', 'murge']
                p = subprocess.Popen(cmd, env=self._env,
                                     stdout=subprocess.PIPE,
                                     cwd=directory)
                p.wait()

                git_diff = u''

            else:
                cmd = [u'cp', u'-rf',
                       os.path.join(rootdir),
                       os.path.join(export_path, branchname)]
                if __debug__:
                    print str(cmd)
                subprocess.Popen(cmd, env=self._env,
                                 cwd=os.path.dirname(export_path)).wait()
                cmd = [u'make', u'clean']
                if __debug__:
                    print str(cmd)
                subprocess.Popen(cmd, env=self._env,
                                 cwd=os.path.join(export_path,
                                                  branchname)).wait()
                cmd = [u'rm', u'-rf',
                       os.path.join(export_path,
                                    branchname, self._config_name),
                       os.path.join(export_path, branchname, '*/bin'),
                       os.path.join(export_path, branchname, '*/obj')]
                if __debug__:
                    print str(cmd)
                subprocess.Popen(u' '.join(cmd), shell=True,
                                 env=self._env,
                                 cwd=os.path.dirname(export_path)).wait()

                if os.path.isdir(os.path.join(rootdir, u'.git')):
                    try:
                        p_git_diff = subprocess.Popen([u'git', u'diff'],
                                                      env=self._env,
                                                      cwd=rootdir,
                                                      stdout=subprocess.PIPE)
                        git_diff = p_git_diff.communicate()[0]
                        git_diff = git_diff.decode(u'utf8', u'replace')
                    except:
                        git_diff = u'error in  git diff'
                    try:
                        cmd = [u'git', u'rev-parse', 'HEAD']
                        directory = os.path.join(export_path,branchname)
                        p = subprocess.Popen(cmd, env=self._env,
                                             stdout=subprocess.PIPE,
                                             cwd=directory)
                        revision = p.communicate()[0].strip()
                    except:
                        revision = -1


            f=open(diffFile, u'w')
            f.write(git_diff)
            f.close()
            f = open(os.path.join( os.path.join(export_path,
                                                branchname) ,
                                   u"Revision"), u"w")
            f.write(str(revision))
            f.close()
            f = open(os.path.join( os.path.join(export_path, branchname) ,
                                   u"__REVISION__"), u"w")
            f.write(str(revision))
            f.close()

        if branches:
            cmd = [u'tar', u'-czf',
                   archive,
                   os.path.basename(export_path)]
            
            subprocess.Popen(cmd, env=self._env,
                             cwd=os.path.dirname(export_path)).wait()

    def _checkError(self):
        """ Correct warning or errors number 
        
        If MEM_END != 0 increase warnings
        If NB_ITER > 10 we increase the number of warnings.
        If NB_ITER = 250 we increase the number of errors.
        If NB_ITER is not defined we increase the number of errors.

        """
        list_count = [u'ERRORS', u'WARNINGS']
        list_var   = {
            u'MEM_END' :  {
                None : None,
                'not' : {u"0" : u'WARNINGS'},
                'is'  : {},
                'gt'  : {},
                'lt'  : {}}, 
            u'NB_ITER' : {
                None : u'ERRORS',
                'not' : {},
                'is'  : {u"250"  : u'ERRORS'},
                'gt'  : {u"10"   : u'WARNINGS'},
                'lt'  : {} 
                }
            }

        list_nodes = {}
        for run_xml in self.doc_xml.getiterator(u'Run'):
            if ( run_xml.find(u'BINARY_NAME').text in self._s_binaries or
                 run_xml.find(u'BINARY_NAME').text in self._c_binaries ):
                list_var[u'PRECISION'] = {
                    None : u'ERRORS',
                    'not' : {},
                    'is'  : {},
                    'gt'  : {"5e-6"   : u'WARNINGS'},
                    'lt'  : {} 
                    }
            elif ( run_xml.find(u'BINARY_NAME').text in self._d_binaries or
                   run_xml.find(u'BINARY_NAME').text in self._z_binaries ):
                list_var[u'PRECISION'] =  {
                    None : u'ERRORS',
                    'not' : {},
                    'is'  : {},
                    'gt'  : {"1e-12"   : u'WARNINGS'},
                    'lt'  : {} 
                    }
            else:
                list_var[u'PRECISION'] = {
                    None : u'ERRORS',
                    'not' : {},
                    'is'  : {},
                    'gt'  : {"5e-6"   : u'WARNINGS'},
                    'lt'  : {} 
                    }

            for key in list_count:
                list_nodes[key] = None
            for key in list_var.keys():
                list_nodes[key] = None

            #Search for keys in list_count/list_var and update 
            # list_nodes
            for output in run_xml.getiterator(u'OUTPUT'):
                
                for key in list_count:
                    if output.attrib[u'key'] == key:
                        list_nodes[key] = output
                for key in list_var.keys():
                    if output.attrib[u'key'] == key:
                        list_nodes[key] = output

            re_small = re.compile(u"small.rsa")
            re_struc = re.compile(u"structural zeros found on the diagonal.")
            if ( re_small.match(run_xml.find(u'CASE_NAME').text) and
                 list_nodes[u'WARNINGS'] != None):
                log = run_xml.find(u'LOG')
                if (log != None):
                    log = log.text
                    list_nodes[u'WARNINGS'].text = str(
                        int(list_nodes[u'WARNINGS'].text) -
                        len(re_struc.findall(log)))

            for key in list_var.keys():
                increase = None
                if list_nodes[key] == None: 
                    increase = list_var[key][None]
                else:
                    for key2 in list_var[key][u'not'].keys():
                        if key2 != list_nodes[key].text:
                            increase = list_var[key][u'not'][key2]
                    
                    for key2 in list_var[key][u'is'].keys():
                        if key2 == list_nodes[key].text:
                            increase = list_var[key][u'is'][key2]

                    for key2 in list_var[key][u'gt'].keys():
                        if float(key2) < float(list_nodes[key].text):
                            increase = list_var[key][u'gt'][key2]

                    for key2 in list_var[key][u'lt'].keys():
                        if float(key2) > float(list_nodes[key].text):
                            increase = list_var[key][u'lt'][key2]

                if increase != None:
                    if list_nodes[increase] == None:
                        child = xml.Element(u"OUTPUT")
                        child.attrib[u"key"]  = increase
                        child.attrib[u"type"] = u"INT"
                        child.text = u'1'
                        run_xml.append(child)
                    else:
                        newval = str(int(list_nodes[increase].text) + 1)
                        list_nodes[increase].text = newval
                
