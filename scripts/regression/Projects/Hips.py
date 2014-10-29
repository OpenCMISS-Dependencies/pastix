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



class Hips(Project):
    """Class defining PaStiX project."""

    _config =  u"""
###
###  HIPS specific compilation flags
###

##   Arithmetic
##    - default is -DTYPE_REAL in double precision (-DPREC_DOUBLE)
##    - use -DTYPE_COMPLEX to build Complex version of HIPS
##    - use -DPREC_SIMPLE to compute in single precision

COEFTYPE     = __PRECISION__ __COMPLEX__

#COEFTYPE    = -DTYPE_REAL
#COEFTYPE    = -DTYPE_COMPLEX

#COEFTYPE    = -DTYPE_REAL    -DPREC_SIMPLE
#COEFTYPE    = -DTYPE_COMPLEX -DPREC_SIMPLE


##   Partitionner
##    - default partitioner is METIS
##    - use -DSCOTCH_PART to use SCOTCH

PARTITIONER  = __ORDERING__
#PARTITIONER = -DSCOTCH_PART

##   Integer size
##    - default int type is    : INTS = INTL = int (C type length)
##    - use -DINTSIZE32 to set : INTS = INTEGER(4) and INTL = INTEGER(4)
##    - use -DINTSIZE64 to set : INTS = INTEGER(4) and INTL = INTEGER(8)

INTSIZE      = __INTSIZE__
#INTSIZE     = -DINTSIZE64
#INTSIZE     = -DINTSIZE32


###
###  Compiler
###

ARCH       = __ARCH__

CC         = __CC__    # C compiler
MPICC      = __MPCC__
FC         = __FC__    # Fortran compiler
MPIFC      = __MPFC__
LD         = $(FC)     # Linker
MPILD      = $(MPIFC)

CFLAGS     = __CFLAGS__ __OPT__# Additional C compiler flags
FFLAGS     = __FFLAGS__# Additional Fortran compiler flags
LFLAGS     = __LFLAGS__# Additional linker flags




COPTFLAGS  = __COPTFLAGS__  # Optimization flags
FOPTFLAGS  = __FOPTFLAGS__  #

###
###  Library
###

IBLAS      = __IBLAS__    # BLAS include path
LBLAS      = __LBLAS__    # BLAS linker flags

IMPI       = __IMPI__     # Additional MPI include path
LMPI       = __LMPI__     # Additional MPI linker flags

##   METIS_DIR : path to METIS
METIS_DIR  = __METIS_DIR__
IMETIS     = -I$(METIS_DIR)/include
LMETIS     = -L$(METIS_DIR)/lib -lmetis

##   SCOTCH_DIR : path to SCOTCH
SCOTCH_DIR = __SCOTCH_DIR__
ISCOTCH    = -I$(SCOTCH_DIR)/include
LSCOTCH    = -L$(SCOTCH_DIR)/lib -lscotch -lscotcherr

###
###  Misc
###

MAKE       = __MAKE__
AR         = __AR__
ARFLAG     = __ARFLAGS__ # -crs
LN         = __LN__      # ln
CP         = __CP__      # cp

###

##   Uncomment that to run in DEBUG mode
#DEBUG     = -g -DDEBUG_M
"""
    _inputs = """__MATRIX_NAME__ #0#  Matrix name and driver
__ISSYM__       #1#  0=unsymmetric pattern, 1=symmetric pattern, 2=symmetric matrix in RSA (lower triangular part stored)
__RHSNAME__     #2#  RHS file name (0=make a rhs (sol = [1]))
__METHOD__      #3#  Method (HYBRID, ITERATIVE)
__PREC__        #4#  Relative residual norm
__FILLIN__      #5#  Fill-in : Strictly=0, Locally=ALL
__MAXITER__     #6#  GMRES maximum iteration
__RESTART__     #7#  GMRES restart
__THRES_INT__   #8#  Numerical threshold in ILUT for interior domain (important : set 0.0 in HYBRID)
__THRES_SCHUR__ #9#  Numerical threshold in ILUT for Schur preconditioner
__THRES_COUPL__ #10# Numerical threshold for coupling between the interior level and Schur
__VERBOSE__     #11# Verbose level [0-4] (the higher the more it prints informations)
"""

    _config_name = u"makefile.inc"
    _name        = u'Hips'
    _libraries = {
        u"libhips.a" : [u"__MAKE__ clean",
                        u"__MAKE__ all"]}
    _optionsComp = {
        u'INTEGER' : {
            u'database_id' : u'INTEGER',
            u'default' : {
                u'replace'      : { u'__INTSIZE__': u'' },
                u'database_val' : u'int'
                },
            u'int32'   : {
                u'replace'      : { u'__INTSIZE__': u'-DINTSIZE32' },
                u'database_val' : u'int32',
                u'searchInName' : re.compile(u'_int32')
                },
            u'int64'   : {
                u'replace'      : { u'__INTSIZE__': u'-DINTSIZE64' },
                u'database_val' : u'int64',
                u'searchInName' : re.compile(u'_int64')
                }
            },

        u'PRECISION' : {
            u'database_id' : u'PRECISION',
            u'default'     : {
                u'replace'      : { u'__PRECISION__' : u''},
                u'database_val' : u'double'},
            u'simple'      : {
                u'replace'      : { u'__PRECISION__' : u'-DPREC_SIMPLE'},
                u'database_val' : u'simple',
                u'searchInName' : re.compile(u'_simple')}
            },

        u'COMPLEX'  : {
            u'database_id'  : u'COMPLEX',
            u'default'      : {
                u'replace'      : { u'__COMPLEX__'  : u""},
                u'database_val' : u'real'},
            True            : {
                u'replace'      : { u'__COMPLEX__'  : u"-DTYPE_COMPLEX"},
                u'database_val' : u'complex',
                u'searchInName' : re.compile(u'_complex')}
            },

        u'ORDERING' : {
            u'database_id'  : u'ORDERING',
            u'default'      : {
                u'replace'      : { u'__ORDERING__'  : u""},
                u'database_val' : u'metis'},
            True            : {
                u'replace'      : { u'__ORDERING__'  : u"-DSCOTCH_PART"},
                u'database_val' : u'scotch',
                u'searchInName' : re.compile(u'_scotch')}
            }
        }

    _configFlags = {        
        u'__IBLAS__' : {
            u'machine_key' : u'INCBLAS',
            u'default'     : u''
            },
        u'__LBLAS__' : {
            u'machine_key' : u'LIBBLAS',
            u'default'     : u'-lblas'
            },
        u'__LMPI__' : {
            u'machine_key' : u'LIBMPI',
            u'default'     : u''
            },
        u'__IMPI__' : {
            u'machine_key' : u'INCMPI',
            u'default'     : u''
            },
        u'__LFLAGS__' : {
            u'machine_key' : u'LDFLAGS',
            u'default'     : u''
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
        u'__FOPTFLAGS__' : {
            u'machine_key' : u'FOPT',
            u'default'     : u''
            },
        u'__AR__' : {
            u'machine_key' : u'AR',
            u'default'     : u'ar'
            },
        u'__MAKE__' : {
            u'machine_key' : u'MAKE',
            u'default'     : u'make'
            },
        u'__ARFLAGS__' : {
            u'machine_key' : u'ARFLAGS',
            u'default'     : u'-crs'
            },
        u'__ARCH__' : {
            u'machine_key' : u'ARCH_HIPS',
            u'default'     : u'-DLINUX'
            },
        u'__LN__' : {
            u'machine_key' : u'LN',
            u'default'     : u'ln'
            },
        u'__CP__' : {
            u'machine_key' : u'CP',
            u'default'     : u'cp'
            },
        }

    _parse_count = {
        u"ERRORS" : {
            u'parse_key' : re.compile(u"ERROR:", flags=re.IGNORECASE),
            u'count'     : 0
            },
        u"WARNINGS" : {
            u'parse_key' : re.compile(u"WARNING:", flags=re.IGNORECASE),
            u'count'     : 0}
        }

    def __init__(self):
        """ Initialise Hips"""

        ## Binaries
        BINARIES = {}
        BINARIES[u'DBMATRIX'] = [
            {
                u'names'      : [u'testSAVE.ex',
                                 u'testDirect.ex',
                                 u'testDirect2.ex'],
                u'parameters' : []},
            {
                u'names'      : [u'testCut.ex',
                                 u'testDBMatrix.ex',
                                 u'testDBMatrixGEMM1',
                                 u'testDBMatrixGEMM2',
                                 u'testDBMatrixLOAD.ex',
                                 u'testDirectPhidal.ex',
                                 u'testFillIn.ex',
                                 u'testFind.ex',
                                 u'testGEMM.ex'],
                u'parameters' : [u'domsize']}]

        BINARIES[u'DEBUG_PARALLEL'] = [
            {
                u'names'     : [u'testHIPS_DEBUG.ex'],
                u'parameters': [u'domsize']}]

        BINARIES[u'GRID'] = [
            {
                u'names'     : [u'testGRID.ex',
                                u'testILUT_GRID.ex'],
                u'parameters': [u'domsize']}]

        BINARIES[u'MISC_DEVEL'] = [
            {
                u'names'      : [u'testGEMM.ex'],
                u'parameters' : [u'ndom', u'ndom2']},
            {
                u'names'      : [u'testICCT_GRID.ex'],
                u'parameters' : [u'domsize', u'levelfor']},
            {
                u'names'      : [u'testICCT_ML.ex',
                                 u'testILUT_ML.ex'],
                u'parameters' : [u'ndom', u'levelfor']},
            {
                u'names'      : [u'testILUK.ex'],
                u'parameters' : [u'levelk']},
            {
                u'names'      : [u'testILUT_GRID.ex',
                                 u'testPHIDAL.ex'],
                u'parameters' : [u'domsize']},
            {
                u'names'      : [u'testPHIDAL_DD.ex'],
                u'parameters' : [u'ndom']},
            {
                u'names'      : [u'testIOHB.ex'],
                u'parameters' : [u'matrix', u'matrix2', u'type']},
            {
                u'names'      : [u'testDrivers.ex'],
                u'parameters' : []}]

        BINARIES[u'MISC_PARALLEL']  = [
            {
                u'names'      : [u'testDBDistrMatrix.ex',
                                 u'testALLREAD_DBDistrMatrix.ex',
                                 u'testALLREAD_BROADCAST_HIPS.ex',
                                 u'testBROADCAST_HIPS.ex',
                                 u'testUSERPART_HIPS.ex'],
                u'parameters' : [u'domsize']},
            {
                u'names'      : [u'GenereProcessorsLocalData.ex'],
                u'parameters' : [u'domsize', u'nproc']},
            {
                u'names'      : [u'testDBDistrMatrix_debug.ex',
                                 u'testALLREAD_HIPS.ex'],
                u'parameters' : []}]

        BINARIES[u'ORDERING'] = [
            {
                u'names'      : [u'testSizedDomains.ex'],
                u'parameters' : [u'maxdomsize']},
            {
                u'names'      : [u'testXPartition.ex'],
                u'parameters' : [u'']},
            {
                u'names'      : [u'testORDERING.ex'],
                u'parameters' : [u'ndom']}]

        BINARIES[u'PARALLEL'] = [
            {
                u'names'      : [u'testHIPS1-Fortran.ex',
                                 u'testHIPS2-Fortran.ex',
                                 u'testHIPS3-Fortran.ex',
                                 u'testHIPS.ex',
                                 u'testHIPS1.ex',
                                 u'testHIPS2.ex',
                                 u'testHIPS3.ex',
                                 u'testHIPS4.ex',
                                 u'testHIPS-RUN.ex',
                                 u'testHIPS-RUN2.ex',
                                 u'testHIPS-NR.ex',
                                 u'testHIPS-DEBUG.ex',
                                 u'testHIPS-Assembly.ex',
                                 u'testHIPS-CheckSol.ex',
                                 u'testHIPS-Mem.ex',
                                 u'testHIPS2-Mem.ex'],
                u'parameters' : [u'domsize']},
            {
                u'names'      : [u'testHIPS-Save.ex'],
                u'parameters' : [u'nbpproc', u'domsize']},
            {
                u'names'      : [u'testHIPS-Grid.ex'],
                u'parameters' : [u'domsize', u'nc', u'cube']},
            {
                u'names'      : [u'testHIPS-Laplace1-Fortran.ex',
                                 u'testHIPS-Laplace2-Fortran.ex',
                                 u'testHIPS4-Fortran.ex',
                                 u'testHIPS-Load.ex'],
                u'parameters' : []}]

        BINARIES[u'PETSc'] = [
            {
                u'names'      : [u'testPETSc-LOAD.ex',
                                 u'testPETSc-Mtx.ex',
                                 u'testPETSc-GRID.ex',
                                 u'testPETSc-Mtx-SAVE.ex',
                                 u'testPETSc-GRID-SAVE.ex'],
                u'parameters' : []}]

        BINARIES[u'SEQUENTIAL'] = [
            {
                u'names'      : [u'testPHIDAL.ex'],
                u'parameters' : [u'domsize']}]

        BINARIES[u'SOLVERMATRIX'] = [
            {
                u'names'      : [u'testDirect.ex'],
                u'parameters' : []}]

        for dirname in BINARIES:
            for nameparam in BINARIES[dirname]:
                for name in nameparam[u'names']:
                    self._binaries[u'%s_%s' % (dirname, name)] = {
                        u'filename'  : name,
                        u'build_dir' : u'TESTS/%s' % dirname,
                        u'make_cmd'  : name,
                        u'binary'    : u'TESTS/%s/%s' % (dirname,name),
                        u'parameters': nameparam[u'parameters']}

        Project.__init__(self)
        return

    def Configure(self,  machine, options, compilations, ident):
        """Configure the project on the given machine"""
        config_xml = Project.Configure(self,  machine, options,
                                       compilations, ident)

    def Run(self, machines, run_dict, ident):
        Project.Run(self,  machines, run_dict, ident)

    def _runOne(self, binary, case, parameters,
                machine, options, revision, ident):

        launch = u""
        filename = re.compile(".*/([^/]*)")
        binary = self._binaries[binary][u'filename']

        for nproc in parameters[u'PROC']:
            run_xml = xml.Element(u"Run")
            arguments = [u'%s/%s' % (self.install, binary)]
            my_inputs = self._inputs.replace(u'__MATRIX_NAME__', case)
            for key in parameters[u'InputsKeywords'].keys():
                value = parameters[u'InputsKeywords'][key]
                my_inputs = my_inputs.replace(key, str(value))
                self._addTextChild(run_xml, u'option',
                                   str(value),
                                   { u'id'   : key.replace(u'__',''),
                                     u'type' : 'CHAR(20)'})

            self._addTextChild(run_xml, u"BINARY_NAME",  binary)
            date = self._idToDate(ident)
            self._addTextChild(run_xml, u"PROJECT",      self._name)
            self._addTextChild(run_xml, u"DATE",         date)
            self._addTextChild(run_xml, u"MD5CONF",      self.md5conf)
            self._addTextChild(run_xml, u"MD5DIFF",      self.md5diff)
            self._addTextChild(run_xml, u"OPT_SET_NAME", parameters[u"NAME"])
            mycase = filename.search(case)
            if mycase: mycase = mycase.group(1)
            else : mycase = case
            self._addTextChild(run_xml, u"CASE_NAME",    mycase)
            self._addTextChild(run_xml, u"MPI_PROC_NBR", unicode(str(nproc)))
            self._addTextChild(run_xml, u"THREAD_NBR",   u'1')
            self._addTextChild(run_xml, u"USER",         self.user)
            self._addTextChild(run_xml, u"MACHINE",      machine["NAME"])
            self._addTextChild(run_xml, u"BUILD",        options["NAME"])
            self._addTextChild(run_xml, u"REVISION",     revision)


            name = u"-".join([parameters["NAME"],
                              unicode(os.path.basename(case)),
                              unicode(str(nproc))])
            script = machine["script"].replace(u"__PROCNBR__", str(nproc))
            script = script.replace(u"__THRDNBR__", str(1))

            script = script.replace(u"__CMD__", u" ".join(arguments))
            script = script.replace(u"__DIR__", self.install)
            script = script.replace(u"__LOG__", name+u".log")

            self._addTextChild(run_xml, u"logfile",
                               u'%s/%s.log' % (self.install,
                                               name))
            self.doc_xml.append(run_xml)

            f=open(u'%s/%s.sh' % (self.install, name) , u'w')
            f.write(script)
            f.close()
            f=open(u'%s/Inputs' % (self.install) , u'w')
            f.write(my_inputs)
            f.close()
            launch += '%s %s/%s.sh\n' % (machine[u"soumission"],
                                         self.install, name)
        return launch


    def export(self, archive, rootdir, branches, ident):
        export_path = u'/tmp/regression-local-%s/hips' % ident
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
                             env=self._env, stdout=open('/dev/null', 'w')).wait()
            diffFile    = os.path.join( os.path.join(export_path, branchname) , u'__DIFF__' )
            revFile     = os.path.join( os.path.join(export_path, branchname) , u"__REVISION__")
            cmd = [u'rm', u'-rf',
                   os.path.join(export_path, branchname)]
            if __debug__:
                print str(cmd)
            subprocess.Popen(cmd, env=self._env).wait()
            if (re.search(u'svn\+ssh', rootdir)):
                if __debug__:
                    print u"Downloading last revision"
                cmd = [u'svn', u'export',
                       u'%s/%s' %(rootdir, branchname),
                       os.path.join(export_path,branchname)]
                process = subprocess.Popen(cmd, env=self._env, stdout=open('/dev/null', 'w')).wait()
                cmd = [u'svn', u'info', u'%s/%s' %(rootdir, branchname)]
                p1 = subprocess.Popen(cmd, env=self._env, stdout=subprocess.PIPE)
                p2 = subprocess.Popen([u'grep', u'Revision:'],env=self._env, stdin=p1.stdout, stdout=subprocess.PIPE)
                p3 = subprocess.Popen([u'awk', u'{ print $2}'],env=self._env, stdin=p2.stdout, stdout=subprocess.PIPE)
                revision = p3.communicate()[0]
                svn_diff = u''

            else:

                cmd = [u'cp', u'-rf',
                       os.path.join(rootdir, branchname),
                       os.path.join(export_path, branchname)]
                if __debug__:
                    print str(cmd)
                subprocess.Popen(cmd, env=self._env, cwd=os.path.dirname(export_path)).wait()
                cmd = [u'make', u'clean']
                if __debug__:
                    print str(cmd)
                subprocess.Popen(cmd, env=self._env, cwd=os.path.join(export_path, branchname)).wait()
                cmd = [u'rm', u'-rf',
                       os.path.join(export_path, branchname, self._config_name)]
                if __debug__:
                    print str(cmd)
                subprocess.Popen(u' '.join(cmd), shell=True, env=self._env, cwd=os.path.dirname(export_path)).wait()

                if os.path.isdir(os.path.join(rootdir, u'.svn')):
                    try:
                        svn_diff = subprocess.Popen([u'svn', u'diff'],
                                                    env=self._env,
                                                    cwd=rootdir,
                                                    stdout=subprocess.PIPE).communicate()[0]
                        svn_diff = svn_diff.decode(u'utf8', u'replace')
                    except:
                        svn_diff = u'error in  svn diff'

                    try:
                        svn_log = subprocess.Popen([u'svn', u'log', u'-l', u'1'],
                                                   env=self._env,
                                                   cwd=rootdir,
                                                   stdout=subprocess.PIPE,
                                                   stderr=subprocess.STDOUT).communicate()[0]

                        svn_log = svn_log.decode(u'utf8',
                                                 u'replace')
                        res = rev_search.search(svn_log)
                        if (res):
                            revision = int(res.group(1))
                    except:
                        revision = -1

                else:
                    
                    try:
                        svn_diff = subprocess.Popen([u'git',u'diff',
                                                     u'remotes/git-svn', u'.'],
                                                    env=self._env,
                                                    cwd=os.path.join(rootdir, branchname),
                                                    stdout=subprocess.PIPE).communicate()[0]
                        svn_diff = svn_diff.decode(u'utf8', u'replace')
                    except:
                        svn_diff = u"Could not get diff"
                    try:
                        svn_log = subprocess.Popen([u'git', u'svn',
                                                    u'log', u'-n', u'1'],
                                                   env=self._env,
                                                   cwd=os.path.join(rootdir,
                                                                    branchname),
                                                   stdout=subprocess.PIPE,
                                                   stderr=subprocess.STDOUT).communicate()[0]
                        svn_log = svn_log.decode(u'utf8', u'replace')
                        res = rev_search.search(svn_log)
                        if (res):
                            revision = int(res.group(1))
                    except:
                        revision = -1

            f=open(diffFile, u'w')
            f.write(svn_diff)
            f.close()
            f = open(revFile, u"w")
            f.write(str(revision))
            f.close()

        cmd = [u'tar', u'-czf',
               archive,
               os.path.basename(export_path)]
        subprocess.Popen(cmd, env=self._env,
                         cwd=os.path.dirname(export_path)).wait()

    def getText(nodelist):
        """ Return text contained in a xml nodelist"""
        rc = []
        for node in nodelist:
            if node.nodeType == node.TEXT_NODE:
                rc.append(unicode(node.data))
        return (''.join(rc)).replace("'","\\'")

    def _checkError(self):
        return
        #runs = self.doc_xml.getElementsByTagName(u"Run")
        #for run_xml in runs:
        #    try:
        #        nb_iter = int(getText(run_xml.getElementsByTagName(u"NB_ITER")[0].childNodes))
        #        if nb_iter == 250:
        #            try:
        #                nb_warn = int(getText(run_xml.getElementsByTagName(u"WARNINGS")[0].childNodes))
        #                run_xml.getElementsByTagName(u"WARNINGS")[0].childNodes[0].data = str(nb_warn + 1)
        #            except:
        #                self.__addTextChild(run_xml, u"WARNINGS", u'1')
        #    except:
        #        try:
        #            nb_errors = int(getText(run_xml.getElementsByTagName(u"ERRORS")[0].childNodes))
        #            run_xml.getElementsByTagName(u"ERRORS")[0].childNodes[0].data = str(nb_errors + 1)
        #        except:
        #            self.__addTextChild(run_xml, u"ERRORS", u'1')
