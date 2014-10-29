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



class Scotch(Project):
    """Class defining Scotch project."""

    _liste_default_yes = [u'COMMON_FILE_COMPRESS_GZ',
                          u'COMMON_PTHREAD',
                          u'COMMON_RANDOM_FIXED_SEED',
                          u'SCOTCH_RENAME',
                          u'SCOTCH_RENAME_PARSER',
                          u'SCOTCH_PTHREAD']
    _liste_default_not = []

    _config =  u"""
EXE		=
LIB		= .a
OBJ		= .o

MAKE		= __MAKE__
AR		= __AR__
ARFLAGS		= __ARFLAGS__
CAT		= cat
CCS		= __MPCC__
CCP		= __MPCC__
CCD		= __MPCC__
CFLAGS		= __CFLAGS__ __COPTFLAGS__ __OPT__ __INTSIZE__ -Drestrict=__RESTRICT__
CLIBFLAGS	=
LDFLAGS		= __LFLAGS__ __LIBGZ__ -lm -lrt
CP		= __CP__
LEX		= flex -Pscotchyy -olex.yy.c
LN		= __LN__
MKDIR		= mkdir
MV		= mv
RANLIB		= ranlib
YACC		= bison -pscotchyy -y -b y
"""

    _config_name = u"Makefile.inc"
    _name        = u'Scotch'
    _libraries = {
        u"libscotch.a" : [u"__MAKE__ clean",
                          u"__MAKE__ scotch"],
        u"libptscotch.a" : [u"__MAKE__ clean",
                            u"__MAKE__ ptscotch"]}
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
            }
        }
    _configFlags = {
        u'__LFLAGS__' : {
            u'machine_key' : u'LDFLAGS',
            u'default'     : u''
            },
        u'__CFLAGS__' : {
            u'machine_key' : u'CCFLAGS',
            u'default'     : u''
            },
        u'__CC__' : {
            u'machine_key' : u'CC',
            u'default'     : u'gcc'
            },
        u'__COPTFLAGS__' : {
            u'machine_key' : u'COPT',
            u'default'     : u''
            },
        u'__MPCC__' : {
            u'machine_key' : u'MPCC',
            u'default'     : u'mpicc'
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

    _binaries = {}
    _parse_util = {}

    _parse_count = {
        u"ERRORS" : {
            u'parse_key' : re.compile(u"ERROR:", flags=re.IGNORECASE),
            u'count'     : 0},
        
        u"WARNINGS" : {
            u'parse_key' : re.compile(u"WARNING:", flags=re.IGNORECASE),
            u'count'     : 0}
        }
    
    def __init__(self):
        Project.__init__(self)
        # TODO: Pour bien faire il faudrait ajouter la libgz à la machine.
        self._optionsComp[u"COMMON_FILE_COMPRESS_GZ"][u'default'][u'replace']['__LIBGZ__'] = u'-lz'
        self._optionsComp[u"COMMON_FILE_COMPRESS_GZ"][False][u'replace']['__LIBGZ__'] = u''
        
        self._optionsComp[u'RESTRICT'] = {
            u'database_id' : u'RESTRICT',
            u'default'     : {
                u'replace'       : {u'__RESTRICT__' : u'__restrict'},
                u'database_val'  : u'__restrict'},
            False          : {
                u'replace'       : {u'__RESTRICT__' : u''},
                u'database_val'  : u'',
                u'searchInName'  : re.compile(u'_NO_RESTRICT')}
            }

    def Configure(self,  machine, options, compilations, ident):
        """Configure the project on the given machine"""
        config_xml = Project.Configure(self,
                                       machine,
                                       options,
                                       compilations,
                                       ident)

    def Run(self, machines, run_dict, ident):
        Project.Run(self,  machines, run_dict, ident)

    def _runOne(self, binary, case, parameters,
                machine, options, revision, ident):
        return launch


    def export(self, archive, rootdir, branches, ident):
        export_path = u'/tmp/regression-local-%s/scotch' % ident
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
            if (re.search(u'svn\+ssh', rootdir)):
                if __debug__:
                    print u"Downloading last revision"
                cmd = [u'svn', u'export',
                       u'%s/%s' %(rootdir, branchname),
                       os.path.join(export_path,branchname)]
                process = subprocess.Popen(cmd, env=self._env,
                                           stdout=open('/dev/null', 'w')).wait()
                cmd = [u'svn', u'info', u'%s/%s' %(rootdir, branchname)]
                p1 = subprocess.Popen(cmd, env=self._env,
                                      stdout=subprocess.PIPE)
                p2 = subprocess.Popen([u'grep', u'Revision:'],
                                      env=self._env,
                                      stdin=p1.stdout, stdout=subprocess.PIPE)
                p3 = subprocess.Popen([u'awk', u'{ print $2}'],
                                      env=self._env,
                                      stdin=p2.stdout, stdout=subprocess.PIPE)
                revision = p3.communicate()[0]
                svn_diff = u''

            else:

                cmd = [u'cp', u'-rf',
                       os.path.join(rootdir, branchname),
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

                if os.path.isdir(os.path.join(rootdir, u'.svn')):
                    try:
                        p_svn_diff = subprocess.Popen([u'svn', u'diff'],
                                                      env=self._env,
                                                      cwd=rootdir,
                                                      stdout=subprocess.PIPE)
                        svn_diff = p_svn_diff.communicate()[0]
                        svn_diff = svn_diff.decode(u'utf8', u'replace')
                    except:
                        svn_diff = u'error in  svn diff'
                    try:
                        p_svn_log = subprocess.Popen([u'svn', u'log',
                                                      u'-l', u'1'],
                                                     env=self._env,
                                                     cwd=rootdir,
                                                     stdout=subprocess.PIPE,
                                                     stderr=subprocess.STDOUT)
                        svn_log = p_svn_log.communicate()[0]
                        svn_log = svn_log.decode(u'utf8',
                                                 u'replace')
                        res = rev_search.search(svn_log)
                        if (res):
                            revision = int(res.group(1))
                    except:
                        revision = -1

                else:
                    try:
                        p = subprocess.Popen([u'git',u'diff',
                                              u'remotes/git-svn', u'.'],
                                             env=self._env,
                                             cwd=os.path.join(rootdir,
                                                              branchname),
                                             stdout=subprocess.PIPE)
                        svn_diff = p.communicate()[0]
                        svn_diff = svn_diff.decode(u'utf8', u'replace')
                    except:
                        svn_diff = u"Could not get diff"
                    try:
                        p = subprocess.Popen([u'git', u'svn',
                                              u'log', u'-n', u'1'],
                                             env=self._env,
                                             cwd=os.path.join(rootdir,
                                                              branchname),
                                             stdout=subprocess.PIPE,
                                             stderr=subprocess.STDOUT)
                        svn_log = p.communicate()[0]
                        svn_log = svn_log.decode(u'utf8', u'replace')
                        res = rev_search.search(svn_log)
                        if (res):
                            revision = int(res.group(1))
                    except:
                        revision = -1

            f=open(diffFile, u'w')
            f.write(svn_diff)
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

    def getText(nodelist):
        """ Return text contained in a xml nodelist"""
        rc = []
        for node in nodelist:
            if node.nodeType == node.TEXT_NODE:
                rc.append(unicode(node.data))
        return (''.join(rc)).replace("'","\\'")

    def _checkError(self):
        return
