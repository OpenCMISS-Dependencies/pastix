'''
Configuration file example 
To use it rename it conf.py
'''

import os
import sys
from copy import deepcopy
sys.path.append(u'..')
from Projects.PaStiX import PaStiX
from Projects.Hips import Hips
from Projects.Scotch import Scotch

### Database parameters 
_database_data = {
    u"host"   : u"localhost",
    u"user"   : u"pastix",
    u"passwd" : u"PaStiX",
    u"db"     : u"pastix-python"}

## Used when no -export is set
_pastix_dir = u"/path/to/ricar/directory  or None" 
_hips_dir   = u"/path/to/hips/directory   or None"
_scotch_dir = u"/path/to/scotch/directory or None"

_projects = {
    u'PaStiX' : { 
        u'func'        : PaStiX,
        u'dir'         : _pastix_dir,
        u'svnurl'      : u'git+ssh://%s@scm.gforge.inria.fr//gitroot/ricar/ricar.git',
        u"export_name" : u'/tmp/%s/ricar.tgz',
        u"export_dir"  : u'/tmp/%s/ricar'},
    u'Hips'   : { 
        u'func'        : Hips,
        u'dir'         : _hips_dir,
        u'svnurl'      : u'svn+ssh://%s@scm.gforge.inria.fr/svn/hips',
        u"export_name" : u'/tmp/%s/hips.tgz',
        u"export_dir"  : u'/tmp/%s/hips'},
    u'Scotch'   : { 
        u'func'        : Scotch,
        u'dir'         : _scotch_dir,
        u'svnurl'      : u'svn+ssh://%s@scm.gforge.inria.fr/svn/scotch',
        u"export_name" : u'/tmp/%s/scotch.tgz',
        u"export_dir"  : u'/tmp/%s/scotch'}
}
### Definitions des parametres d'execution
cas_classique = {  
    'NAME'       : u"cas_classique",
    'PROC'       : [1,2,4,8], 
    'THREAD'     : [1,2,4,8],
    'iparms'          : {
        "IPARM_ABS"                 : 4,
        "IPARM_MATRIX_VERIFICATION" : 0,
        "IPARM_THREAD_COMM_MODE"    : 0,
        "IPARM_NB_THREAD_COMM"      : 1,
        "IPARM_INCOMPLETE"          : 0,
        "IPARM_OOC_LIMIT"           : 50 }
    }
cas_incomplet_loa5_lof1 =  {  
    'NAME'       : u"cas_incomplet_loa5_lof1",
    'PROC'       : [1,2,4,8], 
    'THREAD'     : [1,2,4,8],
    'iparms'          : {
        "IPARM_ABS"                 : 4,
        "IPARM_MATRIX_VERIFICATION" : 0,
        "IPARM_THREAD_COMM_MODE"    : 0,
        "IPARM_NB_THREAD_COMM"      : 1,
        "IPARM_INCOMPLETE"          : 1,
        "IPARM_OOC_LIMIT"           : 50,
        "IPARM_AMALGAMATION_LEVEL"  : 5,
        "IPARM_LEVEL_OF_FILL"       : 1
        }
    }

### Definition des sets d'options de compilation
base_run = {
    u'project' : u'project_name',
    u'machines':[],
    u'remote_machines' : [],
    u'compilations':[],
    u'binaries':[],
    u'parameters':[],
    u'cases':[]
}


### VERSIONS CENTRALISEES
liste_builds = [u'T_int64',                                 
                u"T_int64_nompi_nosmp",                     
                u"T_int64_nompi",                           
                u"T_int64_nompi_PASTIX_DYNSCHED",           
                u"T_int64_nosmp",                           
                u"T_int64_nosmp_TEST_IRECV",                
                                                             
                u"T_NO_NUMA_ALLOC_int64",                   
                u"T_int64_TEST_IRECV",                      
                u"T_int64_PASTIX_DYNSCHED",     
                u"T_int64_TRACE_SOPALIN",                   
                u"T_int64_MULT_SMX",

                u'T_int32',                                 
                u"T_int32_nompi_nosmp",                     
                u"T_int32_nompi",                           
                u"T_int32_nosmp",                           
                u"T_int32_nosmp_TEST_IRECV",                
                u"T_int32_TEST_IRECV",                      
                u"T_int64_ptscotch_nosmp",           # Full-MPI Recv
                u"T_int64_ptscotch_nosmp_TEST_IRECV",     # Full-MPI IRecv
                u"T_int64_ptscotch",                 # MPI-SMP
                u"T_int64_ptscotch_TEST_IRECV",           # En IRecv
                u"T_int64_ptscotch_PASTIX_DYNSCHED",   # avec 1 thread de comm
                u"T_int64_ptscotch_TRACE_SOPALIN",           # trace

                u"T_int32_ptscotch_nosmp",           # Full-MPI Recv
                u"T_int32_ptscotch_nosmp_TEST_IRECV",     # Full-MPI IRecv
                u"T_int32_ptscotch",                 # MPI-SMP
                u"T_int32_ptscotch_TEST_IRECV"           # En IRecv
                ]

#not available on plafrim yet
liste_builds_metis = []
#liste_builds_metis = [u"T_int_metis"]

### PLAFRIM REEL
p = PaStiX()
_liste_run = []
run = deepcopy(base_run)
run[u'project'] = u'PaStiX'
run[u'machines'].append(u"plafrim-gnu-openmpi-mkl-ramet")
run[u'remote_machines'].append(u'devel09.plafrim.cluster-ramet')
for name in liste_builds:
    run[u'compilations'].append(p.buildCompilationFromName(name, u"master"))
run[u'binaries'].append(u"simple")
run[u'binaries'].append(u"ssimple")
run[u'binaries'].append(u"dsimple")
run[u'parameters'].append(cas_classique)
run[u'cases'].append(u"/home/ramet/pastix/trunk/matrix/orsirr.rua")
run[u'cases'].append(u"/home/ramet/pastix/trunk/matrix/small.rsa")
_liste_run.append(run)        


#### PLAFRIM COMPLEXE
p = PaStiX()
run = deepcopy(base_run)
run[u'project'] = u'PaStiX'
run[u'machines'].append(u"plafrim-gnu-openmpi-mkl-ramet")
run[u'remote_machines'].append(u'devel09.plafrim.cluster-ramet')
for name in liste_builds:
    run[u'compilations'].append(p.buildCompilationFromName(name, u"master"))
run[u'binaries'].append(u"csimple")
run[u'binaries'].append(u"zsimple")
run[u'parameters'].append(cas_classique)
run[u'cases'].append(u"/home/ramet/pastix/trunk/matrix/young4c.mtx")
_liste_run.append(run)    

### VULCAIN REEL

run = deepcopy(base_run)
run[u'project'] = u'PaStiX'
run[u'machines'].append(u"vulcain-gnu-mpich2-mkl-pastix-user")
#run[u'machines'].append(u"vulcain-intel-mpich2-mkl")
run[u'remote_machines'].append(u'vulcain.bordeaux.inria.fr-pastix-user')
for name in liste_builds:
    run[u'compilations'].append(p.buildCompilationFromName(name, u"master"))
for name in liste_builds_metis:
    run[u'compilations'].append(p.buildCompilationFromName(name, u"master"))
run[u'binaries'].append(u"simple")
run[u'binaries'].append(u"ssimple")
run[u'binaries'].append(u"dsimple")
run[u'parameters'].append(cas_classique)
run[u'cases'].append(u"/matrices/orsirr.rua")
run[u'cases'].append(u"/matrices/small.rsa")
_liste_run.append(run) 


### VULCAIN COMPLEXE

run = deepcopy(base_run)
run[u'project'] = u'PaStiX'
run[u'machines'].append(u"vulcain-gnu-mpich2-mkl-pastix-user")
#run[u'machines'].append(u"vulcain-intel-mpich2-mkl")
run[u'remote_machines'].append(u'vulcain.bordeaux.inria.fr-pastix-user')
for name in liste_builds:
    run[u'compilations'].append(p.buildCompilationFromName(name, u"master"))
for name in liste_builds_metis:
    run[u'compilations'].append(p.buildCompilationFromName(name, u"master"))
run[u'binaries'].append(u"csimple")
run[u'binaries'].append(u"zsimple")
run[u'parameters'].append(cas_classique)
run[u'cases'].append(u"/matrices/young4c_.mtx")
_liste_run.append(run) 

