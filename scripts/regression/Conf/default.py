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

### Database parameters 
_database_data = {
    u"host"   : u"localhost",
    u"user"   : u"pastix",
    u"passwd" : u"PaStiX",
    u"db"     : u"pastix-python"}
    
_pastix_dir = None
_hips_dir   = None
    
#### Definitions des parametres d'execution
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
                u"T_int64_nosmp_THREAD_COMM",               
                
                u"T_NO_NUMA_ALLOC_int64",                   
                u"T_int64_TEST_IRECV",                      
                u"T_int64_THREAD_COMM",                     
                u"T_int64_PASTIX_FUNNELED",                 
                u"T_int64_PASTIX_DYNSCHED_THREAD_COMM",     
                u"T_int64_PASTIX_DYNSCHED_PASTIX_FUNNELED", 
                u"T_int64_TRACE_SOPALIN",                   
                
                u'T_int32',                                 
                u"T_int_metis",                           
                u"T_int32_nompi_nosmp",                     
                u"T_int32_nompi",                           
                u"T_int32_nosmp",                           
                u"T_int32_nosmp_TEST_IRECV",                
                u"T_int32_nosmp_THREAD_COMM",               
                u"T_int32_TEST_IRECV",                      
                u"T_int32_THREAD_COMM",                     
                u"T_int32_PASTIX_FUNNELED"]
    
### PLAFRIM REEL
p = PaStiX()
_liste_run = []
run = deepcopy(base_run)
run[u'project'] = u'PaStiX'
run[u'machines'].append(u"plafrim-gnu-openmpi-mkl")
run[u'remote_machines'].append(u'devel09.plafrim.cluster-ramet')
for name in liste_builds:
    run[u'compilations'].append(p.buildCompilationFromName(name))
run[u'binaries'].append(u"simple")
run[u'binaries'].append(u"ssimple")
run[u'binaries'].append(u"dsimple")
run[u'parameters'].append(cas_classique)
run[u'cases'].append(u"/home/ramet/pastix/trunk/matrix/orsirr.rua")
run[u'cases'].append(u"/home/ramet/pastix/trunk/matrix/small.rsa")
_liste_run.append(run)        


#### PLAFRIM COMPLEXE
run = deepcopy(base_run)
run[u'project'] = u'PaStiX'
run[u'machines'].append(u"plafrim-gnu-openmpi-mkl")
run[u'remote_machines'].append(u'devel09.plafrim.cluster-ramet')
for name in liste_builds:
    run[u'compilations'].append(p.buildCompilationFromName(name))
run[u'binaries'].append(u"csimple")
run[u'binaries'].append(u"zsimple")
run[u'parameters'].append(cas_classique)
run[u'cases'].append(u"/home/ramet/pastix/trunk/matrix/young4c.mtx")
_liste_run.append(run)    

### VULCAIN REEL
run = deepcopy(base_run)
run[u'project'] = u'PaStiX'
run[u'machines'].append(u"vulcain-gnu-mpich2-mkl")
#run[u'machines'].append(u"vulcain-intel-mpich2-mkl")
run[u'remote_machines'].append(u'vulcain.bordeaux.inria.fr-pastix-user')
for name in liste_builds:
    run[u'compilations'].append(p.buildCompilationFromName(name))
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
run[u'machines'].append(u"vulcain-gnu-mpich2-mkl")
#run[u'machines'].append(u"vulcain-intel-mpich2-mkl")
run[u'remote_machines'].append(u'vulcain.bordeaux.inria.fr-pastix-user')
for name in liste_builds:
    run[u'compilations'].append(p.buildCompilationFromName(name))
run[u'binaries'].append(u"csimple")
run[u'binaries'].append(u"zsimple")
run[u'parameters'].append(cas_classique)
run[u'cases'].append(u"/matrices/young4c_.mtx")
_liste_run.append(run) 

##### PTSCOTCH

liste_builds = [u"T_int64_ptscotch_nosmp",                # Full-MPI Recv
                u"T_int64_ptscotch_nosmp_TEST_IRECV",     # Full-MPI IRecv
                u"T_int64_ptscotch_nosmp_THREAD_COMM",    # Full-MPI a 2 threads (calcul+Comm)
                u"T_int64_ptscotch",                      # MPI-SMP
                u"T_int64_ptscotch_TEST_IRECV",           # En IRecv
                u"T_int64_ptscotch_THREAD_COMM",          # avec 1 thread de comm
                u"T_int64_ptscotch_PASTIX_FUNNELED",      # Thread Funneled
                u"T_int64_ptscotch_DYNSCHED_THREAD_COMM", # avec 1 thread de comm
                u"T_int64_ptscotch_DYNSCHED_THREAD_PASTIX_FUNNELED", # Thread Funneled
                u"T_int64_ptscotch_TRACE_SOPALIN",        # trace
                u"T_int32_ptscotch_nosmp",                # Full-MPI Recv
                u"T_int32_ptscotch_nosmp_TEST_IRECV",     # Full-MPI IRecv
                u"T_int32_ptscotch_nosmp_THREAD_COMM",    # Full-MPI a 2 threads (calcul+Comm)
                u"T_int32_ptscotch",                      # MPI-SMP
                u"T_int32_ptscotch_TEST_IRECV",           # En IRecv
                u"T_int32_ptscotch_THREAD_COMM",          # avec 1 thread de comm
                u"T_int32_ptscotch_PASTIX_FUNNELED"]

### PLAFRIM REEL
run = deepcopy(base_run)
run[u'project'] = u'PaStiX'
run[u'machines'].append(u"plafrim-gnu-openmpi-mkl")
run[u'remote_machines'].append(u'devel09.plafrim.cluster-ramet')
for name in liste_builds:
    run[u'compilations'].append(p.buildCompilationFromName(name))
run[u'binaries'].append(u"simple")
run[u'binaries'].append(u"ssimple")
run[u'binaries'].append(u"dsimple")
run[u'parameters'].append(cas_classique)
run[u'cases'].append(u"/home/ramet/pastix/trunk/matrix/orsirr.rua")
run[u'cases'].append(u"/home/ramet/pastix/trunk/matrix/small.rsa")
_liste_run.append(run)        


#### PLAFRIM COMPLEXE
run = deepcopy(base_run)
run[u'project'] = u'PaStiX'
run[u'machines'].append(u"plafrim-gnu-openmpi-mkl")
run[u'remote_machines'].append(u'devel09.plafrim.cluster-ramet')
for name in liste_builds:
    run[u'compilations'].append(p.buildCompilationFromName(name))
run[u'binaries'].append(u"csimple")
run[u'binaries'].append(u"zsimple")
run[u'parameters'].append(cas_classique)
run[u'cases'].append(u"/home/ramet/pastix/trunk/matrix/young4c.mtx")
_liste_run.append(run)    

### VULCAIN REEL
run = deepcopy(base_run)
run[u'project'] = u'PaStiX'
run[u'machines'].append(u"vulcain-gnu-mpich2-mkl")
#run[u'machines'].append(u"vulcain-intel-mpich2-mkl")
run[u'remote_machines'].append(u'vulcain.bordeaux.inria.fr-pastix-user')
for name in liste_builds:
    run[u'compilations'].append(p.buildCompilationFromName(name))
run[u'binaries'].append(u"simple_dist")
run[u'binaries'].append(u"ssimple_dist")
run[u'binaries'].append(u"dsimple_dist")
run[u'parameters'].append(cas_classique)
run[u'cases'].append(u"/matrices/orsirr.rua")
run[u'cases'].append(u"/matrices/small.rsa")
_liste_run.append(run) 

### VULCAIN COMPLEXE
run = deepcopy(base_run)
run[u'project'] = u'PaStiX'
run[u'machines'].append(u"vulcain-gnu-mpich2-mkl")
#run[u'machines'].append(u"vulcain-intel-mpich2-mkl")
run[u'remote_machines'].append(u'vulcain.bordeaux.inria.fr-pastix-user')
for name in liste_builds:
    run[u'compilations'].append(p.buildCompilationFromName(name))
run[u'binaries'].append(u"csimple_dist")
run[u'binaries'].append(u"zsimple_dist")
run[u'parameters'].append(cas_classique)
run[u'cases'].append(u"/matrices/young4c_.mtx")
_liste_run.append(run) 
