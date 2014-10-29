#!/usr/bin/perl

###############################################################################
#
#   Module: Common.pm
#
###############################################################################
#
#   This file contains common datas used by PaStiX scripts.
#
#   It include Conf.pm which can contain user's configuration
#
#   Authors:
#     Mathieu Faverge - faverge@labri.fr
#     Xavier  Lacoste - lacoste@labri.fr
#
###############################################################################
package Modules::Common;
use strict;
use warnings;
use vars qw(@ISA @EXPORT);
use Exporter;

our @ISA    = qw(Exporter);
our @EXPORT = qw( 
    $RESULT_PATH
    $FILL_DB_PATH
    $PASTIX_ROOT
    $SCRIPT_PATH
    $TPL_ROOT
    $METIS_HOME
    $SCOTCH_HOME
    $NATURAL_DOCS_DIR
    $PASTIX_DOC
    $PASTIX_PROJECT 
    $PASTIX_NDSCRIPTS
    $USERMAIL
    %PASTIX_SRC
    $SITE_PATH
    @liste_runs
    %liste_cas
    %machines
    $db
    $host
    $userid
    $passwd
    $connectionInfo
    $MAX_IPARM
    $MAX_DPARM
    %metriqueconf
    %keywords
    $memkeyo
    $memkeyKo
    $memkeyMo
    $memkeyGo
    $errorkey
    $warningkey
    @OPT_KEY
    @execlist
    @optlist
    @negoptlist
    @posoptlist
    @paramlist
    @paramlist2
    @colors
    $machine
    $system
    $date
    %data
    PrintOption
    MakeTemplate
    );

###############################################################################
# Group: Default configuration 
#   
#  Configuration variables which can be redefined in Conf.pm


##
## Options Par defaut si conf.pm n'existe pas
##

#
# Strings: PATHS
#
#   PASTIX_ROOT      - PaStiX sources directory.
#   RESULT_PATH      - Root of the runs tree.
#   FILL_DB_PATH     - Directory used to fill database.
#   SCRIPT_PATH      - Directory containing PaStix scripts.
#   TPL_ROOT         - Directory containing templates.
#   METIS_HOME       - Directory containing metis librairies and headers.
#   SCOTCH_HOME      - Directory containing Scotch librairies and headers.
#   NATURAL_DOCS_DIR - Directory containing NaturalDocs executable.
#   PASTIX_DOC       - Output for NaturalDocs.
#   PASTIX_PROJECT   - Project directory for NaturalDocs.
#   PASTIX_NDSCRIPTS - PaStiX Scripts for NaturalDocs directory.
#   SITE_PATH        - Directory which will contain the website associated with the database.
#
our $PASTIX_ROOT      = $ENV{ HOME}."/Pastix";
our $RESULT_PATH      = $PASTIX_ROOT."/results";
our $FILL_DB_PATH     = $PASTIX_ROOT."/results";
our $SCRIPT_PATH      = $PASTIX_ROOT."/Scripts"; 
our $TPL_ROOT         = $SCRIPT_PATH."/templates";
our $METIS_HOME       = $ENV{ METIS_HOME};
our $SCOTCH_HOME      = $ENV{ SCOTCH_HOME};
our $NATURAL_DOCS_DIR = "~/bin/NaturalDocs";
our $PASTIX_DOC       = $PASTIX_ROOT."/Pastix-Doc";
our $PASTIX_PROJECT   = $PASTIX_ROOT."/Pastix-Doc-Project";
our $PASTIX_NDSCRIPTS = $SCRIPT_PATH."/NaturalDocs";
our $SITE_PATH        = $ENV{ HOME}."/pastix/Scripts/site";

#
# String: USERMAIL
#    Mail to send runs summary.
#
our $USERMAIL         = '';
#
# Hashtable: PASTIX_SRC
#   Root of each branches sources directory.
#
our %PASTIX_SRC       = (
    "T" => "$PASTIX_ROOT"."/trunk",
    "O" => "$PASTIX_ROOT"."/branches/ooc",
    "F" => "$PASTIX_ROOT"."/branches/FaxDist"
    );

#
# Array: liste_cas
#   list of hashtables describing parameters to the runs to execute.
#
our %liste_cas = (
    "classique" =>  
    {  
	'MIN_PROC'       => 1, 
	'MAX_PROC'       => 8, 
	'MIN_THREAD'     => 1, 
	'MAX_THREAD'     => 8,
	'MAX_PROCTOTAL'  => 32,
	'iparm'         => {
	    "IPARM_ABS"                 => 4,
	    "IPARM_MATRIX_VERIFICATION" => 0,
	    "IPARM_THREAD_COMM_MODE"    => 0,
	    "IPARM_NB_THREAD_COMM"      => 1,
	    "IPARM_INCOMPLETE"          => 0,
	    "IPARM_OOC_LIMIT"           => 50,
	}
    },
    "incomplet_loa5_lof1" =>  
    {  
	'MIN_PROC'       => 1, 
	'MAX_PROC'       => 8, 
	'MIN_THREAD'     => 1, 
	'MAX_THREAD'     => 8,
	'MAX_PROCTOTAL'  => 32,
	'iparm'         => {
	    "IPARM_ABS"                 => 4,
	    "IPARM_MATRIX_VERIFICATION" => 0,
	    "IPARM_THREAD_COMM_MODE"    => 0,
	    "IPARM_NB_THREAD_COMM"      => 1,
	    "IPARM_INCOMPLETE"          => 1,
	    "IPARM_OOC_LIMIT"           => 50,
	    "IPARM_AMALGAMATION_LEVEL"  => 5,
	    "IPARM_LEVEL_OF_FILL"       => 1
	}
    }
    );

#
# Array: liste_runs
#   list of hashtables describing runs to execute.
#
our @liste_runs = ( 
    { 'libnames'    => ["T_Numa_Recv"],
      'executables' => ["example/src/simple"],
      'parametres'  => ["classique"],
      'matrices'    => ["$PASTIX_SRC{'T'}/matrix/small.rsa" ],
    }
    );

#
# Hashtable: machines
#   Description of all machines
#
our %machines = (
    'default' => {
	'nbcores'    => 1,
	'mempernode' => 0,
	'execcmd'    => '',
	'template'   => '',
	'script'     => '',
	'submit'     => '',
	'time'       => "",
        'args'       => "",
        'argsbe'     => "",
        'bits'       => "_32bits",
	'hostarch'   => "i686_pc_linux",
        'ccprog'     => "gcc -Wall",
        'f77prog'    => "gfortran",
	'mpccprog'   => "mpicc -Wall",
	'mcfprog'    => "mpif90",
	'ccfopt'     => "-O3",
	'ccfdeb'     => "-O0 -g3",
        'f90flags'   => "-ffree-form -x f95-cpp-input",
        'ldflags'    => "-lblas -lgfortran -lm -lrt",
        'ccflags'    => "",
        'lkfopt'     => "-s",
        'arprog'     => "ar",
	'arflags'    => "-ruv"}
    );

# M3PEC
require "Modules/Machines/decrypthon.pm";
# PLAFRIM
require "Modules/Machines/plafrim.pm";
# G5000
require "Modules/Machines/gommette.pm";
require "Modules/Machines/zed.pm";
require "Modules/Machines/bordereau.pm";
require "Modules/Machines/chinqchint.pm";
# CCRT
require "Modules/Machines/titane.pm";
require "Modules/Machines/platine.pm";
require "Modules/Machines/argent.pm";
# CINES
require "Modules/Machines/r41.pm";
require "Modules/Machines/jade.pm";
# CREMI
require "Modules/Machines/hagrid.pm";
# Machines Feux
require "Modules/Machines/vulcain.pm";
require "Modules/Machines/hephaistos.pm";
require "Modules/Machines/loki.pm";
require "Modules/Machines/agni.pm";
# IDRIS
require "Modules/Machines/vargas.pm";
# Daltons
require "Modules/Machines/bertha.pm";


#
# Strings: Database information
#   
#   db             - Database name.
#   host	   - Host of the database.
#   userid 	   - User for the database.
#   passwd	   - Password for the database.
#   connectionInfo - Information to connect the database.
#
our $db="";
our $host="";
our $userid="";
our $passwd="";
our $connectionInfo="dbi:mysql:$db;$host";

#
# Option perso
#
require Modules::Conf or print "Create Modules/Conf.pm to configure scripts";

###############################################################################
# Group: Common options
#
#   Common options for PaStiX Scripts.
#

# 
# Integers: Maximum indexes
#
#   MAX_IPARM - Maximum index for integer parameters.
#   MAX_DPARM - Maximum index for double parameters.
#
our $MAX_IPARM = 64;
our $MAX_DPARM = 64;
#
# Arrays: List of keywords
#
#   optlist    - Compilation option keyword list
#   negoptlist -
#   posoptlist -
#   paramlist  - 
#   execlist   -
#
our @optlist = ('smp', 'mpi', 'stat', 'napa', 'thco', 'ircv', 'isnd', 
		'tag', 'fco', 'rfob', 'scotch',  'met',  'int',  'flt',
                'vers'); 

our @paramlist2 = ( 'size', 'nnza', 'nnzl', 'opc', 'esp', 'tord', 'tana', 'fillin', 'pfact', 'smsize', 'csize', 'ntsk');

our @negoptlist = ('smp', 'mpi', );
our @posoptlist = ('stat', 'napa', 'thco', 'ircv', 'isnd', 
		   'fco', 'rfob', 'scotch',  'met');

our @paramlist = ('mach',  'matr', 'exec', 'vers', 'day', 'hour', 
		  'fact', 'solv', 'refin', 'iter', 'norm', 
		  'nproc', 'nthrd', 'nbtask', 'mem', 'err', 'warn', 'sym' );

our @execlist  = ('day', 'hour', 'nproc', 'nthrd', 'nbtask', 'mach', 'matr', 'fact', 
		  'solv', 'refin', 'mem', 'mend', 'ovh', 'ovhm', 'ovhl', 'iter',
		  'err', 'warn');

#
# Hashtable: metriqueconf
#    Description of all metriques.
#
our %metriqueconf = (
    'date'  => { 'key' => '([0-9]{8}-[0-9]{4})',
		 'nom' => 'Date',
		 'fmt' => '%-13s',
                 'db'  => "DATE"},
    'hour'  => { 'key' => '[0-9]{8}-([0-9]{2})([0-9]{2})',
		 'nom' => 'Hour',
		 'fmt' => '%-4s',
                 'db'  => 'HEURE',
                 'desc'=> "Execution hour" },
    'day'   => { 'key' => '([0-9]{4})([0-9]{2})([0-9]{2})-[0-9]{4}',
		 'nom' => 'Day',
		 'fmt' => '%-10s',
                 'db'  => 'DATE',
                 'desc'=> "Execution date" },
    'exec'  => { 'key' => '',
		 'nom' => 'Executable',
 		 'fmt' => "%-40s",
                 'db'  => 'NOMPROGRAMME',
                 'desc'=> "Executable name" },
    'exemple' => { 'key' => '.*\/[0-9]{8}-[0-9]{4}\/(.*)\.log',
		   'nom' => 'Exemple',
		   'fmt' => '%-13s',
		   'db'  => "NOMEXEMPLE"},
    'nproc' => { 'key' => "",
		 'nom' => 'NbProc',
		 'fmt' => "%-6s",
		 'db'  => 'NBPROC',
                 'desc'=> "Processor number" },
    'nthrd' => { 'key' => "",
		 'nom' => 'NbThrd',
		 'fmt' => "%-6s",
                 'db'  => "NBTHRD",
                 'desc'=> "Thread number" },
    'nbtask'=> { 'key' => "",
		 'nom' => '',
		 'fmt' => "%-6s",
		 'db'  => "NBTOTALTASK",
                 'desc'=> "Total number of thread (Nbproc*Nbthread)"  },
    'cas'   => { 'key' => "",
		 'nom' => 'Cas',
		 'fmt' => "%-23s",
		 'db'  => ""},
    'mach'  => { 'key' => "",
		 'nom' => 'Machine',
		 'fmt' => "%-10s",
		 'db'  => "MACHINE",
                 'desc'=> "Cluster name" },
    'matr'  => { 'key' => "",
		 'nom' => 'Matrice',
		 'fmt' => "%-15s",
		 'db'  => "MATRICE",
                 'desc'=> "Matrice name" },
    'solv'  => { 'key' => "Time to solve[ ]*([-.e+0-9]*) s",
		 'nom' => 'Solve',
		 'fmt' => "%-8s",
                 'db'  => "TEMPS_SOLV",
                 'desc'=> "Solve time" },
    'fact'  => { 'key' => "Time to factorize[ ]*([-.e+0-9]*) s",
		 'nom' => 'Facto',
		 'fmt' => "%-8s" ,
		 'db'  => "TEMPS_FACTO",
                 'desc'=> "Numerical factorization time" },
    'refin' => { 'key' => "Time for refinement[ ]*([-.e+0-9]*) s",
		 'nom' => 'Raff',
		 'fmt' => "%-8s",
		 'db'  => "TEMPS_RAFF",
                 'desc'=> "Reffinment time" },
    'mem'   => { 'key' => "Max memory used after clean[ ]*([.e+0-9]* .o)",
		 'nom' => 'MaxMem',
		 'fmt' => "%-8s",
                 'db'  => "MEMOIRE",
                 'desc'=> "Maximum memory used" },
    'mend'  => { 'key' => "Memory used after clean[ ]*([.e+0-9]* [KMG]*o)",
		 'nom' => 'MenEnd',
		 'fmt' => "%-8s" ,
                 'db'  => "MEMOIRE_FIN",
                 'desc'=> "Memory used at the end" },
    'ovh'   => { 'key' => "Total.*overhead : ([.e+0-9]*)",
		 'nom' => 'OVH',
		 'fmt' => "%-8s",
                 'db'  => "OVERHEAD_TOTAL",
                 'desc'=> "Overhead" },
    'ovhm'  => { 'key' => "Maximum.*overhead : ([.e+0-9]*)",
		 'nom' => 'OVHm',
		 'fmt' => "%-8s",
                 'db'  => "OVERHEAD_MAX",
                 'desc'=> "Overhead maximum" },
    'ovhl'  => { 'key' => "Local.*overhead : ([.e+0-9]*)",
		 'nom' => 'OVHl',
		 'fmt' => "%-8s",
                 'db'  => "OVERHEAD_LOCAL",
                 'desc'=> "Overhead local" },
    'fr2'   => { 'key' => "fillrate2 ([.e+0-9]*)",
		 'nom' => 'FillRate',
		 'fmt' => "%-8s",
                 'db'  => "",
                 'desc'=> "Fill rate" },
    'iter'  => { 'key' => "Refinement[ ]*([0-9]+) iterations",
		 'nom' => 'NbIter',
		 'fmt' => "%-8s",
                 'db'  => "NBITER",
                 'desc'=> "Number of reffinment iteration" },
    'norm'  => { 'key' => ", norm=([0-9e\.-]+)",
		 'key2'=> "Precision : .* = ([0-9e\.-]+)",
		 'nom' => 'Norm',
		 'fmt' => "%-12s",
                 'db'  => "PRECISION",
                 'desc'=> "Precision"  },
    'vers'  => { 'key' => "Version[ ]*:[ ](.*)\n",
		 'nom' => 'Version',
		 'fmt' => "%-8s",
                 'db'  => "VERSION",
                 'desc'=> "SVN Revision number"  },
    'err'   => { 'key' => "^ERROR:",
		 'nom' => 'Error',
		 'fmt' => "%-8s",
                 'db'  => "ERROR",
                 'desc'=> "Number of errors"  }, 
    'warn'  => { 'key' => "^WARNING:",
		 'nom' => 'warning',
		 'fmt' => "%-8s",
                 'db'  => "WARNING",
                 'desc'=> "Number of warning"  }, 
    'smp'   => { 'key' => "  SMP_SOPALIN         : ([a-zA-Z ]*)",
		 'k2'  => "DFORCE_NOSMP",
		 'nom' => 'smp_sopalin',
		 'fmt' => "%-20s",
                 'db'  => "VERSION_SMP" }, 
    'mpi'   => { 'key' => "  VERSION MPI         : ([a-zA-Z ]*)",
		 'k2'  => "DFORCE_NOMPI",
		 'nom' => 'version_mpi',
		 'fmt' => "%-20s",
                 'db'  => "VERSION_MPI" }, 
    'stat'  => { 'key' => "  STATS_SOPALIN       : ([a-zA-Z ]*)",
		 'k2'  => "DSTATS_SOPALIN",
		 'nom' => 'stats_sopalin',
		 'fmt' => "%-20s",
                 'db'  => "STATS_SOPALIN" }, 
    'napa'  => { 'key' => "  NAPA_SOPALIN        : ([a-zA-Z ]*)",
		 'k2'  => 'DNAPA_SOPALIN',
		 'nom' => 'napa_sopalin',
		 'fmt' => "%-20s",
                 'db'  => "NAPA_SOPALIN" }, 
    'thco'  => { 'key' => "  THREAD_COMM         : ([a-zA-Z ]*)",
		 'k2'  => "DTHREAD_COMM",
		 'nom' => 'thread_comm',
		 'fmt' => "%-20s",
                 'db'  => "THREAD_COMM" }, 
    'ircv'  => { 'key' => "  TEST_IRECV          : ([a-zA-Z ]*)",
		 'k2'  => 'DTEST_IRECV',
		 'nom' => 'test_irecv',
		 'fmt' => "%-20s",
                 'db'  => "TEST_IRECV" }, 
    'isnd'  => { 'key' => "  TEST_ISEND          : ([a-zA-Z ]*)",
		 'k2'  => 'DTEST_ISEND',
		 'nom' => 'test_isend',
		 'fmt' => "%-20s",
                 'db'  => "TEST_ISEND" }, 
    'tag'   => { 'key' => "  TAG                 : ([a-zA-Z ]*)",
		 'k2'  => 'DEXACT_T',	
		 'nom' => 'tag',
		 'fmt' => "%-20s",
                 'db'  => "TAG" }, 
    'fco'   => { 'key' => "  FORCE_CONSO         : ([a-zA-Z ]*)",
		 'k2'  => 'DFORCE_CONSO',
		 'nom' => 'forceconso',
		 'fmt' => "%-20s",
                 'db'  => "FORCE_CONSO" }, 
    'rfob'  => { 'key' => "  RECV_FANIN_OR_BLOCK : ([a-zA-Z ]*)",
		 'k2'  => 'DRECV_FANNIN_OR_BLOCK',
		 'nom' => 'recv_fanin_or_block',
		 'fmt' => "%-20s",
                 'db'  => "RECV_FANIN_OR_BLOCK" }, 
    'scotch'=> { 'key' => "  WITH_SCOTCH            : ([a-zA-Z ]*)",
		 'k2'  => 'DWITH_SCOTCH',
		 'nom' => 'scotch',
		 'fmt' => "%-20s",
                 'db'  => "SCOTCH" }, 
    'met'   => { 'key' => "  METIS               : ([a-zA-Z ]*)",
		 'k2'  => 'DMETIS',
		 'nom' => 'metis',
		 'fmt' => "%-20s",
                 'db'  => "METIS" }, 
    'int'   => { 'key' => "  INTEGER TYPE        : ([a-zA-Z0-9_ ]*)",
		 'k2'  => '',
		 'nom' => 'integer_type',
		 'fmt' => "%-20s",
                 'db'  => "VERSION_INT" }, 
    'flt'   => { 'key' => "  FLOAT TYPE          : ([a-zA-Z ]*)",
		 'k2'  => '',
		 'nom' => 'float_type',
		 'fmt' => "%-20s",
                 'db'  => "VERSION_FLOAT" },
    'size'  => { 'key' => "  Matrix size  [ ]*([0-9]*) x [0-9]*",
		 'k2'  => '',
		 'nom' => 'N',
		 'fmt' => "%-9s",
                 'db'  => "",
                 'desc'=> "Matrix size N" },

    'nnza'  => { 'key' => '[0-9]* : Number of nonzeros \(local block structure\)[ ]*([0-9]*)',
		 'k2'  => '',
		 'nom' => 'NNZ_A',
		 'fmt' => "%-15s",
                 'db'  => "",
                 'desc'=> "Matrix number of nonzeros (Only lower part in symmetric graph)" },
    'tord'  => { 'key' => "Time to compute ordering [ ]*([-.e+0-9]*) s",
		 'k2'  => '',
		 'nom' => 'Ord',
		 'fmt' => "%-6s",
                 'db'  => "",
                 'desc'=> "Ordering Time" },
    'tana'  => { 'key' => "Time to analyze  [ ]*([-.e+0-9]*) s",
		 'k2'  => '',
		 'nom' => 'Blend',
		 'fmt' => "%-6s",
                 'db'  => "",
                 'desc'=> "Analyse Time" },
    'nnzl'  => { 'key' => 'Number of nonzeros in factorized matrice[ ]*([0-9]*)',
		 'k2'  => '',
		 'nom' => 'NNZ_L',
		 'fmt' => "%-15s",
                 'db'  => "",
                 'desc'=> "Number of nonzeros of L" },
    'nnzlbl'=> { 'key' => 'Number of nonzeros (block structure)[ ]*([0-9]*)',
		 'k2'  => '',
		 'nom' => 'NNZ_L_BLOCK',
		 'fmt' => "%-15s",
                 'db'  => "",
                 'desc'=> "Number of nonzeros of L (block structure)" },
    'opc'   => { 'key' => 'Number of operations \(L[ULt]*\) [ ]*([-.e+0-9]*)',
		 'k2'  => '',
		 'nom' => 'OPC',
		 'fmt' => "%-15s",
                 'db'  => "",
                 'desc'=> "Number of OPC (LLt or LU)" },
    'sym'    => {'key' => 'Number of operations \((L[ULt]*)\) [ ]*[-.e+0-9]*',
		 'k2'  => '',
		 'nom' => 'Symmetric',
		 'fmt' => "%-1s",
                 'db'  => "",
                 'desc'=> "Matrix symmetry" },
    'fillin' => {'key' => "Fill-in [ ]*([-.e+0-9]*)",
		 'k2'  => '',
		 'nom' => 'Fill-in',
		 'fmt' => "%-8s",
                 'db'  => "",
                 'desc'=> "Fill-in" },
    'fillinb'=> {'key' => "Fill-in (block structure)[ ]*([-.e+0-9]*)",
		 'k2'  => '',
		 'nom' => 'Fill-in-block',
		 'fmt' => "%-8s",
                 'db'  => "",
                 'desc'=> "Fill-in (block structure)" },
    'pfact' => { 'key' => 'Prediction Time to factorize \(IBM PWR5 ESSL\)[ ]*([-.e+0-9]*) s',
		 'k2'  => '',
		 'nom' => 'PFact',
		 'fmt' => "%-8s",
                 'db'  => "",
                 'desc'=> "Prediction time for the numerical fatorization" },
    'smsize' => { 'key' => 'SolverMatrix size \(without coefficient\)[ ]*([-.e+0-9]*) [KMG]*o',
		 'k2'  => '',
		 'nom' => 'SMsize',
		 'fmt' => "%-8s",
                 'db'  => "",
                 'desc'=> "Size of the structure SolverMatrix" },    
    'csize' => { 'key' => 'Maximum coeftab size \(cefficients\)[ ]*([-.e+0-9]* [KMG]*o)',
		 'k2'  => '',
		 'nom' => 'CoefSize',
		 'fmt' => "%-8s",
                 'db'  => "",
                 'desc'=> "Size of coeftab" },    
    'esp'   => { 'key' => "Number of tasks added by esp[ ]*([-0-9]*)",
		 'k2'  => '',
		 'nom' => 'ESP',
		 'fmt' => "%-8s",
                 'db'  => "",
                 'desc'=> "Number of tasks added by the esp option" },
    'ntsk'  => { 'key' => "Number of tasks  [ ]*([-0-9]*)",
		 'k2'  => '',
		 'nom' => 'NTask',
		 'fmt' => "%-8s",
                 'db'  => "",
                 'desc'=> "Number of tasks" },
   );

#
# Hashtable: keywords
#   Keywords to parse logs.
#
our %keywords = (
    'defined'  => { 'Not defined'     => 0,
		    'Defined'         => 1},
    'float'    => { 'double complex'  => 3,
		    'simple complex'  => 2,
		    'double'          => 1,
		    'simple'          => 0},
    'integer'  => { 'int64_t'         => 3,
		    'int32_t'         => 2,
		    'long'            => 1,
		    'int'             => 0} ,
    'tag'      => { 'Exact Thread'    => 0,
                    'Exact Tag'       => 1,
                    'Tag Fanin or block' => 2}
    );
#
# Strings: Parse keys
#
#   Keys to parse logs
#
#   errorkey   - Key to catch errors.
#   warningkey - Key to catch warnings.
#   memkeyo    - Key to catch memory usage in octets.
#   memkeyKo   - Key to catch memory usage in kilo-octets.
#   memkeyMo   - Key to catch memory usage in mega-octets.
#   memkeyGo   - Key to catch memory usage in giga-octets.
#   OPT_KEY    - Keys to catch compilation options.
#
our $errorkey    = "^ERROR:";
our $warningkey  = "^WARNING:";
our $memkeyo     = "Max memory used after clean[ ]*([.e+0-9]*) o";
our $memkeyKo    = "Max memory used after clean[ ]*([.e+0-9]*) Ko";
our $memkeyMo    = "Max memory used after clean[ ]*([.e+0-9]*) Mo";
our $memkeyGo    = "Max memory used after clean[ ]*([.e+0-9]*) Go";

our @OPT_KEY;
$OPT_KEY[0]  = "  SMP_SOPALIN         : ([a-zA-Z ]*)";
$OPT_KEY[1]  = "  VERSION MPI         : ([a-zA-Z ]*)";
$OPT_KEY[2]  = "  STATS_SOPALIN       : ([a-zA-Z ]*)";
$OPT_KEY[3]  = "  NAPA_SOPALIN        : ([a-zA-Z ]*)";
$OPT_KEY[4]  = "  THREAD_COMM         : ([a-zA-Z ]*)";
$OPT_KEY[5]  = "  TEST_IRECV          : ([a-zA-Z ]*)";
$OPT_KEY[6]  = "  TEST_ISEND          : ([a-zA-Z ]*)";
$OPT_KEY[7]  = "  TAG                 : ([a-zA-Z ]*)";
$OPT_KEY[8]  = "  FORCE_CONSO         : ([a-zA-Z ]*)";
$OPT_KEY[9]  = "  RECV_FANIN_OR_BLOCK : ([a-zA-Z ]*)";
$OPT_KEY[10] = "  FLUIDBOX            : ([a-zA-Z ]*)";
$OPT_KEY[11] = "  METIS               : ([a-zA-Z ]*)";
$OPT_KEY[12] = "  INTEGER TYPE        : ([a-zA-Z ]*)";
$OPT_KEY[13] = "  FLOAT TYPE          : ([a-zA-Z ]*)";

#
# Array: colors
#   List of colors
#
our @colors = ( "black", "red", "green", "blue", "cyan", 
	       "magenta", "yellow" , "gris25", "rose", 
	       "violet", "bordeaux", "turquoise", "ciel",
	       "mer","orange","vertfonce","violetpastel");

#
# Strings: Running information.
#
#   Few information about the machine and 
#   the running environnement.
#
#   machine - Machine name.
#   system  - Operating system name.
#   date    - Date when the Script are launched.
#

# Definition du nom de la machine
our $machine;
$machine = `uname -n`; chop $machine;
$machine =~ s/([a-zA-Z]*)[-0-9]*/$1/ if (!($machine =~ /r41/ )) ; 
$machine =~ s/\.[a-z]*\.grid5000\.fr//;; # Suppression Grid5000
$machine =~ s/\.dalton//;;       # Suppression dalton
$machine = "argent" if ($machine =~ /platine/); 

# Definition du system d'exploitation.
our $system  = `uname`;

# Définition de la date de lancement
our $date = `date +%Y%m%d-%H%M`; chop $date;

#
# Hashtable: data
#   Hashtable used to store metrics values.
#

# Données pour la récupération des temps de calculs
our %data = ();

# our %nbproc_by_machine = ( 'decrypthon' => 16,
# 			  'eowin'      => 2,
# 			  'hagrid'     => 16 );


###############################################################################
# Group: Functions

# Function: PrintOption
#
#   Print few options
#
sub PrintOption {
    my $i;
    
    print "\n";
    print "Répertoire de stockage des résultats : $RESULT_PATH\n";
    print "Répertoire de stockage de PaStiX     : $PASTIX_ROOT\n";
    print "Répertoire des templates             : $TPL_ROOT\n";
    print "\n";
}

# Function: MakeTemplate
#
# Creates a file from a template.
#   Replaces all keys found by the associated values.
#
# Parameters: 
#   tpl  - Template file.
#   file - Output file.
#   vars - Couples (key, value) to create the *file*.
#
sub MakeTemplate {
    
    my @inp;
    my ($tpl, $file, %vars) = @_ ;

    ###### Creation fichier #####
    open(CMD, "<$tpl") or die("File doesn't exist");
    @inp = <CMD>;
    close CMD;
 
    foreach (@inp)
    {
	foreach my $key (keys %vars)
	{
	    s/$key/$vars{$key}/g;
	}
    }
    
    open(FOUT, ">$file");
    print FOUT @inp;
    close FOUT;
}

1;
