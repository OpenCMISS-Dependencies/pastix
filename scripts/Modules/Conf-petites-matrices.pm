$RESULT_PATH      = $ENV{ HOME}."/results";
$FILL_DB_PATH     = $ENV{ HOME}."/results/pastix-user";
$PASTIX_ROOT      = $ENV{ HOME}."/ricar";
$SCRIPT_PATH      = $PASTIX_ROOT."/Scripts"; 
$TPL_ROOT         = $SCRIPT_PATH."/templates";
$METIS_HOME       = "/opt/metis";
$SCOTCH_HOME      = '/opt/scotch/int64';
$NATURAL_DOCS_DIR = "~/NaturalDocs";
$PASTIX_DOC       = $PASTIX_ROOT."/Pastix-Doc";
$PASTIX_PROJECT   = $PASTIX_ROOT."/Pastix-Doc-Project";
$PASTIX_NDSCRIPTS = $SCRIPT_PATH."/NaturalDocs";

$USERMAIL         = 'pastix-log@localhost';

%PASTIX_SRC       = (
    "T" => "$PASTIX_ROOT"."/trunk",
    "O" => "$PASTIX_ROOT"."/ooc",
    "F" => "$PASTIX_ROOT"."/FaxDist",
    );
$SITE_PATH        = $PASTIX_ROOT."/site2";

@liste_cas = (
    "classique" =>  {  
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
    "incomplet_loa5_lof1" =>  {  
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
@liste_runs = ( 
    # Version Pretraitement multi-sequentiel
    # Reel ; classique 
    {'libnames'      => ["T_Numa_mem_int64_nompi_nosmp",          # sequentielle
			 "T_Numa_mem_int64_nompi",                # Full-SMP
			 "T_Numa_mem_int64_nompi_dyn",         # Full-SMP
			 "T_Numa_mem_int64_nosmp",                # Full-MPI Recv
			 "T_Numa_mem_int64_nosmp_IRecv",          # Full-MPI IRecv
			 "T_Numa_mem_int64_nosmp_ThComm",         # Full-MPI a 2 threads (calcul+Comm)
			 "T_Numa_mem_int64",                      # MPI-SMP
			 "T_Numa_mem_int32_metis",                # MPI-SMP Metis
			 "T_mem_int64",                           # MPI-SMP
			 "T_Numa_mem_int64_IRecv",                # En IRecv
			 "T_Numa_mem_int64_ThComm",               # avec 1 thread de comm
			 "T_Numa_mem_int64_Funneled",             # Thread Funneled
			 "T_Numa_mem_int64_dyn_ThComm",        # avec 1 thread de comm
			 "T_Numa_mem_int64_dyn_Funneled",      # Thread Funneled
			 "T_Numa_mem_int64_trace",                # trace
			 "T_Numa_mem_int32_nompi_nosmp",          # sequentielle
			 "T_Numa_mem_int32_nompi",                # Full-SMP
			 "T_Numa_mem_int32_nosmp",                # Full-MPI Recv
			 "T_Numa_mem_int32_nosmp_IRecv",          # Full-MPI IRecv
			 "T_Numa_mem_int32_nosmp_ThComm",         # Full-MPI a 2 threads (calcul+Comm)
			 "T_Numa_mem_int32",                      # MPI-SMP
			 "T_Numa_mem_int32_IRecv",                # En IRecv
			 "T_Numa_mem_int32_ThComm",               # avec 1 thread de comm
			 "T_Numa_mem_int32_Funneled"],            # Thread Funneled
     
     'executables'   => ["example/bin/simple",
			 "example/bin/ssimple",
			 "example/bin/dsimple"],

     'parametres'    => ["classique"],
     'matrices'    => [$PASTIX_SRC{"T"}."/matrix/orsirr.rua",
                       $PASTIX_SRC{"T"}."/matrix/small.rsa"], 
 

    },    
    # Version Pretraitement multi-sequentiel
    # Complex ; classique
    {'libnames'      => ["T_Numa_mem_int64_nompi_nosmp",          # sequentielle
			 "T_Numa_mem_int64_nompi",                # Full-SMP
			 "T_Numa_mem_int64_nompi_dyn",         # Full-SMP
			 "T_Numa_mem_int64_nosmp",                # Full-MPI Recv
			 "T_Numa_mem_int64_nosmp_IRecv",          # Full-MPI IRecv
			 "T_Numa_mem_int64_nosmp_ThComm",         # Full-MPI a 2 threads (calcul+Comm)
			 "T_Numa_mem_int64",                      # MPI-SMP
			 "T_mem_int64",                           # MPI-SMP
			 "T_Numa_mem_int64_IRecv",                # En IRecv
			 "T_Numa_mem_int64_ThComm",               # avec 1 thread de comm
			 "T_Numa_mem_int64_Funneled",             # Thread Funneled
			 "T_Numa_mem_int64_dyn_ThComm",        # avec 1 thread de comm
			 "T_Numa_mem_int64_dyn_Funneled",      # Thread Funneled
			 "T_Numa_mem_int64_trace",                # trace
			 "T_Numa_mem_int32_nompi_nosmp",          # sequentielle
			 "T_Numa_mem_int32_nompi",                # Full-SMP
			 "T_Numa_mem_int32_nosmp",                # Full-MPI Recv
			 "T_Numa_mem_int32_nosmp_IRecv",          # Full-MPI IRecv
			 "T_Numa_mem_int32_nosmp_ThComm",         # Full-MPI a 2 threads (calcul+Comm)
			 "T_Numa_mem_int32",                      # MPI-SMP
			 "T_Numa_mem_int32_IRecv",                # En IRecv
			 "T_Numa_mem_int32_ThComm",               # avec 1 thread de comm
			 "T_Numa_mem_int32_Funneled"],            # Thread Funneled
     
     'executables'   => ["example/bin/csimple",
			 "example/bin/zsimple",],

     'parametres'    => ["classique"],
     'matrices'    => [$PASTIX_SRC{"T"}."/matrix/young4c.mtx"], 
 

    },
    # Version ditribuee
    # Reel ; classique 
    {'libnames'    => ["T_Numa_mem_int64_dist_nosmp",           # Full-MPI Recv
		       "T_Numa_mem_int64_dist_nosmp_IRecv",     # Full-MPI IRecv
		       "T_Numa_mem_int64_dist_nosmp_ThComm",    # Full-MPI a 2 threads (calcul+Comm)
		       "T_Numa_mem_int64_dist",                 # MPI-SMP
		       "T_Numa_mem_int64_dist_IRecv",           # En IRecv
		       "T_Numa_mem_int64_dist_ThComm",          # avec 1 thread de comm
		       "T_Numa_mem_int64_dist_Funneled",        # Thread Funneled
		       "T_Numa_mem_int64_dist_dyn_ThComm",   # avec 1 thread de comm
		       "T_Numa_mem_int64_dist_dyn_Funneled", # Thread Funneled
		       "T_Numa_mem_int64_dist_trace",           # trace
		       "T_Numa_mem_int32_dist_nosmp",           # Full-MPI Recv
		       "T_Numa_mem_int32_dist_nosmp_IRecv",     # Full-MPI IRecv
		       "T_Numa_mem_int32_dist_nosmp_ThComm",    # Full-MPI a 2 threads (calcul+Comm)
		       "T_Numa_mem_int32_dist",                 # MPI-SMP
		       "T_Numa_mem_int32_dist_IRecv",           # En IRecv
		       "T_Numa_mem_int32_dist_ThComm",          # avec 1 thread de comm
		       "T_Numa_mem_int32_dist_Funneled"],       # Thread Funneled
  
     'executables'   => ["example/bin/simple_dist",
			 "example/bin/ssimple_dist",
			 "example/bin/dsimple_dist",],

     'parametres'    => ["classique"],
     'matrices'    => [$PASTIX_SRC{"T"}."/matrix/orsirr.rua",
                       $PASTIX_SRC{"T"}."/matrix/small.rsa"], 
 
    },
    # Version complexe
    # Reel ; classique 
    {'libnames'    => ["T_Numa_mem_int64_dist_nosmp",           # Full-MPI Recv
		       "T_Numa_mem_int64_dist_nosmp_IRecv",     # Full-MPI IRecv
		       "T_Numa_mem_int64_dist_nosmp_ThComm",    # Full-MPI a 2 threads (calcul+Comm)
		       "T_Numa_mem_int64_dist",                 # MPI-SMP
		       "T_Numa_mem_int64_dist_IRecv",           # En IRecv
		       "T_Numa_mem_int64_dist_ThComm",          # avec 1 thread de comm
		       "T_Numa_mem_int64_dist_Funneled",        # Thread Funneled
		       "T_Numa_mem_int64_dist_dyn_ThComm",   # avec 1 thread de comm
		       "T_Numa_mem_int64_dist_dyn_Funneled", # Thread Funneled
		       "T_Numa_mem_int64_dist_trace",           # trace
		       "T_Numa_mem_int32_dist_nosmp",           # Full-MPI Recv
		       "T_Numa_mem_int32_dist_nosmp_IRecv",     # Full-MPI IRecv
		       "T_Numa_mem_int32_dist_nosmp_ThComm",    # Full-MPI a 2 threads (calcul+Comm)
		       "T_Numa_mem_int32_dist",                 # MPI-SMP
		       "T_Numa_mem_int32_dist_IRecv",           # En IRecv
		       "T_Numa_mem_int32_dist_ThComm",          # avec 1 thread de comm
		       "T_Numa_mem_int32_dist_Funneled"],       # Thread Funneled
  
     'executables'   => ["example/bin/csimple_dist",
			 "example/bin/zsimple_dist",],

     'parametres'    => ["classique"],
     'matrices'    => [$PASTIX_SRC{"T"}."/matrix/young4c.mtx"], 
 
    },
    );



# database information
$db="pastix";
$host="localhost";
$userid="pastix";
$passwd="PaStiX";
$connectionInfo="dbi:mysql:$db;$host";
