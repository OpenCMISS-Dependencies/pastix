#!/usr/bin/perl

###############################################################################
#
#   Module: Compilation.pm
#
###############################################################################
#
#   This file contains routine to build compilation options.
#
#   Authors:
#     Mathieu Faverge - faverge@labri.fr
#     Xavier  Lacoste - lacoste@labri.fr
#
###############################################################################
package Modules::Compilation;
use Term::ANSIColor;
use strict;
use warnings;
use vars qw(@ISA @EXPORT);
use Exporter;
use Modules::Common;
use File::Spec;

our @ISA    = qw(Exporter);
our @EXPORT = qw( Exec_GetParam );

# 
# Function: Exec_GetParam
#
# Build the list of compilation option in fonction of the name
# of the library and of the machine option, defined in 
# Modules/Machines/_Machine_.pm.
#
# Here is the list of all keywords 
#  
#   _nompi   - Add the flad -DFORCE_NOMPI
#   _mpi     - Nothing, opposite to _nompi
#   _nosmp   - Add the flag -DFORCE_NOSMP
#   _smp     - Nothing, opposite to _nosmp
#   _static  - Nothing, opposite to _dyn
#   _dyn     - Add the flag -DPASTIX_DYNSCHED
#   _int     - Uses default integer type
#   _long    - Uses long integers
#   _int32   - Uses 32 bits integers
#   _int64   - Uses 64 bits integers
#   _simple  - Uses simple float precision
#   _double  - Uses double float precision
#   _real    - Uses real float
#   _complex - Uses complex float
#   _scotch  - Uses scotch
#   _metis   - Uses metis
#   _dist    - Uses distributed interface
#   _conso   - Uses -DFORCE_CONSO
#
# Parameters:
#   cas - The name of the librairy
#
# Returns:
#   path       - path to the branche
#   param_exec - Compilation parameters
#
sub Exec_GetParam #($cas)
{
    my ($cas) = @_ ;
    my %param_exec = ();
    my %option = () ;
    
    #####################################################
    #              Configuration par défaut
    #####################################################
    
    $param_exec{_MACH_NBNODES_} = 0;
    $param_exec{_MACH_NBCORES_} = 2;
    
    $param_exec{_HOSTARCH_}   = $machines{'default'}{'hostarch'};
    $param_exec{_VERSIONBIT_} = $machines{'default'}{'bits'};
    $param_exec{_CCPROG_}     = $machines{'default'}{'ccprog'};  
    $param_exec{_F77PROG_}    = $machines{'default'}{'f77prog'};
    $param_exec{_MCFPROG_}    = $machines{'default'}{'mcfprog'};	 
    $param_exec{_F90FLAGS_}   = $machines{'default'}{'f90flags'};  
    $param_exec{_CCFLAGS_}    = $machines{'default'}{'ccflags'};
    $param_exec{_LDFLAGS_}    = $machines{'default'}{'ldflags'};	 
    $param_exec{_CCFOPT_}     = $machines{'default'}{'ccfopt'};  
    $param_exec{_CCFDEB_}     = $machines{'default'}{'ccfdeb'};
    $param_exec{_LKFOPT_}     = $machines{'default'}{'lkfopt'};	 
    $param_exec{_MPCCPROG_}   = $machines{'default'}{'mpccprog'};  
    $param_exec{_ARFLAGS_}    = $machines{'default'}{'arflags'};
    $param_exec{_ARPROG_}     = $machines{'default'}{'arprog'};	 

    $option{mpi} = 1; # 0: nompi, 1: mpi
    $option{smp} = 1; # 0: nosmp, 1: smp
    $option{dyn} = 0; # 0: static, 1: dynamic
    $option{int} = 3; # 0: int, 1: long, 2: int32, 3: int64
    $option{prc} = 1; # 0: simple, 1: double
    $option{flt} = 0; # 0: real, 1: complexe
    $option{ord} = 0; # 0: Scotch, 1: Metis
    $option{dis} = 0; # 0: classical pastix only, 1 : -DDISTRIBUTED (2 executables)
    $option{fco} = 0; # 0: Froce conso 1 : -DFORCE_CONSO            

    #
    # Choix de la version a compiler
    #
    my $branche = $cas;
    $branche =~ s/([A-Z])_.*/$1/ ;
    my $path = $PASTIX_SRC{$branche};
    
    #print $branche." : ".$path."\n";

    #
    # Adaptation des options de compilation à chaque machine
    #

    if ($system =~ /Darwin/ ) #Mac
    {
	$param_exec{_MACH_NBCORES_} = 0;
	$param_exec{_HOSTARCH_} = "i686_mac";
	$param_exec{_LDFLAGS_}  = "-lblas -lg2c -lm";
	$param_exec{_F77PROG_}  = "g77 ";
    }
    elsif ( $machine =~ /borderline/ )
    {
	# Avec les modules load j'ai pas trop su comment me démerder avec la hastable machine
	my $modulecmd;
	$param_exec{_MACH_NBCORES_} = 16;
 	$param_exec{_VERSIONBIT_} = "_64bit";
        $param_exec{_LDFLAGS_}    = " -lm -lrt";

	$modulecmd ="/cvos/local/apps/environment-modules/3.2.6/bin/modulecmd ";
	$modulecmd.="bash -l list 2>&1 ";
	
	# Compilo Intel
	if ( `$modulecmd | grep "intel/cce/9.1.051"` && `$modulecmd | grep "intel/fce/9.1.051"`)
	{
	    print "Compilation avec icc\n";
	    $param_exec{_CCPROG_}   = "icc ";
	    $param_exec{_F77PROG_}  = "ifort";
	    $param_exec{_F90FLAGS_} = "-fpp";
	    $param_exec{_LDFLAGS_} .= " -L/cvos/shared/apps/intel/fce/9.1.051/lib -lifcore";
	}
	else # Compilo Gcc
	{
	    print "Compilation avec gcc\n";
	    $param_exec{_F77PROG_}  = "gfortran";
	    $param_exec{_F90FLAGS_} = "-ffree-form -x f95-cpp-input";
	    $param_exec{_LDFLAGS_} .= " -lgfortran";
	}
	    
	# Blas
	if ( `$modulecmd | grep "acml/gcc-int64/64/3.6.0"` )
	{
	    print "Utilisation de la libacml\n";
	    $param_exec{_LDFLAGS_}   .= " -lgfortran";
	    $param_exec{_LDFLAGS_}   .= " -L/cvos/shared/apps/acml/3.6.0/gfortran64_int64/lib -lacml";
	} 
	else 
	{
	    print "Il faut charger une librairie Blas avec module load\n";
	    print "Si c'est déjà le cas, il faut l'ajouter dans le fichier :\n";
	    print "$SCRIPT_PATH/Modules/Compilation.pm\n";
	}
	#$param_exec{_LDFLAGS_}    = " -L/cvos/shared/apps/gotoblas/opteron/64/1.11/ -lgoto";
    }
    elsif ( (defined($machines{$machine})) && (defined($machines{$machine}{'ccprog'})) )  
    {
	$param_exec{_HOSTARCH_}   = $machines{$machine}{'hostarch'};
	$param_exec{_VERSIONBIT_} = $machines{$machine}{'bits'};
	$param_exec{_CCPROG_}     = $machines{$machine}{'ccprog'};
	$param_exec{_F77PROG_}    = $machines{$machine}{'f77prog'};
	$param_exec{_MCFPROG_}    = $machines{$machine}{'mcfprog'};
	$param_exec{_F90FLAGS_}   = $machines{$machine}{'f90flags'};
	$param_exec{_CCFLAGS_}    = $machines{$machine}{'ccflags'};
	$param_exec{_LDFLAGS_}    = $machines{$machine}{'ldflags'};
	$param_exec{_CCFOPT_}     = $machines{$machine}{'ccfopt'};
	$param_exec{_CCFDEB_}     = $machines{$machine}{'ccfdeb'};
	$param_exec{_LKFOPT_}     = $machines{$machine}{'lkfopt'};
	$param_exec{_MPCCPROG_}   = $machines{$machine}{'mpccprog'};
	$param_exec{_ARFLAGS_}    = $machines{$machine}{'arflags'};
	$param_exec{_ARPROG_}     = $machines{$machine}{'arprog'};	 

	$param_exec{_MACH_NBCORES_} = $machines{$machine}{'nbcores'};
    }
    else 
    {
	print color 'yellow';
	print "WARNING : configuration par defaut ( machine = $machine )\n";
	print color 'reset';
    }

    #
    # Etude du nom de l'executable pour positionner les options
    # 
    $option{mpi} = 0 if ( $cas =~ /.*_nompi.*/ );
    $option{mpi} = 1 if ( $cas =~ /.*_mpi.*/ );
    $option{smp} = 0 if ( $cas =~ /.*_nosmp.*/ );
    $option{smp} = 1 if ( $cas =~ /.*_smp.*/ );
    $option{dyn} = 0 if ( $cas =~ /.*_static.*/ );
    $option{dyn} = 1 if ( $cas =~ /.*_dyn.*/ );
    $option{int} = 0 if ( $cas =~ /.*_int.*/ );
    $option{int} = 1 if ( $cas =~ /.*_long.*/ );
    $option{int} = 2 if ( $cas =~ /.*_int32.*/ );
    $option{int} = 3 if ( $cas =~ /.*_int64.*/ );
    $option{prc} = 0 if ( $cas =~ /.*_simple.*/ );
    $option{prc} = 1 if ( $cas =~ /.*_double.*/ );
    $option{flt} = 0 if ( $cas =~ /.*_real.*/ );
    $option{flt} = 1 if ( $cas =~ /.*_complex.*/ );
    $option{ord} = 0 if ( $cas =~ /.*_scotch.*/ );
    $option{ord} = 1 if ( $cas =~ /.*_metis.*/ );
    $option{dis} = 1 if ( $cas =~ /.*_dist.*/ );
    $option{fco} = 1 if ( $cas =~ /.*_conso.*/ );
    #
    # Si le $SCOTCH_HOME ne correspond pas au $option{int} on le corrige
    #

    my $type_entier = $SCOTCH_HOME;
    # suppression du / final si il y en a un
    $type_entier =~ s/(.*)\/$/$1/;
    # selection du dernier mot 
    $type_entier =~ s/.*\/([^\/]*)/$1/;
    my $type_voulu = "";
    $type_voulu = "int"   if (  $option{int} == 0);
    $type_voulu = "long"  if (  $option{int} == 1);
    $type_voulu = "int32" if (  $option{int} == 2);
    $type_voulu = "int64" if (  $option{int} == 3);
    if ($type_entier eq "int" || 
	$type_entier eq "long" || 
	$type_entier eq "int32" || 
	$type_entier eq "int64")
    {
	if ($type_voulu ne $type_entier)
	{
	    print "Changement de type  : " . $type_entier . "->" . $type_voulu    . "\n";
	    $SCOTCH_HOME =~ s/$type_entier/$type_voulu/;
	}
    }

    $param_exec{_CCTYPES_} = "";

    # Utilisation de MPI ou non 
    if ( $option{mpi} ){
	$param_exec{_VMPI_} = "VERSIONMPI  = _mpi";
    } else {
	$param_exec{_CCTYPES_} .= " -DFORCE_NOMPI";
	$param_exec{_VMPI_}  = 'VERSIONMPI  = _nompi'."\n";
	$param_exec{_VMPI_} .= 'MPCCPROG   := $(CCPROG)';
    }

    # Utilisation de SMP ou non
    if ( $option{smp} ){ 
	$param_exec{_SMP_} = "VERSIONSMP  = _smp";
    }
    else {
	$param_exec{_SMP_} = 'VERSIONSMP  = _nosmp';
	$param_exec{_CCTYPES_} .= " -DFORCE_NOSMP";
    }

    # Définition du type d'entiers utilisés
    if ( $option{int} == 1) {
	$param_exec{_CCTYPES_} .= " -DFORCE_LONG -DLONG";
	$param_exec{_INT_} = "VERSIONINT  = _long";
    } elsif ( $option{int} == 2) {
	$param_exec{_CCTYPES_} .= " -DFORCE_INT32 -DINTSIZE32";
	$param_exec{_INT_} = "VERSIONINT  = _int32";
    } elsif ( $option{int} == 3) {
	$param_exec{_CCTYPES_} .= " -DFORCE_INT64 -DINTSSIZE64";
	$param_exec{_INT_} = "VERSIONINT  = _int64";
    } else {
	$param_exec{_INT_} = "VERSIONINT  = _int";
    }

    # Définition de la précision utilisée
    if ( $option{prc} ) {
	$param_exec{_CCTYPESFLT_} .= " -DFORCE_DOUBLE -DPREC_DOUBLE";
	$param_exec{_PRC_} = 'VERSIONPRC  = _double';
    } else {
	$param_exec{_PRC_} = 'VERSIONPRC  = _simple';
    }

    # Utilisation de réel ou de complexe 
    if ( $option{flt} ) {
	$param_exec{_CCTYPESFLT_} .= " -DFORCE_COMPLEX -DTYPE_COMPLEX";
	$param_exec{_FLT_} = 'VERSIONFLT  = _complex';
    } else {
	$param_exec{_FLT_} = "VERSIONFLT  = _real";
    }
        
    # Utilisation de Metis ou de Scotch
    if ( $option{ord} ){
	$param_exec{_ORD_}      = 'VERSIONORD  = _metis';
	$param_exec{_CCFLAGS_} .= " -DMETIS";
	$param_exec{_CCFLAGS_} .= ' -I$(METIS_HOME)/include';
	$param_exec{_ORDPATH_}  = "METIS_HOME  = $METIS_HOME";  
	$param_exec{_LDFLAGS_} .= ' -L$(METIS_HOME)/lib -lmetis';  
    } else {
	$param_exec{_ORD_}      = "VERSIONORD  = _scotch";
	$param_exec{_CCFLAGS_} .= " -DWITH_SCOTCH";
	$param_exec{_CCFLAGS_} .= ' -I$(SCOTCH_HOME)/include';
	$param_exec{_ORDPATH_} .= "SCOTCH_HOME = $SCOTCH_HOME";
	
	if ( $option{dis} == 1)
	{
	    if ( $option{mpi} == 0 )
	    {
		print "L'option distribuée est incompatible avec NoMPI\n";
		exit;
	    }
	    $param_exec{_LDFLAGS_} .= ' -L$(SCOTCH_HOME)/lib -lptscotch -lscotcherrexit';
	    $param_exec{_CCFLAGS_} .= " -DDISTRIBUTED";
	}
	else
	{
	    $param_exec{_LDFLAGS_} .= ' -L$(SCOTCH_HOME)/lib -lscotch -lscotcherrexit';
	}
    }



    # Gestion de hwloc
    my $lstopodir = "";

    # If the installed directory is specified, we just take it
    if ( (defined($machines{$machine}{'lstopodir'})) &&
	 !($machines{$machine}{'lstopodir'} eq "" ))
    {
	$lstopodir = $machines{$machine}{'lstopodir'};
    }
    # We try to find where hwloc is installed
    else 
    {
	my @path = File::Spec->path();
	for my $file (map { File::Spec->catfile($_, "lstopo") } @path) {
	    if (-x $file)
	    {
		$lstopodir = $file;
		last;
	    }
	}
	if (length($lstopodir) > 0)
	{
	    $lstopodir = `dirname $lstopodir`;
	    chop $lstopodir;
	    $lstopodir =~ s/\/bin//;
	}

    }
    
    $param_exec{_HWLPATH_} = "";
    if (!($lstopodir eq ""))
    { 
	$param_exec{_HWLPATH_}  = "HWLOC_HOME  = $lstopodir";
	$param_exec{_CCFLAGS_} .= ' -I$(HWLOC_HOME)/include -DWITH_HWLOC';
	$param_exec{_LDFLAGS_} .= ' -L$(HWLOC_HOME)/lib -lhwloc';
    }

    if ( $option{fco} == 1)
    {
        $param_exec{_CCFLAGS_} .= " -DFORCE_CONSO";
    }
   
    # Allocation NoSMPRaff
    if ( $cas =~ /.*_noSMPRaff.*/ ){
	$param_exec{_CCFLAGS_} .= " -DNOSMP_RAFF";
    }

    # Allocation OOC
    if ( $cas =~ /.*_ooc.*/ ){
	$param_exec{_CCFLAGS_} .= " -DOOC_FTGT -DOOC_PERCENT_COEFTAB -DOOC_CLOCK -DOOC_DETECT_DEADLOCKS";
    }

    # Allocation Numa
    if ( $cas =~ /.*_Numa.*/ ){
	$param_exec{_CCFLAGS_} .= " -DNUMA_ALLOC";
    }

    # Affichage de la memoire utilisee
    if ( $cas =~ /.*_mem.*/ ){
	$param_exec{_CCFLAGS_} .= " -DMEMORY_USAGE";
    }
    if ( $cas =~ /.*_ovh.*/ ){
	$param_exec{_CCFLAGS_} .= " -DSTATS_SOPALIN";
    }

    if ( $cas =~ /.*_symbmtx.*/ ){
	$param_exec{_CCFLAGS_} .= " -DONLY_LOAD_SYMBMTX";
    }
    if ( $cas =~ /.*_compact.*/ ){
	$param_exec{_CCFLAGS_} .= " -DCOMPACT_SMX";
    }
    if ( $cas =~ /.*_ddt.*/ ){
	$param_exec{_CCFLAGS_} .= " -DDOUBLE_DOWN_TIME";
    }
    if ( $cas =~ /.*_diagdom.*/ ){
	$param_exec{_CCFLAGS_} .= " -DFORCE_DIAGDOM";
    }
    
    # Réception Non Bloquante 
    if ( $cas =~ /.*_IRecv.*/ ){
	$param_exec{_CCFLAGS_} .= " -DTEST_IRECV";
    }

    # Suppression des types MPI 
    if ( $cas =~ /.*_notypes.*/ ){
	$param_exec{_CCFLAGS_} .= " -DNO_MPI_TYPE";
    }

    # Thread de réception séparé
    if ( $cas =~ /.*_ThComm.*/ ){
	$param_exec{_CCFLAGS_} .= " -DTHREAD_COMM";
    }

    # Thread Funneled pour MPI n'acceptant pas le Multiple
    if ( $cas =~ /.*_ComOrd.*/ ){
	$param_exec{_CCFLAGS_} .= " -DCOMM_REORDER";
    }

    # Thread Funneled pour MPI n'acceptant pas le Multiple
    if ( $cas =~ /.*_Funneled.*/ ){
	$param_exec{_CCFLAGS_} .= " -DPASTIX_FUNNELED";
    }

    # 2D Dynamique
    if ( $cas =~ /.*_esc.*/ ){
      $param_exec{_CCFLAGS_} .= " -DPASTIX_ESC";
    }

    # Ordo Dynamique
    if ( $cas =~ /.*_dyn.*/ ){
      $param_exec{_CCFLAGS_} .= " -DPASTIX_DYNSCHED";
    }

    # Trace
    if ( $cas =~ /.*_trace.*/ ){
	$param_exec{_CCFLAGS_} .= " -DTRACE_SOPALIN";
    }

    # Définition de la flavor à utiliser
    $param_exec{_PM2_FLAVOR_} = "";
    if ( $cas =~ /.*_FL([-a-z0-9]*)FL.*/ ){
	$param_exec{_PM2_FLAVOR_} = "$1";
    }    

    # Utilisation des thread Marcels à la place des threads POSIX
    if ( $cas =~ /.*_Marcel.*/ || $cas =~ /.*_dynSched.*/ ){
	if ( $param_exec{_PM2_FLAVOR_} eq ""){
	    print "ERROR : Flavor non définie";
	    exit;
	}
	$param_exec{_LDFLAGS_} .= " `pm2-config --libs`";
        $param_exec{_CCFLAGS_} .= " `pm2-config --cflags` -I$ENV{ PM2_ROOT}/marcel/include/pthread";
    } else {
        $param_exec{_LDFLAGS_} .= " -lpthread";
    }

    # Utilisation de l'ordo à bulles
    if ( $cas =~ /.*dyn3.*/ ){
	$param_exec{_CCFLAGS_}   .= " -DTEST_REPART";
    }

    # Utilisation des bulles avec les ordo de Marcel
    if ( $cas =~ /.*_dynSched.*/ ){
	$param_exec{_DYN_}      = "VERSIONSCH  = _dyn"; 
        $param_exec{_CCFLAGS_} .= " -DPASTIX_DYNSCHED -DPASTIX_BUBBLESCHED";
    }
    else
    {
	$param_exec{_DYN_} = "VERSIONSCH = _static"; 
    }
    return ($path, %param_exec);
}
    
1;
