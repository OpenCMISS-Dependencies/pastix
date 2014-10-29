#!/usr/bin/perl

###############################################################################
#
#   Script: Pastix.pl
#
###############################################################################
#
#   The main run launching script.
#   It allows to compile executables and runs them.
#   It is also possible to generate Latex output.
#
#   Usage:
#    Pastix.pl [ -h ] [ -A | -cfe ]
#       -h    Affiche cette aide
#       -H    Affiche les options du scripts
#       -A    Equivalent a -cfe
#       -c    Compile les différentes versions
#       -f    Crée les répertoires et les fichiers des différents cas
#       -e    Exécute l'ensemble des cas
#       -o    Spécifie le fichier de sortie de la liste de soumission
#       -b    Soumet en 'besteffort' pour OAR
#       -g    Génère un tableau de résultat en latex
#       -p    Permet de spécifier le reseau pour Madeleine
#       -m    Mail errors
#       -d    Fill/Check database
#       -v    Stoppe la compilation après chaque compilation avec erreur
#       -V    Passe IPARM_VERBOSE à 2
#       -s x  IO Strategy (1 : Load / 2 : Save / 3 : Save on first, 
#             Load on next one)
#
#   Authors:
#     Mathieu Faverge - faverge@labri.fr
#     Xavier  Lacoste - lacoste@labri.fr
#
###############################################################################

use POSIX;
use strict;
#use warnings;
use Getopt::Std;
use Term::ANSIColor;
use FindBin qw($Bin);
use lib "$Bin";
use Modules::Common;
use Modules::Compilation;
use Modules::GenLatex;
use Modules::APIparser;
use Modules::LogParser;
use Modules::PastixDatabase;

###############################################################################
# Group: Variables

# hashtable: opts
#   commande line argument options    
my %opts;
# boolean: besteffort
#   indicate if we want to use besteffort
my $besteffort = 0;
# boolean: multiprocbynode
#   indicate if we have to indicate the number of task by node
my $multiprocbynode = 0;
# boolean: verbosemode
#   indicate if we want to interact with the script, stopping for errors and 
#   warning, etc.
my $verbosemode = 0;
# boolean: mailingmode
#   indicate if we want to send the result of the rune by mail
my $mailingmode = 0;
# boolean: databasemode
#   indicate if we want to fill/check the database
my $databasemode = 0;
# string: mailtxt
#   the text of the sended mail
my $mailtxt     = "";
# integer: comp_tot
#   the number of compilation to perform
my $comp_tot    = 0;
# integer: comp_err
#   the number of compilation which failed
my $comp_err    = 0;
# integer: exec_tot
#   the number of execution to perform
my $exec_tot    = 0;
# integer: exec_err
#   the number of execution which failed
#
my $exec_err    = 0;
# integer: exec_war
#   the number of execution which we have to care about
#
my $exec_war    = 0;
# array: lstmetr
#   The list of the metrique that we want in output
my @lstmetr = ( 'exemple', 'nproc', 'nthrd', 'fact', 'solv', 'refin', 'iter', 'norm', 'mem', 'mend');
# string: sendmail
#   The sendmail command
my $sendmail = "/usr/sbin/sendmail -t";


###############################################################################
# Group: Functions

#
# Function: Usage
#
# Prints usage.
#
sub Usage #()
{

    print "Pastix.pl [ -h ] [ -A | -cfe ]\n";
    print "   -h    Affiche cette aide\n";
    print "   -H    Affiche les options du scripts\n";
    print "   -A    Equivalent a -cfe\n";
    print "   -c    Compile les différentes versions\n";
    print "   -f    Crée les répertoires et les fichiers des différents cas\n";
    print "   -e    Exécute l'ensemble des cas\n";
    print "   -o    Spécifie le fichier de sortie de la liste de soumission\n";
    print "   -b    Soumet en 'besteffort' pour OAR\n";
    print "   -g    Génère un tableau de résultat en latex\n";
    print "   -p    Permet de spécifier le reseau pour Madeleine\n";
    print "   -m    Mail errors\n";
    print "   -d    Fill/Check database\n";
    print "   -v    Stoppe la compilation après chaque compilation avec erreur\n";
    print "   -V    Passe IPARM_VERBOSE à 2\n";
    print "   -s x  IO Strategy (1 : Load / 2 : Save / 3 : Save on first, Load on next one)\n";
}

# 
# Function: contain
#
# Returns 1 if a string appear in an array of string, 0 otherwise.
#
# Parameters:
#   string - The string searched for.
#   array  - The array in which we search.
#
sub contain #($string, @array)
{
    my ($string, @array) = @_;
    my $count = 0;
    foreach my $string2 (@array)
    {
	return 1 if ($string eq $string2);
    }
    return 0;
}
#
# Function: min
#
# Returns the minimum value between two values.
#
# Parameters: 
#   a - First value to compare
#   b - Second value to compare
sub min #($a, $b)
{
    my ($a, $b) = @_ ;
    return (($a > $b) ? $b : $a);
}


#
# Function: Compilation
#
# Compilate the different executables
#
sub Compilation()
{
    my $user = `whoami`;
    chomp $user;
    my $logfile;
    my $ligne;
    my $stop      = 0;
    my $configtpl = "$TPL_ROOT/config.tpl";
    my $branche;
    my $comp_done;
    my %liste_compilations_globale = ();

    $mailtxt   .= "<h3>Compilation</h3>\n";
    # Boucle sur les cas
    if ($databasemode)
    {
	checkMachinePresence($machine);
    }

    for my $key (keys(%PASTIX_SRC))
    {
	my $path_src = $PASTIX_SRC{$key};
	if (-f "$path_src/config.in")
	{
	    system("mv $path_src/config.in $path_src/config.in.saved.$date");
	}
    }


    # construction de la hastable de compilation
    foreach my $run_ref (@liste_runs)
    {
	my %run = %$run_ref;
	my $list_lib_ref  = $run{'libnames'};
	my $list_exec_ref = $run{'executables'};
	my @array;
	foreach my $lib (@$list_lib_ref)
	{	
	    foreach my $executable (@$list_exec_ref)
	    {
		
		my $liste_ref = $liste_compilations_globale{$lib};
		my $found = contain($executable, @$liste_ref);
		if ($found == 0)
		{
		    push(@$liste_ref, $executable);
		    $comp_tot = $comp_tot + 1;
		}
		$liste_compilations_globale{$lib} = $liste_ref;
	    }
	}
	
    }
    $comp_done = 0;
    for my $libname (keys(%liste_compilations_globale))
    {

      repeat:
	# Compilation
	print "\n";
	print "--- Compil $libname (".($comp_done + 1)."/$comp_tot)---\n";
	print "\n";
	$mailtxt .= "<p>Compilation de $libname (".($comp_done + 1)."/$comp_tot) :</p>\n";
	
	# Positionnement des variables 
	my ($path_src, %param_exec) = Exec_GetParam($libname);
	my $version = `cd $path_src; svnversion`;
	
	my $err = 0;
	
	my $path_exec = $RESULT_PATH.'/'.$user.'/'.$machine.'/'.$libname;
	
	my $MAKE    = "PM2_FLAVOR=$param_exec{_PM2_FLAVOR_} $machines{$machine}{'makecmd'} ";
	#my $sorties = "2> $path_exec/$date.err > $path_exec/$date.out";
	my $sorties = ">> $path_exec/$date.log 2>&1";
	
	# Creation du repertoire d'exécution 
	`mkdir -p $path_exec`;
	chdir "$path_exec";
	
	# Recuperation de la version utilisée
	`(cd $path_src && svn info > $path_exec/revision)`;
	
	# création du fichier config.in 
	MakeTemplate($configtpl, "$path_exec/config.in", %param_exec);
	
	# mise en place du lien dans le rep de compilation
	`(cd $path_src && ln -sf $path_exec/config.in)`;
	
	# Compilation
	print "-- make clean   --\n";
	$err = system("cd $path_src && $MAKE  -j 8 clean $sorties");
	#$err = system("cd $path_src && $MAKE clean $sorties");
	print "-- make print_options   --\n";
	$err = system("cd $path_src && $MAKE  -j 8 utils/bin/".
		      $machines{$machine}{'hostarch'}.
		      "/print_options $sorties") if ($err == 0);
	system("svn diff $path_src > /$path_exec/patch_$date;");
	$err = system("$path_src/utils/bin/*/print_options > $path_exec/log_po_$date") if ($err == 0);
	system("cp $path_exec/config.in config.in_$date");

	my $lib_list_ref = $liste_compilations_globale{$libname};
	if ($databasemode)
	{
	    $branche = `basename $path_src`;
	    chomp $branche;
	    
	    foreach my $execname (@$lib_list_ref)
	    {
		my $execname_base = `basename $execname`; chop $execname_base;
		fillExecutable($execname_base,  $branche, "$path_exec/patch_$date",
			       "$path_exec/log_po_$date", "", "$path_exec/config.in_$date");
	    }
	}
	
	if ($err != 0)
	{
	    $mailtxt .= " <font color=red>&eacute;chou&eacute;e, erreur print_options</font></p>\n";	    
	    $comp_err = $comp_err + scalar(@$lib_list_ref);
	    next;
	}
	    
	print "-- make expor   --\n";
	$err = system("echo \"-- make expor    --\" $sorties")           if ($err == 0);
	$err = system("echo \"cd $path_src && $MAKE -j 8 \" $sorties")   if ($err == 0);
	$err = system("cd $path_src && $MAKE -j 8 $sorties")            if ($err == 0);
	#$err = system("cd $path_src && $MAKE expor $sorties")            if ($err == 0);
	print "-- make install --\n";
	$err = system("echo \"-- make install  --\" $sorties")           if ($err == 0);
#	$err = system("echo \"cd $path_src && $MAKE install\" $sorties") if ($err == 0);
#	$err = system("cd $path_src && $MAKE install $sorties")          if ($err == 0);
	print "-- make drivers --\n";
	$err = system("echo \"-- make drivers  --\" $sorties")           if ($err == 0);
	$err = system("echo \"cd $path_src && $MAKE -j 8  drivers\" $sorties") if ($err == 0);
	$err = system("cd $path_src && $MAKE  -j 8 drivers $sorties")          if ($err == 0);
	#$err = system("cd $path_src && $MAKE drivers $sorties")          if ($err == 0);
	$err = system("echo \"-- make examples --\" $sorties")           if ($err == 0);
	print "-- make examples    --\n";
	
	# Boucle sur les executables et leur parametres
	foreach my $execname (@$lib_list_ref)
	{
	    print "$execname ";
	    my $compil;
	    $execname = `basename $execname`;
	    chop $execname;

	    $err = system("echo \"cd $path_src/example/src && $MAKE ../bin/$execname\" $sorties") if ($err == 0);
	    $err = system("cd $path_src/example/src && $MAKE ../bin/$execname $sorties") if ($err == 0);
	    
	    # Copie de l'exécutable dans le rep
	    
	    my $nberrEN = `cat $path_exec/$date.log | grep -i "error[ ]*:" | wc -l`;
	    chop $nberrEN;
	    my $nberrFR = `cat $path_exec/$date.log | grep -i "erreur[ ]*:" | wc -l`;
	    chop $nberrFR;
	    my $nberr = $nberrFR + $nberrEN;
	    
	    my $nbwarnEN = `cat $path_exec/$date.log | grep -i "warning[ ]*:" | grep -v /usr/include/bits/string3.h| wc -l`;
	    chop $nbwarnEN;
	    my $nbwarnFR = `cat $path_exec/$date.log | grep -i "attention[ ]*:" | grep -v /usr/include/bits/string3.h| wc -l`;
	    chop $nbwarnFR;
	    my $nbwarn = $nbwarnFR + $nbwarnEN;
	    
	    $err = system("cp $path_src/example/bin/$execname  $path_exec/ $sorties") if ($err == 0);
	    
	    if (  ($nberr+$nbwarn+$err) == 0 ){
		print color 'green';
		$compil = "réussie";
		if ($databasemode)
		{
		    $logfile = fillCompiler($execname,  $err, $nberr, $nbwarn, 
					    "$path_exec/$date.log", 
					    $date, $user);
		    $mailtxt .= "<p><a href=\"http://vulcain.local/regression/".
			$logfile."\"> $execname </a>";
		}
		else
		{
		    $mailtxt .= "<p> $execname ";
		}
		
		$mailtxt .= " <font color=green>r&eacute;ussie</font></p>\n";
		
	    }
	    elsif ($nberr+$err == 0)
	    {
		print color 'yellow';
		$compil = "warnings";
		$stop   = 1;
		if ($databasemode)
		{
		    $logfile = fillCompiler($execname,  0, $nberr, $nbwarn,
					    "$path_exec/$date.log", $date, $user);
		    $mailtxt .= "<p><a href=\"http://vulcain.local/regression".
			$logfile."\"> $execname </a>";
		}
		else
		{
		    $mailtxt .= "<p> $execname ";
		}
		$mailtxt .= " <font color=orange>$nbwarn Warnings</font></p>\n";
	    }
	    else
	    {
		print color 'red';
		$compil = "échouée";
		$stop   = 1;
		$comp_err = $comp_err +1;
		if ($databasemode)
		{
		    $logfile = fillCompiler($execname,  0,  $nberr, $nbwarn,
					    "$path_exec/$date.log", 
					    $date, $user);
		    $mailtxt .= "<p><a href=\"http://vulcain.local/regression/".
			$logfile."\"> $execname </a>";
		}
		else
		{
		    $mailtxt .= "<p> $execname ";
		}
		$mailtxt .= " <font color=red>&eacute;chou&eacute;e</font></p>\n";
	    }
	    
	    print "Compilation $compil (W : $nbwarn, E : $nberr)\n";
	    print "$path_exec/$date.log\n";
	    print color 'reset';
	    if ( $verbosemode )
	    {
		if ( $nberr > 0 )
		{
		    print "--- Erreurs ---\n";
		    system("cat $path_exec/$date.log | grep -i 'error[ ]*:' 2> /dev/null");
		    system("cat $path_exec/$date.log | grep -i 'erreur[ ]*:' 2> /dev/null");
		}
		if ( $nbwarn > 0 )
		{
		    print "--- Warnings ---\n";
		    system("cat $path_exec/$date.log | grep -i 'warning[ ]*:' | grep -v /usr/include/bits/string3.h 2> /dev/null");
		    system("cat $path_exec/$date.log | grep -i 'attention[ ]*:' | grep -v /usr/include/bits/string3.h 2> /dev/null");
		}
		
		while ( $stop )
		{
		    print '(q : quit, c : continue, r : repeat, o : show output) : ';
		    $ligne = <STDIN>;
		    chop $ligne;
		    if ( $ligne eq "q" )
		    {
			exit;
		    } 
		    elsif ( $ligne eq "c" )
		    {
			$stop = 0;
			next;
		    }
		    elsif ( $ligne eq "r" )
		    {
			system("/bin/rm -f $path_exec/$date.log");
			$stop = 0;
			goto repeat;
		    }
		    elsif ( $ligne eq "o" )
		    {
			system("less $path_exec/$date.log");
			next;
		    }
		    else
		    {
			next;
		    }
		}
	    }
	    # Print errors messages in the mail
	    if ($err > 0)
	    {
		$mailtxt .= "<p><font color=red>Compilation returned : $err</font></p>\n"; 
	    }
	    if ( $nberr > 0 )
	    {
		$mailtxt .= "<h5><font color=red>Erreurs</font></h5>\n";
		my $errorlist = `cat $path_exec/$date.log | grep -i 'error[ ]*:' 2> /dev/null`;
		$errorlist .= `cat $path_exec/$date.log | grep -i 'erreur[ ]*:' 2> /dev/null`;
		$errorlist =~ s/\n/<br>/g;
		$mailtxt .= "<p>$errorlist<br></p>\n";
	    }
	    
	    if ( $nbwarn > 0 )
	    {
		$mailtxt .= "<h5> <font color=orange>Warnings</font></h5>\n";
		my $warninglist = `cat $path_exec/$date.log | grep -i 'warning[ ]*:' 2> /dev/null`;
		$warninglist .= `cat $path_exec/$date.log | grep -i 'attention[ ]*:' 2> /dev/null`;
		$warninglist =~ s/\n/<br>/g;
		$mailtxt .= "<p>$warninglist<br></p>\n";
	    }
	    $comp_done++;
	    
	} # For each execname
    } #for each libname

    for my $key (keys(%PASTIX_SRC))
    {
	my $path_src = $PASTIX_SRC{$key};
	if (-f "$path_src/config.in.saved.$date")
	{
	    system("mv $path_src/config.in.saved.$date $path_src/config.in");
	}
    }
    
   
}

####################################################
#               Matrices                           #
####################################################

## Modifie les paramètres de facto et de symmétrie en fonction 
## de la matrice et initialise les autres parametres
sub Matrice_GetParam {

    my ($cas, @iparms) = @_ ;
    my %param_matrix = ();

    my $link    = "rsaname";
    my $driver  = "-rsa $cas";
    my $complex = 0;

    if ( $cas =~ /\.rua/ or $cas =~ /matr5/ or $cas =~ /orsirr/ or $cas =~ /ultrasound/ or $cas =~ /MHD/ ){
	@iparms = setiparm("IPARM_FACTORIZATION",    2, @iparms); # LU
	@iparms = setiparm("IPARM_SYM",              1, @iparms); # Non-Symmétrique
    }

    if ( $cas =~ /\.mm/ || $cas =~ /\.mtx/ ){
	$link   = "mmname";
	$driver = "-mm $cas";
    }	
    if ( $cas =~ /\.ijv/ ) {
	$link   = "ijvname";
	$driver = "-ijv $cas";
    }	

    if ( $cas =~ /Haltere/ or $cas =~ /Amande/ or $cas =~ /new_15/  or $cas =~ /empty/ or $cas =~ /young4c/ or $cas =~ /10Millions/  or 
	 $cas =~ /vfem/ or $cas =~ /3Dspectralwave/ or $cas =~ /fem_hifreq_circuit/ or $cas =~ /RFdevice/ or $cas =~ /mono_500Hz/)
    {
	$complex = 1;
    }

    return ( $link, $driver, $complex, @iparms );
}

#
# Function: CreateDirAndFiles
#
# Create the result tree
#
sub CreateDirAndFiles {

    my $matrice;
    my $user = `whoami`;
    chomp $user;
    my @iparm = ();
    my @dparm = ();
    my @api   = ();
    my $lib_list_ref;
    # Boucle sur les cas
    foreach my $run_ref (@liste_runs)
    {
	my %run = %$run_ref;

	$lib_list_ref = $run{'libnames'};
	foreach my $libname (@$lib_list_ref)
	{
	    print " ---- $libname \n";
	    
	    my ($path_src, %param_exec) = Exec_GetParam($libname);
	    my $filename = $path_src."/common/src/api.h";
	    if ( ! -f $filename )
	    {
		print "$filename n'existe pas !!!\n";
		exit;
	    }

            # Boucles sur les Matrices
	    #
            # Pour chaque matrice on construit le repertoire 
	    # $RESULT_PATH/$user/$machine/$libname/$mat
	    my $list_mat_ref = $run{'matrices'};
	    foreach my $matrice (@$list_mat_ref)
	    {
		my $mat = `basename $matrice`; chop $mat;
				
		# Suppression du suffixe 
		$mat =~ s/\..*//;
		
		print " ------ $mat \n";
		
		# Creation du repertoire de la matrice
		`mkdir -p $RESULT_PATH/$user/$machine/$libname/$mat`;
		
		
		# pour chaque jeu de parametres,
		# on construit iparm et dparm
		# On crée pour chaque combianison de thread et de proc
		# un dossier  $RESULT_PATH/$user/$machine/$libname/$mat/$cas
		# ou $cas depend de T de P et de la liste de parametres
		# et on sauve iparm et dparm dedans...
		my $list_params_ref = $run{'parametres'};
		foreach my $param_list_name (@$list_params_ref)
		{
		    print " -------- $param_list_name \n";
		    my %liste_param = %{$liste_cas{$param_list_name}};

		    # Création des [id]parm / A faire obligatoirement apres le choix de matrice
		    # sinon il y a un risque de passer en LU sur une symmetrique
		    my ($refiparm, $refdparm, $refapi) = API_parser($filename);
		    @iparm = @$refiparm;
		    @dparm = @$refdparm;
		    @api   = @$refapi; 
		    
		    if (defined($liste_param{'iparm'}))
		    {
			my %liste_iparm = %{$liste_param{'iparm'}};
			foreach my $liparm (keys %liste_iparm)
			{
			    @iparm = setiparm($liparm, $liste_param{'iparm'}{$liparm}, @iparm);
			}
		    }
		    
		    if (defined($liste_param{'dparm'}))
		    {
			my %liste_dparm = %{$liste_param{'dparm'}};
			foreach my $ldparm (keys %liste_dparm)
			{
			    @dparm = setiparm($ldparm, $liste_param{'dparm'}{$ldparm}, @dparm);
			}
		    }
		    my $withricar = getiparm("IPARM_INCOMPLETE", @iparm);
		    if ($withricar == 1)
		    {
			@dparm = setdparm("DPARM_EPSILON_REFINEMENT", "1e-7", @dparm);
		    }
		    # Fin creation [id]parm
		    
		    # Recupération des parametres liés à la matrice
		    my ($link, undef, $complex, @iparm) = Matrice_GetParam($matrice, @iparm);
		    if ( $libname =~ /.*_symbmtx.*/ ) {
			$link  = "symbmtx";
		    }
		    
		    my $level = getiparm("IPARM_LEVEL_OF_FILL", @iparm);
		    if ( ( $libname =~ /.*_metis.*/ ) && ( $level!=-1 ) )
		    {
			print color 'red' if (!defined($opts{c}));
			print "METIS sans KASS n'est pas possible.\n";
			print color 'reset' if (!defined($opts{c}));
			next;
		    }
		    # Boucle sur les cas en thread et proc
		    for ( my $nbthread = $liste_param{'MIN_THREAD'} ;  
			  ($nbthread < $liste_param{'MAX_THREAD'}+1);
			  $nbthread = $nbthread * 2){
			next if ( ( $libname =~ /.*_nosmp.*/ ) && ($nbthread != 1) );
			
			# Affectation du nombre de thread dans iparm.txt
			@iparm = setiparm("IPARM_THREAD_NBR", $nbthread, @iparm);
			for ( my $nbproc = $liste_param{'MIN_PROC'}; 
			      ($nbproc < $liste_param{'MAX_PROC'}+1);
			      $nbproc = $nbproc * 2){
			    next if ( ( $libname =~ /.*_nompi.*/ ) && ($nbproc != 1) );
			    next if ($nbthread*$nbproc > $liste_param{'MAX_PROCTOTAL'});
			    my $cas = sprintf("%02dP_%02dT_%s", 
					      $nbproc, $nbthread,$param_list_name);
			    print $cas." ";
			    
			    #Creation du répertoire du cas
			    my $path="$RESULT_PATH/$user/$machine/$libname/$mat/$cas";
			    `mkdir -p $path`;
			    
			    # Creation du fichier dparm.txt
			    printDparm("$path/dparm.txt",@dparm);
			    
			    # Creation du fichier iparm.txt
			    @iparm = setiparm("IPARM_NB_SMP_NODE_USED", 
					      min(
						  ceil( ($nbproc * min($nbthread, $param_exec{_MACH_NBCORES_}))
							/ $param_exec{_MACH_NBCORES_}),
						  $param_exec{_MACH_NBNODES_}),
					      @iparm);
			    printIparm("$path/iparm.txt",@iparm);
				    
			    
			    #Creation du lien sur l'executable
			    my $list_execs_ref =  $run{'executables'};
			    foreach my $execname (@$list_execs_ref)
			    {
				$execname = `basename $execname`;
				chop $execname;
				`(cd $path && ln -sf $RESULT_PATH/$user/$machine/$libname/$execname)`;
				
				#Creation du lien sur la matrice
				`(cd $path && ln -sf $matrice $link)`;
			    }
			} # Fin Boucle proc
			
		    } #Fin boucle thread
		    
		    print "\n";
		} #Fin boucle params
		print "\n";
	    } #Fin de boucle sur la matrice
	} # Fin boucle sur la lib
    } # fin de boucle sur le jeu de runs
}


#
# Function: ExecuteJobs
#
# Runs the jobs
#
sub ExecuteJobs {

    my ($file) = @_;
    my $user = `whoami`;
    chomp $user;
    my $matrice;
    my %param_job = ();
    my $j;
    my $total = 0;
    my $done  = 0;
    my $timeout;
    my $deadlock = 0;
    my $lib_list_ref; 

    $mailtxt .= "<h3>Execution</h3>\n";
    open(M, ">>".$file);

    # Sauvegarde du rep courant
    my $rep = `pwd`; chop $rep;
			
    $param_job{_USER_} = $USERMAIL;
 
    foreach my $run_ref (@liste_runs)
    {
	my %run = %$run_ref;
	my $list_lib_ref = $run{'libnames'};
	my @list_lib = @$list_lib_ref; 
	$total = $total + $#list_lib+1; 
	
    }

    # Boucle sur les cas
    foreach my $run_ref (@liste_runs)
    {
	my %run = %$run_ref;
	$lib_list_ref = $run{'libnames'};

	foreach my $libname (@$lib_list_ref)
	{
	    print "--- Exec : $libname (".($done+1)."/$total) ---\n";
	    $mailtxt .= "<h4>Exec : $libname (".($done+1)."/$total)</h4>\n";

	    # Positionnement des variables 
	    my ($path_src, %param_exec) = Exec_GetParam($libname);
	    
	    # Boucles sur les Matrices
	    my $liste_matrices_ref = $run{'matrices'};
	    my @liste_matrices = @$liste_matrices_ref;
	    
	    my $total_mat = scalar(@liste_matrices);
	    my $j = -1;
	    foreach my $matrice (@liste_matrices)
	    {
		$j = $j+1;
		my $mat = `basename $matrice`;
		chop $mat;
		
		my $matwoext = $matrice;
		$matwoext =~ s/(.*)\..*/\1/; 

		if ($mat =~ /small/ ||
		    $mat =~ /orsirr/ ||
		    $mat =~ /young4c/)
		{
		    $timeout = 10;
		    if ($libname =~ /.*_ooc.*/)
		    {
			$timeout = 200;
		    }
		    
		}
		else
		{
		    $timeout = 0;
		}
		
		# Suppression du suffixe 
		$mat =~ s/\..*//;
		if ($databasemode)
		{
		    checkMatrixPresence($mat);
		}
		
		print "----- Matrice : $mat (".($j+1)."/$total_mat)\n";
		print " $matwoext\n";
		print " ".LogPrintMetriquesHeader(@lstmetr)."\n" if ($verbosemode);
		$mailtxt .= "<h5>Matrice : $mat (".($j+1)."/$total_mat)</h5>\n";
		
		# Recupération des parametres liés à la matrice
		my ($link, $driver, $complex, undef) = Matrice_GetParam($matrice);

		
		
		# Boucle sur les jeu de parametres
		my $list_params_ref = $run{'parametres'};
		foreach my $param_list_name (@$list_params_ref)
		{
		    my %liste_param = %{$liste_cas{$param_list_name}};
		    
		    $mailtxt .= "<h6>Parametres : $param_list_name </h6>\n";
		    $mailtxt .= "<table><tr><td></td>".LogPrintMetriquesHeaderHtml(@lstmetr)."</tr>\n";
		    my $filename = $path_src."/common/src/api.h";
		    if ( ! -f $filename )
		    {
			print "$filename n'existe pas !!!\n";
			exit;
		    }

		    # Création des [id]parm / A faire obligatoirement apres le choix de matrice
		    # sinon il y a un risque de passer en LU sur une symmetrique
		    my ($refiparm, $refdparm, $refapi) = API_parser($filename);
		    my @iparm = @$refiparm;
		    my @dparm = @$refdparm;
		    my @api   = @$refapi; 
		    
		    if (defined($liste_param{'iparm'}))
		    {
			my %liste_iparm = %{$liste_param{'iparm'}};
			foreach my $liparm (keys %liste_iparm)
			{
			    @iparm = setiparm($liparm, $liste_param{'iparm'}{$liparm}, @iparm);
			}
		    }
		    
		    if (defined($liste_param{'dparm'}))
		    {
			my %liste_dparm = %{$liste_param{'dparm'}};
			foreach my $ldparm (keys %liste_dparm)
			{
			    @dparm = setiparm($ldparm, $liste_param{'dparm'}{$ldparm}, @dparm);
			}
		    }
		    my $withricar = getiparm("IPARM_INCOMPLETE", @iparm);
		    if ($withricar == 1)
		    {
			@dparm = setdparm("DPARM_EPSILON_REFINEMENT", "1e-7", @dparm);
			}
		    # Fin creation [id]parm			
		    my $level        = getiparm("IPARM_LEVEL_OF_FILL",      @iparm);
		    my $amalg        = getiparm("IPARM_AMALGAMATION_LEVEL", @iparm);
		    my $ooc_limit    = getiparm("IPARM_OOC_LIMIT",          @iparm);
		    my $with2d       = getiparm("IPARM_DISTRIBUTION_LEVEL", @iparm);
		    my $verboselevel = getiparm("IPARM_VERBOSE",            @iparm)+1;
		    

		    if ( ( $libname =~ /.*_metis.*/ ) && ( $level!=-1 ) )
		    {
			print "METIS sans KASS n'est pas possible.\n";
			$mailtxt .= "<p><font color=orange>METIS sans KASS n'est pas possible.</font><br></p>";
			next;
		    }
		    # Boucle sur les cas en thread et proc
		    for ( my $nbproc   = $liste_param{'MAX_PROC'}; 
			  ($nbproc >= $liste_param{'MIN_PROC'}); 
			  $nbproc = $nbproc / 2)
		    {
			
			next if ( ( $libname =~ /.*_nompi.*/ ) && ($nbproc != 1) );
			
			for ( my $nbthread = $liste_param{'MIN_THREAD'};  
			      ($nbthread < $liste_param{'MAX_THREAD'}+1) 
			      && ($nbthread*$nbproc < $liste_param{'MAX_PROCTOTAL'}+1);
			      $nbthread = $nbthread * 2)
			{
			    
			    next if ( ( $libname =~ /.*_nosmp.*/ ) && ($nbthread != 1) );
			    
			    my $cas = sprintf("%02dP_%02dT_%s", 
					      $nbproc, $nbthread, $param_list_name);
			    
			    
			    print $cas." " if (! $verbosemode );
			    
			    my $path = "$RESULT_PATH/$user/$machine/$libname/$mat/$cas/";
			    
			    #Creation du répertoire du cas
			    system("mkdir -p $path/$date");
			    system("cp $path/*.txt       $path/$date/");
			    my $list_execs_ref =  $run{'executables'};
			    foreach my $execname (@$list_execs_ref)
			    {
				
				$execname = `basename $execname`;
				chop $execname;
				
				system("ln -fs $path/$execname  $path/$date/");
				
				if ( $libname =~ /.*_symbmtx.*/ ) {
				    system("ln -fs $path/symbmtx $path/$date/");
				}
				else {
				    system("ln -fs $path/*name   $path/$date/");
				}

				if ( getiparm("IPARM_IO_STRATEGY", @iparm) == 1 )
				{
				    system("ln -fs $matwoext/ordername $path/$date/");
				    system("ln -fs $matwoext/graphname $path/$date/");
				    system("ln -fs $matwoext/symbmtx   $path/$date/");				    
				}

				my $options = "";
				#next if (($complex && !($execname =~ /[cz].*/))
				#	 || (!$complex && ($execname =~ /[cz].*/)));
				
				$param_job{_NAME_}        = "Pastix_$cas";
				$param_job{_MULTITHREAD_} = "yes";
				$param_job{_NBPROC_}      = $nbproc;
				$param_job{_NBTHREAD_}    = $nbthread;
				$param_job{_NBNODE_}      = $nbproc;
				$param_job{_TASKPERNODE_} = 1;
				$param_job{_PATH1_}       = $RESULT_PATH;
				$param_job{_PATH2_}       = "$user/$machine/$libname/$mat/$cas/$date";
				$param_job{_PATH_}        = "$path/$date";
				$param_job{_FILE_OUT_}    = "$path/$date/$execname.log";
				$param_job{_FILE_ERR_}    = "$path/$date/$execname.err.log";
				$param_job{_DATE_}        = "$date";
				$param_job{_PM2_FLAVOR_}  = $param_exec{_PM2_FLAVOR_};
				$param_job{_NBPROCTOTAL_} = 
				    $machines{$machine}{'nbcores'} * $param_job{_NBPROC_};
				
				$options  = "-t $nbthread -ooc $ooc_limit";
				$options .= " -v $verboselevel";
				$options .= " -incomp $amalg $level" if ($withricar == 1);
				$options .= " -iparm $path/$date/iparm.txt -dparm $path/$date/dparm.txt" if ($execname =~ /_param/); 
				$param_job{_EXECUTABLE_}  = "./$execname $options";
				
				$param_job{_EXECCMD_}     = "";
				if (!( $param_job{_PM2_FLAVOR_} eq "" ))
				{
				    $param_job{_EXECCMD_}.= "PM2_FLAVOR=$param_job{_PM2_FLAVOR_} ";
				}
				if ( $libname =~ /.*_ooc.*/ )
				{
				    $param_job{_EXECCMD_}     = "ls /tmp/pastix > /dev/null || mkdir /tmp/pastix && ";
				}
				$param_job{_DRIVER_}      = $driver;
				$param_job{_LINK_}        = $link;
				$param_job{_MATRICE_}     = $matrice;
				$param_job{_MAT_}         = $mat;
				$param_job{_CAS_}         = $cas;
				
				$param_job{_PROTOCOL_}    = "tcp";
				$param_job{_NMADSTRAT_}   = "default";
				if ( $libname =~ /ComOrd/ )
				{
				    $param_job{_NMADSTRAT_} = "aggreg";
				}
				if (defined($opts{p}))
				{
				    $param_job{_PROTOCOL_}  = "$opts{p}";
				}

				if ( getiparm("IPARM_IO_STRATEGY", @iparm) == 2 )
				{
				    $param_job{_POSTEXEC_}  = "mkdir $matwoext;\n";
				    $param_job{_POSTEXEC_} .= "cp graphgen $matwoext/graphname\n";
				    $param_job{_POSTEXEC_} .= "cp ordergen $matwoext/ordername\n";
				    $param_job{_POSTEXEC_} .= "cp symbgen  $matwoext/symbmtx\n";
				}
				
				# Execution du job
				if ( defined $machines{$machine} && defined $machines{$machine}{'submit'} )
				{
				    # Creation du fichier de job
				    $param_job{_TIME_}        = $machines{$machine}{'time'};
				    if ( $multiprocbynode ){
					$param_job{_TASKPERNODE_} = 
					    $machines{$machine}{'nbcores'} / $nbthread 
					    if ( $nbthread < 16 );
					$param_job{_NBNODE_}      = 
					    ceil($nbproc / $param_job{_TASKPERNODE_});
					$param_job{_NBPROCTOTAL_} = 
					    $machines{$machine}{'nbcores'} * $param_job{_NBNODE_};
				    }
				    $param_job{_MEM_}         = $param_job{_NBNODE_} *  
					$machines{$machine}{'mempernode'};
				    $param_job{_EXECCMD_}     = $machines{$machine}{'execcmd'};
				    
				    # On force a passer en parallele
				    $param_job{_NBTHREAD_}    = $machines{$machine}{'nbcores'};
				    
				    MakeTemplate($TPL_ROOT."/".$machines{$machine}{'template'},
						 "$path/$date/".$execname.".".$machines{$machine}{'script'},
						 %param_job);
				    
				    system("chmod +x $path/$date/".$execname.".".$machines{$machine}{'script'});
				    if ( $besteffort ) 
				    {
					print M $machines{$machine}{'submit'}.
					    " $path/$date/".$execname.".".
					    $machines{$machine}{'script'}.
					    " ".$machines{$machine}{'argsbe'}."\n";
				    }
				    else
				    {
					print M $machines{$machine}{'submit'}.
					    " $path/$date/".$execname.".".
					    $machines{$machine}{'script'}.
					    " ".$machines{$machine}{'args'}."\n";
					
				    }
				}
				else 
				{ 
				    # Mode interactif  
				    my $ret;
				    my $ligne;
				    my $sortie  = "2>> $param_job{_FILE_ERR_}  >> $param_job{_FILE_OUT_}";
				    
				    chdir $path/$date;
				  seq:
				    $deadlock = 0;
				    eval {
					local $SIG{ALRM} = sub { die "deadlock\n" }; # NB: \n required
					$timeout = $timeout * 5 if ($withricar == 1);
					alarm $timeout;
					$ret = system("rm -f $param_job{_FILE_ERR_}  $param_job{_FILE_OUT_}");
					$ret = system("echo \"cd $path/$date &&".
						      "$param_job{_EXECCMD_} mpirun -np $nbproc ".
						      "$param_job{_EXECUTABLE_} ".
						      "$driver\" $sortie");
					$ret = system("cd $path/$date &&".
						      "$param_job{_EXECCMD_} mpirun -np $nbproc ".
						      "$param_job{_EXECUTABLE_} ".
						      "$driver $sortie");
					alarm 0;
				    };
				    if ($@) {
					die unless $@ eq "deadlock\n";   # propagate unexpected errors
					# timed out
					print color 'red';
					system ("mpirun -np $nbproc pkill $execname");
					print "Deadlock occured\n";
					print color 'reset';
					$deadlock = 1;
				    }
				    else {
				    }
				    $exec_tot = $exec_tot +1;
				    my $path_exec = "$RESULT_PATH/$user/$machine/$libname/";

				    my $patch_name = `ls -altr $path_exec/patch* |tail -n 1`;
				    $patch_name =~ s/.*(patch_.*)/$1/g;
				    my $datecomp = $patch_name;
				    $datecomp =~ s/patch_//;
				    my $patchmd5  = "";
				    $patchmd5 = `md5sum $path_exec/$patch_name`;
				    $patchmd5 =~ s/([^ ]*) .*/$1/;
				    chomp $patchmd5;
				    
				    my $configmd5 = `md5sum $path_exec/config.in`; 
				    $configmd5 =~ s/([^ ]*) .*/$1/;
				    chomp $configmd5;
				    
				    if ($databasemode)
				    {
					my $branche = `basename $path_src`;
					chomp $branche;
					my $execname_base = `basename $execname`;
					chop $execname_base;
					my $log_po = "$path_exec/log_po_$datecomp";
					chomp $log_po;
					my %metriques = LogGetMetriques($log_po, "");
					
					my $version = $metriques{"vers"};
					addExecution($param_job{_FILE_OUT_}, $param_job{_FILE_ERR_},
						     "$path/$date/iparm.txt", "$path/$date/dparm.txt",
						     $execname_base, $branche, $patchmd5, $configmd5, $user, $version); 
				    }
				    
				    if ( $verbosemode )
				    {
					
					my ($str, %metriques) = LogPrintMetriques("$param_job{_FILE_OUT_}", 
										  "$param_job{_FILE_ERR_}", 
										  @lstmetr);
					
					$str    = " ".$str;
					my $stop = 0;					
					if ( ($ret != 0)
					     || ( (!($with2d)) && 
						  (( $metriques{'iter'} == 250 )  ||
						   ( $metriques{'fact'} == -1.0 ) ||
						   ( $metriques{'iter'} == -1 )))
					     || ( $with2d &&
						  ( $metriques{'fact'} == -1.0 )))
					{
					    print color 'red';
					    print $str."\n";
					    print color 'reset';
					    $stop = 1;
					}
					elsif ( (!($with2d)) && ($metriques{'iter'} > 2))
					{
					    print color 'yellow';
					    print $str."\n";
					    print color 'reset';
					    $stop = 1;
					}
					else {
					    print $str."\n";
					}
					while ( $stop )
					{
					    print '(q : quit, c : continue, r : repeat, o : show output, e : show err) : ';
					    $ligne = <STDIN>;
					    chop $ligne;
					    if ( $ligne eq "q" )
					    {
						exit;
					    } 
					    elsif ( $ligne eq "c" )
					    {
						$stop = 0;
						next;
					    }
					    elsif ( $ligne eq "r" )
					    {
						$stop = 0;
						goto seq;
						next;
					    }
					    elsif ( $ligne eq "o" )
					    {
						system("less $param_job{_FILE_OUT_}");
						next;
					    }
					    elsif ( $ligne eq "e" )
					    {
						system("less $param_job{_FILE_ERR_}");
						next;
					    } 
					    else
					    {
						next;
					    }
					}
				    }
				    my ($str2, %metriques) = LogPrintMetriquesHtml("$param_job{_FILE_OUT_}", 
										   "$param_job{_FILE_ERR_}",
										   @lstmetr);
				    
				    
				    if ($deadlock == 1)
				    {
					$mailtxt .= "<tr><font color=red>Deadlocked</font></tr>\n";
				    }
				    if ( ($ret != 0)
					 || ($metriques{'mend'} != 0)
					 || ( (!($with2d)) && 
					      (( $metriques{'iter'} == 250 )  ||
					       ( $metriques{'fact'} == -1.0 ) ||
					       ( $metriques{'iter'} == -1 )))
					 || ( $with2d &&
					      ( $metriques{'fact'} == -1.0 )))
				    {
					$str2     =~ s/<td>/<td><font color=red>/g;
					$str2     =~ s/<\/td>/<\/font><\/td>/g;
					$mailtxt .= "<tr><td></td>$str2</tr>\n";
					$exec_err ++;
				    }
				    elsif ((!($with2d)) && ($metriques{'iter'} > 2))
				    {
					$str2     =~ s/<td>/<td><font color=orange>/g;
					$str2     =~ s/<\/td>/<\/font><\/td>/g;
					$mailtxt .= "<tr><td></td>$str2</tr>\n";
					$exec_war ++;
				    }
				    else
				    {
					$mailtxt .= "<tr><td></td>$str2</tr>\n";
				    }
				}
			    }  # Fin boucle executables
			} # Fin Boucle thread
		    } #fin boucle proc
		    print "\n";
		    $mailtxt .= "</table>\n";
		} # Fin Boucle Paramtres
	    }# Fin Boucle Matrice 
	    print "\n";
	    $done ++;
	} #Fin de cas 
	
    } # Fin Boucle liste de cas
    
    chdir $rep;
    close(M);
}

#
# Function: GenLatex
#
# Generate Latex output
sub GenLatex {

    my $user = `whoami`;
    chomp $user;
    my $libname;
    my $matrice;
    my $j;
    my $total;
    my $libnum = 0;
    
    # Sauvegarde du rep courant
    my $rep = `pwd`; chop $rep;
			
    open(F, ">Table.tex");
    PrintEntete(*F);

    # Boucle sur les cas
    foreach my $run_ref (@liste_runs)
    {
	my %run = %$run_ref;
	my $list_lib_ref = $run{'libnames'};
	my @list_lib = @$list_lib_ref; 
	$total = $total + $#list_lib+1; 
	
    }
    print F '\section{Paramètres utilisée}'."\n\n";
    # Boucle sur les parametres
    foreach my $paramlist (keys %liste_cas)
    {
	print F '\subsection{'.latexifyword($paramlist).'}'."\n\n";
	my %liste_param = %{$liste_cas{$paramlist}};
	print F '\begin{tabular}{|l|c|}'."\n";
	print F '\hline'."\n";
	if (defined($liste_param{'iparm'}))
	{
	    my %liste_iparm = %{$liste_param{'iparm'}};
	    foreach my $liparm (keys %liste_iparm)
	    {
		print F latexifyword($liparm)." & ".$liste_iparm{$liparm}.'\\\\'."\n";
	    }
	}
	
	if (defined($liste_param{'dparm'}))
	{
	    my %liste_dparm = %{$liste_param{'dparm'}};
	    foreach my $ldparm (keys %liste_dparm)
	    {
		print F latexifyword($ldparm)." & ".$liste_dparm{$ldparm}.'\\\\'."\n";
	    }
	}
	print F '\hline'."\n";
	print F '\end{tabular}'."\n";
    }

    # Boucle sur les cas
    foreach my $run_ref (@liste_runs)
    {	    
	my %run = %$run_ref;
	my $liste_matrices_ref = $run{'matrices'};
	my @liste_matrices = @$liste_matrices_ref;
	
	# Boucles sur les Matrices
	for (my $j=0; $j <= $#liste_matrices; $j++)
	{
	    my $matrice = $liste_matrices[$j];
	    my $mat = `basename $matrice`;
	    chop $mat;
	    
	    # Suppression du suffixe 
	    $mat =~ s/\..*//;
	    print F '\section{Matrice '.$mat.'}'."\n\n";
	    print " ------ $mat \n";
	
	    my $lib_list_ref = $run{'libnames'};
	    foreach my $libname (@$lib_list_ref)
	    {
		print " ---- $libname \n";
		
		# Boucle sur les executables et leur parametres
		my $liste_exec_ref = $run{'executables'};
		my @liste_exec     = @$liste_exec_ref;
		my $i = -1;
		foreach my $execname (@liste_exec)
		{
		    $i = $i + 1;
		    $execname = `basename $execname`;
		    chomp $execname;
			
			print F '\subsection{Executable : '.latexifyword($libname).' - '.$execname."}\n";
			# Boucle sur les jeu de parametres
			my $liste_exec_param_ref = $run{'parametres'};
			my @liste_exec_param     = @$liste_exec_param_ref;
			foreach my $param_name (@liste_exec_param)
			{
			    
			    print F '\subsubsection{Cas : '. latexifyword($param_name)."}\n";
			    my $level = 1;
			    
			    my %liste_param = %{$liste_cas{$param_name}};

			    if ( $libname =~ /.*_metis.*/ ) 
			    {
				if (defined($liste_param{'iparm'}))
				{
				    my %liste_iparm = %{$liste_param{'iparm'}};
				    foreach my $liparm (keys %liste_iparm)
				    {
					if ($liparm =~ /IPARM_LEVEL_OF_FILL/)
					{
					    $level = $liste_param{'iparm'}{$liparm};
					}
				    }
				}
			    }
			    
			    if ( ( $libname =~ /.*_metis.*/ ) && ( $level!=-1 ) )
			    {
				print "METIS sans KASS n'est pas possible.\n";
			    }
			    else
			    {
				
				#on liste les dossiers d'execution
				my $cas = sprintf("%02dP_%02dT_%s", 
						  $liste_param{'MIN_PROC'}, $liste_param{'MIN_THREAD'}, $param_name);

				my @listedate = `ls $RESULT_PATH/$user/$machine/$libname/$mat/$cas/ | grep -E "[0-9]{8}-[0-9]{4}"`;
				foreach $date (@listedate)
				{
				    chomp $date;
				    print F '\paragraph{Date : '.$date."}\n\n";
				    print F '\begin{tabular}{|l';
				    my $nbthread = $liste_param{'MIN_THREAD'};
				    for ( my  $nbproc = $liste_param{'MIN_PROC'};
					  $nbproc < $liste_param{'MAX_PROC'}+1;
					  $nbproc = $nbproc * 2){
					next if ($nbthread*$nbproc > $liste_param{'MAX_PROCTOTAL'});
					print F '|c';
				    }
				    print F '|}'."\n";
				    print F '\hline'."\n";
				    print F '\backslashbox{Thread Number}{MPI Processus number}'."\n";
				    for ( my $nbproc = $liste_param{'MIN_PROC'};  
					  $nbproc < $liste_param{'MAX_PROC'}+1;
					  $nbproc = $nbproc * 2){
					next if ($nbthread*$nbproc > $liste_param{'MAX_PROCTOTAL'});
					print F '& '.$nbproc;
				    }
				    print F '\\\\'."\n";
				    print F '\hline'."\n";
				    for ( $nbthread = $liste_param{'MIN_THREAD'};  
					  $nbthread < $liste_param{'MAX_THREAD'}+1;
					  $nbthread = $nbthread * 2){
					
					# Boucle sur les cas en thread et proc
					print F $nbthread;
					for ( my $nbproc = $liste_param{'MIN_PROC'};  
					      $nbproc < $liste_param{'MAX_PROC'}+1;
					      $nbproc = $nbproc * 2){
					    my $cas = sprintf("%02dP_%02dT_%s", 
							      $nbproc, $nbthread, $param_name);
			
					    # On liste les fichiers de log
					    my $log = "$RESULT_PATH/$user/$machine/$libname/$mat/$cas/$date/$execname.log";
					    my $logerr = "$RESULT_PATH/$user/$machine/$libname/$mat/$cas/$date/$execname.err.log";
					    if (! -f $log)
					    {
						#print F ' & \textcolor{'.$colors[$i%17].'}{-} ';
						print F ' & - ';
						next;
					    }
					    my %metriques = LogGetMetriques($log, $logerr);
					    #print F ' & \textcolor{'.$colors[$i%17].'}{'.$metriques{'fact'}.'} ';
					    print F ' & '.$metriques{'fact'}.' ';
					} # Fin Boucle sur les procs
					
					print F '\\\\'."\n";
					print F '\hline'."\n";
				    } # Fin Boucle sur les threads
				    print F '\end{tabular}'."\n";
				} # Fin de boucle sur les dates
			    } # Fin else
			} # Fin param list
		    } #Fin execname
		} #Fin libname
	    print "\n";
	} # Fin boucle sur les matrices
    } # Fin Boucle sur les liste de runs
	    
    print F '\end{document}'."\n";
    chdir $rep;
}

if ( $machine =~ /hagrid/ ){
    `export LD_LIBRARY_PATH=/net/cremi/mfaverge/acml3.6.0/gnu64/lib:`;
}

print " -- Execution sur $machine -- \n";

getopts("hvHAcfegbo:p:mds:V", \%opts);

if ( defined($opts{h}) ){
    Usage();
    exit;
}
if ( defined($opts{H}) ){
    PrintOption();
    exit;
}
if ( defined($opts{v}) ){
    $verbosemode = 1;
}
if ( defined($opts{m}) ){
    $mailingmode = 1;
}
if ( defined($opts{d}) ){
    $databasemode = 1;
}

PrintOption();
if ( defined($opts{A}) || defined($opts{c}) ){
    print " -- Compilation --\n";
    Compilation();
}
if ( defined($opts{A}) || defined($opts{f}) ){
    print " -- Create Dirs --\n";
    CreateDirAndFiles();
}
if ( defined($opts{A}) || defined($opts{e}) ){
    my $file;

    print " --  Execution  --\n";
    $file = "launch";
    $file = $opts{o} if (defined($opts{o}));
    $besteffort = 1 if (defined($opts{b}));
    ExecuteJobs($file);

}
if ($mailingmode)
{
    open(SENDMAIL, "|$sendmail") or die "Cannot open $sendmail: $!";
    print SENDMAIL "Reply-to: do_not_reply\@this_is_not_a_mail.com\n";
    my $comp_ok = $comp_tot - $comp_err;
    my $exec_ok = $exec_tot - $exec_err;
    my $version = `svnversion $PASTIX_ROOT`;
    chomp $version;
    print SENDMAIL "Subject: [Pastix-run] r$version : $comp_ok/$comp_tot compilations OK, $exec_ok/$exec_tot executions OK ($exec_war Warnings)\n";
    print SENDMAIL "To: $USERMAIL\n";
    print SENDMAIL "Content-type: text/html\n\n";
    print SENDMAIL $mailtxt;
    close(SENDMAIL);
}

if ( defined($opts{g}) ){
     GenLatex();
}
