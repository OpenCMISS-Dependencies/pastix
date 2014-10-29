#!/usr/bin/perl

###############################################################################
#
#   Script: GetAllRes.pl
#
###############################################################################
#
#   This script parse all logs and print summary in the console
#
#   Usage:
#   >  GetAllRes.pl [ -hebr ][ -d=Directory ] [ -g=patern ] [-c=metr1,metr2,...]
#   >     -h   Affiche cette aide
#   >     -p   Affiche les résultats en distibue
#   >     -b   Affiche uniquement les meilleurs temps par cas
#   >     -r   Supprie le chemin du repertoire de resultats dans l'affichage
#   >     -e   N'affiche pas les resultats a -1
#   >     -x   N'affiche que les resultats a -1
#   >     -d   Choose directory result
#   >     -g   limits result following patern
#   >     -f   format affichage
#   >     -c   desactivation de la couleur
#   >     -m   choose metrics among the following possibilities :
#   >                    mach   : Cluster name
#   >                    matr   : Matrice name
#   >                    exec   : Executable name
#   >                    vers   : SVN Revision number
#   >                    day    : Execution date
#   >                    hour   : Execution hour
#   >                    fact   : Numerical factorization time
#   >                    solv   : Solve time
#   >                    refin  : Reffinment time
#   >                    iter   : Number of reffinment iteration
#   >                    norm   : Precision
#   >                    nproc  : Processor number
#   >                    nthrd  : Thread number
#   >                    nbtask : Total number of thread (Nbproc*Nbthread)
#   >                    mem    : Maximum memory used
#   >                    err    : Number of errors
#   >                    warn   : Number of warning
#   >                    sym    : Matrix symmetry
#   >                    size   : Matrix size N
#   >                    nnza   : Matrix number of nonzeros (Only lower part in symmetric graph)
#   >                    nnzl   : Number of nonzeros of L
#   >                    opc    : Number of OPC (LLt or LU)
#   >                    esp    : Number of tasks added by the esp option
#   >                    tord   : Ordering Time
#   >                    tana   : Analyse Time
#   >                    fillin : Fill-in
#   >                    pfact  : Prediction time for the numerical fatorization
#   >                    smsize : Size of the structure SolverMatrix
#   >                    csize  : Size of coeftab
#   >                    ntsk   : Number of tasks
#
#   Authors:
#     Mathieu Faverge - faverge@labri.fr
#     Xavier  Lacoste - lacoste@labri.fr
#
###############################################################################

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";
use Getopt::Std;
use Term::ANSIColor;
use Modules::Common;
use Modules::LogParser;

my $formatres      = "(%7.2f, %5.2f, %2d, %4.2f, %4.2f, %3d)";
my $format         = "%-40s %-10s %-23s";
my @units          = ( "Ko", "Mo", "Go", "To", "Po");
my $showerr        = 1;
my $showcolor      = 1;
my $showresultpath = 1;
my $showbest       = 0;
my $formataff      = 0;
my $GREP           = "";
my $prefixdist     = "";
my @listemetriques;

#
# Function: GetResultats
#
# Parse all log files and print information from it into the console.
#
sub GetResultats #()
{

    my $executable;
    my $matrice;
    my $str;
    my %metriques;
    my %old_metriques;
    my $filedate;
    my $cas;
    my $res;
    my $old_log;

    # Sauvegarde du rep courant
    my $rep = `pwd`; chop $rep;

    # Boucle sur les cas
    # On liste les fichiers de log
    my @listelog;
    if (!(`uname` eq "AIX\n"))
    {
         @listelog = `find $RESULT_PATH  -regex '.*/.*\.log' | grep -E '[0-9]{8}-[0-9]{4}/.*\.log' | grep -v "\.err\.log" $GREP | sort`;
    }
    else
    {
         @listelog = `find $RESULT_PATH  -name *.log | grep -E "[0-9]{8}-[0-9]{4}/.*\.log" | grep -v "\.err\.log" $GREP | sort`;
    }
    if(@listelog == 0)
    {
	print "No log file\n";
	exit;
    }

    my $old_cas = "";

    #
    # Affichage de l'entete des colonnes
    #
    if ( !$formataff )
    {
	$str = sprintf("%-95s", "log");
    }

    $str .= LogPrintMetriquesHeader(@listemetriques); 

    #foreach my $metr (@listemetriques)
    #{
    #	$str =sprintf("%s %-9s",$str, $metr)
    #}
    print $str."\n";
    
    #
    # Boucle sur les fichiers de log
    #
    foreach my $log (@listelog)
    {
	
	# Suppresion du retour a la ligne
	chop $log;
	my $logerr = $log;
	$logerr =~ s/\.log/.err.log/; 

	#Recuperation de la date
	$filedate = $log;
	$filedate =~ s/\/([0-9]{8}-[0-9]{4})\/.*\.log/$1/;

#	$data{$executable}{$mat}{$cas}{$filedate} = ();
	
	# Si ce n'est pas un fichier on le passe
	if ( ! -f $log || ! ($log =~/.*\.log/))
	{
	    print "$log\n";
# 	    $data{$executable}{$mat}{$cas}{$filedate}{Facto} = "x";
# 	    $data{$executable}{$mat}{$cas}{$filedate}{Solve} = "x";
# 	    $data{$executable}{$mat}{$cas}{$filedate}{Mem} = "x";
	    next; 
	}
	 
	%metriques = LogGetMetriques($log, $logerr);

	$log = $metriques{'cas'} if ( ! $showresultpath );

	#suppression du $HOME dans le log
	$log =~ s/$ENV{ HOME}/~/;

	if ( $showbest )
	{
	    if ( $cas eq $old_cas)
	    {
		if ( ($old_metriques{'fact'} == -1.0) || 
		     (($old_metriques{'fact'} > $metriques{'fact'}) && 
		      ($metriques{'fact'} != -1.0 )))
		{

		    foreach my $metr (@listemetriques)
		    {
			$old_metriques{$metr}  = $metriques{$metr};

		    }
		}
	    }
	    else
	    {
		#$res = sprintf($formatres, $old_fact, $old_solv, $old_refin, $old_mem, $old_overheadm, $old_fillrate2, $old_iter);
		if (!$formataff)
		{
		    $str = sprintf("%-95s", $log);
		}
		
		foreach my $metr (@listemetriques)
		{
		    $str = sprintf("%s $metriqueconf{$metr}{'fmt'}", $str, $old_metriques{$metr});
		}
		foreach my $metr (@listemetriques)
		{
		    $str = sprintf("%s %-8s", $str,$old_metriques{$metr});
		}

		if ($old_metriques{'iter'} == -1 && $showcolor)
		{
		    print color 'red';
		}
		if ($old_metriques{'iter'} == 250 && $showcolor)
		{
		    print color 'yellow';
		}
		print $str."\n" if ( ($showerr != 2 || ($old_metriques{'fact'} = -1.0)) && 
				     ($showerr || ($old_metriques{'fact'}!= -1.0)) && 
				     $old_cas ne "" && $showcolor);
		print color 'reset';
		foreach my $metr (@listemetriques)
		{
		    $old_metriques{$metr}  = $metriques{$metr};
		    
		}
		$old_cas = $cas;
		$old_log = $log;
	    }
	}
	else
	{
	    $str = " ";
	    if (!$formataff)
	    {
		$str = sprintf("%-95s", $log);
	    }

	    foreach my $metr (@listemetriques)
	    {
		$str = sprintf("%s $metriqueconf{$metr}{'fmt'}", $str, $metriques{$metr});
	    }

	    if ($metriques{'iter'} == -1  && $showcolor)
	    {
		print color 'red';
	    }
	    if ($metriques{'iter'} == 250 && $showcolor)
	    {
		print color 'yellow';
	    }
	    print $str."\n" if ( $showerr == 1 || 
				 ($showerr == 2 && ($metriques{'iter'} == -1)) || 
				 ($showerr == 0 && !($metriques{'iter'} == -1)) );
	    print color 'reset' if ( $showcolor);
	}	

    } # Fin boucle sur les fichiers

    if ( $showbest )
    {
	$str = "";
	if (!$formataff)
	{
	    $str = sprintf("%-95s", $old_log);
	}
	
	foreach my $metr (@listemetriques)
	{
	    $str = sprintf("%s $metriqueconf{$metr}{'fmt'}", $str, $old_metriques{$metr});
	}

	if ($old_metriques{'iter'} == -1 && $showcolor)
	{
	    print color 'red';
	}
	if ($metriques{'iter'} == 250    && $showcolor)
	{
	    print color 'yellow';
	}
	print $str if ( ($showerr != 2 || ($old_metriques{'fact'} = -1.0)) &&
			( $showerr || (!($old_metriques{'fact'} == -1.0))));
	if ($showcolor)
	{
	    print color 'reset';
	}
    }	

    chdir $rep;
}

#
# Function: Usage
#
# Prints usage
sub Usage #()
{

    print "GetAllRes.pl [ -hebr ][ -d=Directory ] [ -g=patern ] [-c=metr1,metr2,...]\n";
    print "   -h   Affiche cette aide\n";
    print "   -p   Affiche les résultats en distibue\n";
    print "   -b   Affiche uniquement les meilleurs temps par cas\n";
    print "   -r   Supprie le chemin du repertoire de resultats dans l'affichage\n";
    print "   -e   N'affiche pas les resultats a -1\n";
    print "   -x   N'affiche que les resultats a -1\n";
    print "   -d   Choose directory result\n";
    print "   -g   limits result following patern\n";
    print "   -f   format affichage\n";
    print "   -c   desactivation de la couleur\n";
    print "   -m   choose metrics among the following possibilities :\n";
    LogPrintMetriquesHelp();

}

my %opts;

getopts("hxebfrpd:g:m:c",\%opts);

if ( defined $opts{h} ){
    Usage();
    exit;
}

$showcolor      = 0 if ( defined $opts{c} );
$showerr        = 0 if ( defined $opts{e} );
$showerr        = 2 if ( defined $opts{x} );
$showbest       = 1 if ( defined $opts{b} );
$formataff      = 1 if ( defined $opts{f} );
$showresultpath = 0 if ( defined $opts{r} || defined $opts{f} );

$RESULT_PATH = $opts{d}    if ( defined $opts{d} );
$GREP = "| grep ".$opts{g} if ( defined $opts{g} );
$prefixdist = "d"          if ( defined $opts{p} );

if (defined $opts{m}){
    my $i = 0;
    if ($formataff)
    {
	$listemetriques[$i] =  'exec';
	$i++;
	$listemetriques[$i] =  'exemple';
	$i++;
	$listemetriques[$i] =  'matr';
	$i++;
	$listemetriques[$i] =  'cas';
	$i++;
    }

#     if ($opts{m} =~ /.*fact.*/)     {$listemetriques[$i] = 'fact';  $i++;}
#     if ($opts{m} =~ /.*solv.*/)     {$listemetriques[$i] = 'solv';  $i++;}
#     if ($opts{m} =~ /.*refin.*/)    {$listemetriques[$i] = 'refin'; $i++;}
#     if ($opts{m} =~ /.*ovh[^m^l]*/) {$listemetriques[$i] = 'ovh';   $i++;}
#     if ($opts{m} =~ /.*ovhm.*/)     {$listemetriques[$i] = 'ovhm';  $i++;}
#     if ($opts{m} =~ /.*ovhl.*/)     {$listemetriques[$i] = 'ovhl';  $i++;}
#     if ($opts{m} =~ /.*mem.*/)      {$listemetriques[$i] = 'mem';   $i++;}
#     if ($opts{m} =~ /.*mend.*/)     {$listemetriques[$i] = 'mend';  $i++;}
#     if ($opts{m} =~ /.*fr2.*/)      {$listemetriques[$i] = 'fr2';   $i++;}
#     if ($opts{m} =~ /.*iter.*/)     {$listemetriques[$i] = 'iter';  $i++;}
#     if ($opts{m} =~ /.*norm.*/)     {$listemetriques[$i] = 'norm';  $i++;}
   
    while($opts{m} =~ /([a-zA-Z0-9]*),.*/)
    {
	$listemetriques[$i] = $1;  $i++;
	$opts{m} =~ s/[a-zA-Z0-9]*,(.*)/$1/;
    }
    $listemetriques[$i] = $opts{m};  $i++;
    
}  
else
{
    if ($formataff)
    {
	@listemetriques = ('exec', 'exemple', 'matr', 'cas', 'fact', 'solv', 'refin', 'mem', 'iter'); 
    }
    else
    {
	@listemetriques = ('fact', 'solv', 'refin', 'mem', 'mend', 'iter'); 
    }
}
GetResultats();

