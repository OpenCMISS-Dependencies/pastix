#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";
use Getopt::Std;
use Term::ANSIColor;
use Modules::Common;

my $formatres = "(%7.2f, %5.2f, %2d, %4.2f, %4.2f, %3d)";
my $format    = "%-40s %-10s %-23s";
my @units     = ( "Ko", "Mo", "Go", "To", "Po");
my $showerr        = 1;
my $showcolor      = 1;
my $showresultpath = 1;
my $showbest       = 0;
my $formataff      = 0;
my $GREP = "";
my @listemetriques;
my %resultats      = ();
my $matrix;
my $metrique;

###############################################################
#     Récupération des résultats dans les fichiers de log     #
###############################################################

sub GetResultats {

    my $executable;
    my $matrice;
    my $str;
    my %metriques;
    my $filedate;
    my $cas;
    my $res;
    my $old_log;
    my $temp;
    my $numproc;
    my $numthrd;

    # Sauvegarde du rep courant
    my $rep = `pwd`; chop $rep;

    # Boucle sur les cas
    # On liste les fichiers de log
    my @listelog = `find $RESULT_PATH  -name "*log" $GREP | sort`;

    # On quitte s'il n'y a pas de resultats
    exit if(@listelog == 0);

    #
    # Boucle sur les fichiers de log
    #
    foreach my $log (@listelog)
    {
	
	# Suppresion du retour a la ligne
	chop $log;
	
	#Recuperation de la date
	$filedate = $log;
	$filedate =~ s/\/([0-9]{8}-[0-9]{4})[\/.]log/$1/;

	# Si ce n'est pas un fichier on le passe
	next if ( ! -f $log || $log =~/\.log/);
	
	$metriques{'fact'} = -1.0;
	$metriques{'solv'} = -1.0;
	$metriques{'refin'}= -1.0;
	$metriques{'mem'}  = -1;
	$metriques{'mend'} = -1;
	$metriques{'ovh'}  = -1;
	$metriques{'ovhm'} = -1;
	$metriques{'ovhl'} = 950;
	$metriques{'fr2'}  = -1;
	$metriques{'iter'} = -1;
	$metriques{'norm'} = -1.0;

	# Ouverture du fichier de log
	open (M, "$log");
	
	# Parcours du fichier de log
	foreach (<M>)
	{
	    $metriques{'fact'}  = $1 if (/$metriqueconf{'fact' }{'key'}/  && $1 > $metriques{'fact'} );
	    $metriques{'solv'}  = $1 if (/$metriqueconf{'solv' }{'key'}/  && $1 > $metriques{'solv'} );
	    $metriques{'refin'} = $1 if (/$metriqueconf{'refin'}{'key'}/ && $1 > $metriques{'refin'});
	    $metriques{'mem'}   = $1 if (/$metriqueconf{'mem'  }{'key'}/  );
	    $metriques{'mend'}  = $1 if (/$metriqueconf{'mend' }{'key'}/ );
	    $metriques{'ovh'}   = $1 if (/$metriqueconf{'ovh'  }{'key'}/   && $1 > $metriques{'ovh'}  );
	    $metriques{'ovhm'}  = $1 if (/$metriqueconf{'ovhm' }{'key'}/  && $1 > $metriques{'ovhm'} );
	    $metriques{'ovhl'}  = $1 if (/$metriqueconf{'ovhl' }{'key'}/  && $1 < $metriques{'ovhl'} && $1 > 0.5);
	    $metriques{'fr2'}   = $1 if (/$metriqueconf{'fr2'  }{'key'}/   && $1 > $metriques{'fr2'}  && $1 > 0);
	    $metriques{'iter'}  = $1 if (/$metriqueconf{'iter' }{'key'}/ );
	    $metriques{'norm'}  = $1 if (/$metriqueconf{'norm' }{'key'}/ );
	}
	close M;

	# Formatage de l'affichage des floatants
	$metriques{'fact'}  = sprintf("%.2f", $metriques{'fact'});
	$metriques{'solv'}  = sprintf("%.2f", $metriques{'solv'});
	$metriques{'refin'} = sprintf("%.2f", $metriques{'refin'});
	$metriques{'ovh'}   = sprintf("%.2f", $metriques{'ovh'});
	$metriques{'ovhm'}  = sprintf("%.2f", $metriques{'ovhm'});
	$metriques{'ovhl'}  = sprintf("%.2f", $metriques{'ovhl'});
	$metriques{'fr2'}   = sprintf("%.2f", $metriques{'fr2'});
	
	$cas = $log;
	#Suppression du répertoire de stockage
	$cas =~ s/$RESULT_PATH\///;
	# Suppression de la date + [/.]log
	$cas =~ s/\/[0-9]{8}-[0-9]{4}[\/.]log//;

	$log = $cas;

	# Récupération du cas XXP_XXT_... (suppression de tout ce qu'il y a avant le dernier /)
	$temp = $log;	    
	$temp =~ s/[-\/_a-zA-Z0-9]*\///;
	
	#print $temp." ";
	$numproc = $temp;
	$numproc =~ s/([0-9]*)P_.*/$1/;
	$numproc = sprintf("%d", $numproc);
	$numthrd = $temp; 
	$numthrd =~ s/.*_([0-9]*)T.*/$1/;
	$numthrd = sprintf("%d", $numthrd);

	# Récupération du nom de la matrice ( Récupération des infos entre le premier / et le suivant)
	$matrice = $log;
	$matrice =~ s/.*\/([_a-zA-Z0-9]*)\/.*/$1/;
	if ( $matrice ne $matrix )
	{
	    print "Wrong matrice : $matrice $matrix\n";
	    next;
	}
	
        #Récupération du nom de l'exécutable ( Suppression de tt ce qu'il y a après le premier /)
	$executable = $log;
	$executable =~ s/\/.*//;
	
	$resultats{$numproc}                        = () if (! defined($resultats{$numproc}                       ) );
	$resultats{$numproc}{$numthrd}              = () if (! defined($resultats{$numproc}{$numthrd}             ) );
	if (! defined($resultats{$numproc}{$numthrd}{$executable}) )
	{
	    $resultats{$numproc}{$numthrd}{$executable} = () ;
	    $resultats{$numproc}{$numthrd}{$executable}{'fact'} = 9999999999999999.0;
	    $resultats{$numproc}{$numthrd}{$executable}{'solv'} = 9999999999999999.0;
	    $resultats{$numproc}{$numthrd}{$executable}{'refin'}= 9999999999999999.0;
	    $resultats{$numproc}{$numthrd}{$executable}{'mem'}  = -1;
	    $resultats{$numproc}{$numthrd}{$executable}{'mend'} = -1;
	    $resultats{$numproc}{$numthrd}{$executable}{'ovh'}  = 9999999999999999;
	    $resultats{$numproc}{$numthrd}{$executable}{'ovhm'} = 9999999999999999;
	    $resultats{$numproc}{$numthrd}{$executable}{'ovhl'} = -1 ;
	    $resultats{$numproc}{$numthrd}{$executable}{'fr2'}  = 9999999999999999;
	    $resultats{$numproc}{$numthrd}{$executable}{'iter'} = 9999999999999999;
	    $resultats{$numproc}{$numthrd}{$executable}{'norm'} = 9999999999999999;
	}		
	
	$resultats{$numproc}{$numthrd}{$executable}{'fact'}  = $metriques{'fact'}  if ($metriques{'fact'}  < $resultats{$numproc}{$numthrd}{$executable}{'fact'} );
	$resultats{$numproc}{$numthrd}{$executable}{'solv'}  = $metriques{'solv'}  if ($metriques{'solv'}  < $resultats{$numproc}{$numthrd}{$executable}{'solv'} );
	$resultats{$numproc}{$numthrd}{$executable}{'refin'} = $metriques{'refin'} if ($metriques{'refin'} < $resultats{$numproc}{$numthrd}{$executable}{'refin'});
	$resultats{$numproc}{$numthrd}{$executable}{'mem'}   = $metriques{'mem'};
	$resultats{$numproc}{$numthrd}{$executable}{'mend'}  = $metriques{'mend'};
	$resultats{$numproc}{$numthrd}{$executable}{'ovh'}   = $metriques{'ovh'}   if ($metriques{'ovh'}   < $resultats{$numproc}{$numthrd}{$executable}{'ovh'}  );
	$resultats{$numproc}{$numthrd}{$executable}{'ovhm'}  = $metriques{'ovhm'}  if ($metriques{'ovhm'}  < $resultats{$numproc}{$numthrd}{$executable}{'ovhm'} );
	$resultats{$numproc}{$numthrd}{$executable}{'ovhl'}  = $metriques{'ovhl'}  if ($metriques{'ovhl'}  > $resultats{$numproc}{$numthrd}{$executable}{'ovhl'} && $metriques{'ovhl'} > 0.5);
	$resultats{$numproc}{$numthrd}{$executable}{'fr2'}   = $metriques{'fr2'}   if ($metriques{'fr2'}   < $resultats{$numproc}{$numthrd}{$executable}{'fr2'}  && $metriques{'fr2'}  > 0);
	$resultats{$numproc}{$numthrd}{$executable}{'iter'}  = $metriques{'itet'};
	$resultats{$numproc}{$numthrd}{$executable}{'norm'}  = $metriques{'norm'};

	#print $numproc." ".$numthrd." ".$executable." ".$resultats{$numproc}{$numthrd}{$executable}{'fact'}."\n";

    }	

    chdir $rep;
}

sub GenTabLatex
{
    my $line;
    my $nb = $#liste_runs + 1;

    #print $metrique."\n";
    $line = '\begin{tabular}{|c|c|';
    for(my $j=$MIN_PROC; $j<$MAX_PROC+1; $j=$j*2)
    {
	$line .= "c|";
    }
    $line .= "}\n";
    print OUT $line;
    

    # Ligne de titre
    $line = '\hline'."\n"."                   &    ";
    for(my $j=$MIN_PROC; $j<$MAX_PROC+1; $j=$j*2)
    {
	$line .= sprintf(" & %7d ", $j);
    }
    $line .= "\\\\\n";
    print OUT $line;

    for (my $i=$MIN_THREAD; $i<$MAX_THREAD+1; $i=$i*2)
    {
	my $numexec = 0;
	$line = '\hline'."\n".'\multirow'."{$nb}{*}{$i}";
	
	for my $executable (@liste_runs)
	{
	    print "V$numexec = $executable\n";
	    $line .= " & V$numexec ";
	    for(my $j=$MIN_PROC; $j<$MAX_PROC+1; $j=$j*2)
	    {
		#print $j." ".$i." ".$executable." ".$resultats{$j}{$i}{$executable}{$metrique}."\n";
		
		if ( defined($resultats{$j}{$i}{$executable}{$metrique}) )
		{
		    $line .= sprintf(" & %7.2f ",$resultats{$j}{$i}{$executable}{'fact'});
		}
		else
		{
		    $line .= " &    -    ";
		}
	    }
	    $line .= "\\\\\n";
	    print OUT $line;
	    $line = "                  ";
	    $numexec++;
	}
    } 
    print OUT '\hline'."\n".'\end{tabular}'."\n";
}

sub Usage {

    print "GetLatex.pl [ -h][ -m metric ] [ -g=patern ] matrixname\n";
    print "   -h   Affiche cette aide\n";
    print "   -g   limits result following patern\n";
    print "   -d   Choose directory result\n";
    print "   -e   N'affiche pas les resultats a -1\n";
    print "   -f   format affichage\n";
    print "   -c   desactivation de la couleur\n";
    print "   -m   Metrique a conserver (par defaut : fact)\n";
    print "        available metrics : fact : facto time\n";
    print "                            solv : solv  time\n";
    print "                            refin: refin  time\n";
    print "                            mem  : memory used\n";
    print "                            mend : memory used at the end\n";
    print "                            ovh  : overhead total\n";
    print "                            ovhm : overhead max \n";
    print "                            ovhl : overhead local\n";
    print "                            fr2  : fillrate2\n";
    print "                            iter : number of iteration\n";
    print "                            norm : precision\n";

}

my %opts;

getopts("hebfrd:g:m:c",\%opts);

if ( defined $opts{h} || $#ARGV < 0){
    Usage();
    exit;
}

open(OUT, ">test.tex");

$showcolor      = 0 if ( defined $opts{c} );
$showerr        = 0 if ( defined $opts{e} );
$showbest       = 1 if ( defined $opts{b} );
$formataff      = 1 if ( defined $opts{f} );
$showresultpath = 0 if ( defined $opts{r} || defined $opts{f} );

$matrix = $ARGV[0];

$RESULT_PATH = $opts{d}    if ( defined $opts{d} );
$GREP .= "| grep ".$opts{g}." " if ( defined $opts{g} );
$GREP .= "| grep ".$matrix." ";

if (defined $opts{m}){
    if ($opts{m} =~ /.*fact.*/)     {$metrique = 'fact';  next;}
    if ($opts{m} =~ /.*solv.*/)     {$metrique = 'solv';  next;}
    if ($opts{m} =~ /.*refin.*/)    {$metrique = 'refin'; next;}
    if ($opts{m} =~ /.*ovh[^m^l]*/) {$metrique = 'ovh';   next;}
    if ($opts{m} =~ /.*ovhm.*/)     {$metrique = 'ovhm';  next;}
    if ($opts{m} =~ /.*ovhl.*/)     {$metrique = 'ovhl';  next;}
    if ($opts{m} =~ /.*mem.*/)      {$metrique = 'mem';   next;}
    if ($opts{m} =~ /.*mend.*/)     {$metrique = 'mend';  next;}
    if ($opts{m} =~ /.*fr2.*/)      {$metrique = 'fr2';   next;}
    if ($opts{m} =~ /.*iter.*/)     {$metrique = 'iter';  next;}
    if ($opts{m} =~ /.*norm.*/)     {$metrique = 'norm';  next;}
}  
else
{
    $metrique = 'fact'; 
}

print "Script pas mis à jour, essayer plutôt Pastix.pl -g\n";
close M;
#GetResultats();
#GenTabLatex();
