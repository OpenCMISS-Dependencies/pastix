#!/usr/bin/perl

use FindBin qw($Bin);
use lib "$Bin";
use Modules::Common;
use File::Path;
use File::Copy;


###############################################################
#     Récupération des résultats dans les fichiers de log     #
###############################################################

sub mvResultats {

    my $executable;
    my $matrice;

    # Sauvegarde du rep courant
    $rep = `pwd`; chop $rep;

    @liste_runs = `ls $RESULT_PATH/`;			
    # Boucle sur les cas
    for($i=0; $i <= $#liste_runs; $i++)
    {
	$executable = $liste_runs[$i];
	chop $executable;
	if ( ! -d "$RESULT_PATH/$executable" )
	{
	    next;
	}

#	print " ---- $executable \n";
	
	$data{$executable} = ();

	# Boucles sur les Matrices
	@liste_matrices = `ls $RESULT_PATH/$executable/`;			
	for ($j=0; $j <= $#liste_matrices; $j++)
	{

	    $matrice = $liste_matrices[$j];
	    chop $matrice;
	    if ( ! -d "$RESULT_PATH/$executable/$matrice" )
	    {
		next;
	    }
	    $mat = $matrice;

	    # Suppression du suffixe 
#	    print " ------ $mat \n";

	    $data{$executable}{$mat} = ();

	    # Boucle sur les cas en thread et proc
	    $nbthread = $MIN_THREAD;
	    $nbproc   = $MIN_PROC;
	    @liste_cas = `ls $RESULT_PATH/$executable/$matrice/`;			
	    for ($k=0; $k <= $#liste_cas; $k++)
	    {
		$cas = $liste_cas[$k];
		chop $cas;
		$data{$executable}{$mat}{$cas} = ();
#		print " --------- $cas \n";
		# On liste les fichiers de log
		@listelog = `ls $RESULT_PATH/$executable/$matrice/$cas/*.log 2>/dev/null`;
		if(@listelog == 0)
		{
# 			        print STDERR "$RESULT_PATH/$executable/$mat/$cas/ est vide : ON l'OUBLIE ! \n";
# 			        system "ls $RESULT_PATH/$executable/$mat/$cas/*";
		    next;
		}
		
		mkpath "$ARCHIVE_PATH/$executable/$matrice/$cas/";
		#Boucle sur les fichiers de log
		foreach $log (@listelog){
		    
		    # Suppresion du retour a la ligne
		    chop $log;
				
		    #Recuperation de la date
		    $filedate = `basename $log`; 
		    chop $filedate;
		    
		    if ($filedate =~ /.*$ARGV[0].*/) {
			move($log,$ARCHIVE_PATH."/".$executable."/".$matrice."/".$cas."/".$filedate);
		    }				    
		} # Fin boucle sur les fichiers
	    } # Fin Boucle cas
	} # Fin boucle sur les matrices
    } # Fin Boucle sur les executables

    chdir $rep;
}

$ARCHIVE_PATH=$ENV{HOME}."/Archive".$ARGV[0]; 


# mvResultats();
print "Le script n'a pas été mis à jour.\n";
