#!/usr/bin/perl

use POSIX;
use Modules::Common;

sub PrintEntete {
    my ($F) = @_ ;
    print $F 
	'\documentclass[10pt,a4wide]{report}'."\n".
	'\usepackage[latin1]{inputenc}'."\n".
	'\usepackage{slashbox}'."\n".
	'\usepackage{a4wide}'."\n".
	'\usepackage{amsmath}'."\n".
	'\usepackage{amsfonts}'."\n".
	'\usepackage{amssymb}'."\n".
	'\usepackage{verbatim}'."\n".
	'\usepackage{multirow}'."\n".
	'\usepackage{color}'."\n".
	'\definecolor{violet}{rgb}{0.5,0,0.5}'."\n".
	'\definecolor{gris25}{gray}{0.75}'."\n".
	'\definecolor{rose}{rgb}{1,0,1}'."\n".
	'\definecolor{bordeaux}{rgb}{0.6,0,0.2}'."\n".
	'\definecolor{turquoise}{rgb}{0.2,1,1}'."\n".
	'\definecolor{ciel}{rgb}{0.4,0.6,0.8}'."\n".
	'\definecolor{mer}{rgb}{0,0,0.4}'."\n".
	'\definecolor{orange}{rgb}{1,0.8,0}'."\n".
	'\definecolor{vertfonce}{rgb}{0,0.4,0.2}'."\n".
	'\definecolor{violetpastel}{rgb}{0.8,0.6,1}'."\n".
	'\title{Experiments}'."\n".
	'\begin{document}'."\n";
#    '\maketitle'."\n";
}

sub PrintFoot {
    my ($F) = @_ ;
    print $F '\end{document}'."\n";
}


sub PrintTabHeadProc {
    my ($F) = @_ ;
    
    if ($MAX_PROC < $MIN_PROC){
	return 0;
    }
    
    print $F '\begin{tabular}{|c||';
    for ($i=$MIN_PROC; $i <  $MAX_PROC+1; $i*=2){
	foreach $level (@levels_of_fill) {
	    foreach $amalg (@amalgamations) {
		print $F 'c|';
	    }
	}
    }
    print $F '} \hline'."\n ";
    print $F "Proc number ";
    for ($i=$MIN_PROC; $i < $MAX_PROC+1; $i*=2){
	print $F "& \\multicolumn{". ($#amalgamations+1) * ($#levels_of_fill+1)."}{|c|}{$i}";
    }
    print $F '\\\\'."\n";
    print $F "Amalgamation ";
    for ($i=$MIN_PROC; $i <  $MAX_PROC+1; $i*=2){
	foreach $amalg (@amalgamations) {
	    print $F "& \\multicolumn{".($#levels_of_fill+1)."}{|c|}{".$amalg."}";
	}
    }
    print $F '\\\\'."\n";

    print $F "Level of fill";
    for ($i=$MIN_PROC; $i <  $MAX_PROC+1; $i*=2){
	foreach $amalg (@amalgamations) {
	    foreach $level (@levels_of_fill) {
		print $F "& ".$level;
	    }
	}
    }

    print $F '\\\\'."\n".'\hline'."\n".'\hline'."\n";
    
}


sub PrintTabFoot {
    my ($F) = @_ ;
    print $F '\end{tabular}'."\n\n\n";
}

sub GenerateLatex {
    my $executable;
    my $matrice;
    
    open(F, ">Table.tex");
    PrintEntete(*F);
    
    # Boucles sur les Matrices
    for ($j=0; $j <= $#liste_matrices; $j++)
    {
	$matrice = $liste_matrices[$j];
	$mat = `basename $matrice`;
	chop $mat;
	
	# Suppression du suffixe 
	$mat =~ s/\..*//;
	print " ------ $mat \n";
	
	print F '\section{Matrice '.$mat.'}'."\n";
	
	# Boucle sur les dates d'exécutions
	foreach $date (keys %{ $data{$liste_runs[0]}{$mat}{$MIN_PROC}{$MIN_THREAD}{$amalgamations[0]}{$levels_of_fill[0]}})
	{
	    
	    print F '\subsection{Exécution '.$date.'}'."\n";
	    
	    PrintTabHeadProc(*F);
	    
            # Boucle sur le nombre de threads (ligne principales)
	    for ($nbthread = $MIN_THREAD ; $nbthread < $MAX_THREAD+1; $nbthread*=2){
		
		print F '\multirow{'.($#liste_runs+1).'}{*}{'.$nbthread.'}';
		
		# Boucle sur les executables (ligne intermediaires) 
		for($i=0; $i <= $#liste_runs; $i++)
		{
		    $executable = $liste_runs[$i];
		    print " ---- $executable \n";
		    
		    # Boucle sur le nb de proc (colonnes)
		    for ($nbproc=$MIN_PROC; $nbproc < $MAX_PROC+1; $nbproc*=2){
			
			# Boucle sur $amalg
			foreach $amalg (@amalgamations){
			    # Boucle sur $level
			    foreach $level (@levels_of_fill){
				if ( $data{$liste_runs[$i]} 
				     && $data{$liste_runs[$i]}{$mat}
				     && $data{$liste_runs[$i]}{$mat}{$nbproc}
				     && $data{$liste_runs[$i]}{$mat}{$nbproc}{$nbthread}
				     && $data{$liste_runs[$i]}{$mat}{$nbproc}{$nbthread}{$amalg}
				     && $data{$liste_runs[$i]}{$mat}{$nbproc}{$nbthread}{$amalg}{$level}
				     && $data{$liste_runs[$i]}{$mat}{$nbproc}{$nbthread}{$amalg}{$level}{$date} ){
				    print F '& \textcolor{'.$colors[$i%17].'}{'.
					$data{$liste_runs[$i]}{$mat}{$nbproc}{$nbthread}{$amalg}{$level}{$date}{Facto}.'} '; 
				} else {
				    print F '& \textcolor{'.$colors[$i%17].'}{-} ';
				}
			    } # Fin boucle sur $level
			} # Fin de boucle sur amalgamations
		    }# Fin boucle Proc (colonnes)
		    
		    print F '\\\\'."\n";
		    
		} # Fin boucle exec (lignes intermediaires)
		
		print F '\hline'."\n";
		
	    } # Fin boucle threads (lignes principales)
	    
	    PrintTabFoot(*F);

	    print F '\begin{itemize}'."\n";
	    for($i=0; $i <= $#liste_runs; $i++)
	    {
		$toto = $liste_runs[$i];
		$toto =~ s/_/\\_/g;
		print F '\item \textcolor{'.$colors[$i%17].'}{'.$toto.'}'."\n";
	    }
	    print F '\end{itemize}'."\n";

	} # Fin Boucle sur dates 
    } # Fin boucle sur Matrices
    
    PrintFoot(*F);
    close F;
}



return 1;
