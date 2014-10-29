#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

###############################################################
#     Récupération des résultats dans les fichiers de log     #
###############################################################

#Options par défaut
my $branch;
my $DIRNAME;
my $BASENAME;
my $svn    = "scm.gforge.inria.fr/svn/ricar";
my $user   = "hydra";
my $header = "header";
my %Rep2include;
my @file2delete;
my $RELEASE_PATH;
my %opts;
my $NUMREL = "";
    
sub MakeRelease {

    my $numversion;
    my $cmd;

    if (!defined $RELEASE_PATH){
	if (!defined $opts{r})
	{
	    #$cmd = "LANG=C svn info | grep 'Revision' | awk \'{ print $2 }\'";
	    $cmd = "LANG=C svn info".' | grep "Revision:" | awk \'{ print $2 }\'';
	    $numversion = `$cmd`;
	    chop $numversion;
	    if ( $numversion eq "")
	    {
		print "Probleme de recuperation du numero de version\n";
		exit;
	    }
	}
	else
	{
	    $numversion = $opts{r};
	}
	$RELEASE_PATH = $ENV{ PWD}."/pastix_release_".$numversion;
    }

    #Compilation de la doc
    print "Compilation de la documentation\n"; 
    system("PASTIX_SRC=$RELEASE_PATH/$Rep2include{$branch} SCRIPT_PATH=$RELEASE_PATH/$Rep2include{'Scripts'} make -C $RELEASE_PATH/$Rep2include{'Doc/refcard'}");
    system("make -C $RELEASE_PATH/$Rep2include{'Doc/refcard'} clean");
    
    #Suppression des fichiers ne devant pas etre présents dans la release
    foreach my $fichier (@file2delete){
	print "Suppression de $fichier\n";
	system("rm -rf $RELEASE_PATH/$fichier");
    }
 
    # Boucle sur les cas
    # On liste les fichiers source
    my @listefichier = `find $RELEASE_PATH/$Rep2include{$branch}/ -name "*\.[ch]"`;
    if(@listefichier == 0)
    {
	exit;
    }
    
    #Boucle sur les fichiers de source pour ajouter les copyright
    print "Adding header";
    foreach my $fichier (@listefichier){
	
	# Suppresion du retour a la ligne
	chop $fichier;

	# Si ce n'est pas un fichier on le passe
	if ( ! -f $fichier)
	{
	    next; 
	}
	if (system("grep -i copyright $fichier 1>/dev/null")) { 
	    #print "adding " . $header . " to " . $fichier . "\n"; 
	    system("cat $header ".$fichier." > ".$user."_tototititata");
	    system("mv ".$user."_tototititata ".$fichier);
	}
    }

    # Mise a jour pastix_fortran.h
    # WARNING : doit etre genere apres l'ajout des entetes pour etre ecrase
    print ("perl APIc2f.pl -f $RELEASE_PATH/$Rep2include{$branch}/common/src/api.h > $RELEASE_PATH/$Rep2include{$branch}/common/src/pastix_fortran.h\n");
    system("perl APIc2f.pl -f $RELEASE_PATH/$Rep2include{$branch}/common/src/api.h > $RELEASE_PATH/$Rep2include{$branch}/common/src/pastix_fortran.h");
    print ("perl APIc2f.pl -f $RELEASE_PATH/$Rep2include{$branch}/common/src/api.h -m > $RELEASE_PATH/$Rep2include{$branch}/common/src/pastix_fortran.inc\n");
    system("perl APIc2f.pl -f $RELEASE_PATH/$Rep2include{$branch}/common/src/api.h -m > $RELEASE_PATH/$Rep2include{$branch}/common/src/pastix_fortran.inc");

    # Mise a jour iparm dparm
    #system("perl APIc2param.pl $RELEASE_PATH/$Rep2include{$branch}/common/src/api.h");
    #system("mv iparm.txt dparm.txt $RELEASE_PATH/$Rep2include{$branch}/test_c/src/");

    #Ajout du fichier de revision
    print "Ajout du fichier Revison\n"; 
    system("echo $numversion > $RELEASE_PATH/$Rep2include{$branch}/Revision");
    
    #Création de l'archive
    print "Création de l'archive\n";
    $DIRNAME=`dirname $RELEASE_PATH`;
    $BASENAME=`basename $RELEASE_PATH`;
    chop $DIRNAME;
    chop $BASENAME;
    system("(cd $DIRNAME && tar cjf ${BASENAME}.tar.bz2 $BASENAME)");
}

sub Usage {
    
    print "MakeRelease.pl [ -h ][ -d Directory ] [ -f headerfile ] [ -u username ] [ -r numrelease ] [ -b branch ] \n";
    print "   -b   branch";
    print "   -h   Affiche cette aide\n";
    print "   -d   Choose directory for release\n";
    print "   -r   Choose svn release number\n";
    print "   -u   username\n";
    print "   -f   Choose header file\n";

}

getopts("hd:f:u:r:b:",\%opts);
$branch = "trunk";

if ( defined $opts{h} ){
    Usage();
    exit;
}
if (defined $opts{b}){
    $branch = $opts{b};
}
%Rep2include = (
    $branch       => "src",
    "Doc/refcard" =>   "doc/refcard",
    "Doc/fonctions" => "doc/fonctions",
    "Scripts"     => "scripts"
    );
@file2delete  = (
    "$Rep2include{$branch}/blend/oldsrc",
    "$Rep2include{$branch}/sopalin/src/mc64a.F",
    "$Rep2include{$branch}/graph",
    "$Rep2include{$branch}/order/src/order_scotch*",
    "$Rep2include{$branch}/order/src/order_grid*",
    "$Rep2include{$branch}/fax/src/main_*",
    "$Rep2include{$branch}/fax/src/symbol_fax*_grid*",
    "$Rep2include{$branch}/fax/src/symbol_fax*_mesh*",
    "$Rep2include{'Scripts'}/header",
    "$Rep2include{'Scripts'}/fillDB.pl",
    "$Rep2include{'Scripts'}/Pastix.sh",
    "$Rep2include{'Scripts'}/MakeRelease.pl",
    "$Rep2include{'Doc/refcard'}/*.tex",
    "$Rep2include{'Doc/refcard'}/*.dvi",
    "$Rep2include{'Doc/refcard'}/*.sty",
    "$Rep2include{'Doc/refcard'}/Makefile",
    );

if (defined $opts{d}){
    $RELEASE_PATH = $opts{d};
}
if (defined $opts{f}){
    $header = $opts{f};
}
if (defined $opts{u}){
    $user = $opts{u}
}

if (defined $opts{r}){
    $NUMREL = "-r $opts{r}";
}

MakeRelease();
