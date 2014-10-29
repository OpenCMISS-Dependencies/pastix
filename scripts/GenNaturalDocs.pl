#!/usr/bin/perl

#use strict;
use POSIX;
use FindBin qw($Bin);
use lib "$Bin";
use Getopt::Std;
use Term::ANSIColor;
use Modules::Common;

my $autoyes = 0;

sub GenDoc {

    # Adding C.pm to NaturalDocs
    if ( 0 == system "test -f $NATURAL_DOCS_DIR/Modules/NaturalDocs/Languages/C.pm" )
    {
	print "Fichier ${NATURAL_DOCS_DIR}/Modules/NaturalDocs/Languages/C.pm existant\n";
	diffFiles($PASTIX_NDSCRIPTS,"$NATURAL_DOCS_DIR/Modules/NaturalDocs/Languages/","C.pm");
    }
    else
    {
	print "Fichier ${NATURAL_DOCS_DIR}/Modules/NaturalDocs/Languages/C.pm inexistant\n";
	print "cp $PASTIX_NDSCRIPTS/C.pm $NATURAL_DOCS_DIR/Modules/NaturalDocs/Languages/C.pm\n";
	`cp $PASTIX_NDSCRIPTS/C.pm $NATURAL_DOCS_DIR/Modules/NaturalDocs/Languages/C.pm`;
    }
    
    # Editing Language.txt
    diffFiles($PASTIX_NDSCRIPTS,"$NATURAL_DOCS_DIR/Config","Languages.txt");
    
    #Editing Language.pm
    diffFiles($PASTIX_NDSCRIPTS,"$NATURAL_DOCS_DIR/Modules/NaturalDocs","Languages.pm");

    $cmd  = "$NATURAL_DOCS_DIR/NaturalDocs -i ".$PASTIX_SRC{"T"};
    $cmd .= " -o HTML $PASTIX_DOC -p $PASTIX_PROJECT";
    $cmd .= " -xi ".$PASTIX_SRC{"T"}."/blend/oldsrc";
    my @lsdir = `find ../trunk/ -name obj -type d`;
    for my $dir (@lsdir)
    {
	chomp $dir;
	$cmd .= " -xi ".$dir;
    }
    $cmd .= " -img ".$SCRIPT_PATH."/logo";# -s Pastix -s sh_bright";
    print $cmd;
    #creatind NaturalDocs Directories
    if ( 0 != system "test -d $PASTIX_DOC"  )
    {
	`mkdir $PASTIX_DOC`;
	`mkdir $PASTIX_PROJECT`;
	system "cp ".$PASTIX_NDSCRIPTS."/Menu.txt ".$PASTIX_PROJECT."/Menu.txt ";
	system "cp ".$PASTIX_NDSCRIPTS."/Topics.txt ".$PASTIX_PROJECT."/Topics.txt ";
    }
    else
    {
	#print "Dossiers présents\n";
    }
    diffFiles($PASTIX_NDSCRIPTS,$PASTIX_PROJECT,"Menu.txt");
    diffFiles($PASTIX_NDSCRIPTS,$PASTIX_PROJECT,"Topics.txt");
    #Running Natural Docs
    system "$cmd";
    system "cp ".$PASTIX_NDSCRIPTS."/sh_c.js ".$PASTIX_DOC."/javascript/sh_c.js ";
    system "cp ".$PASTIX_NDSCRIPTS."/sh_main.js ".$PASTIX_DOC."/javascript/sh_main.js ";
    my @listfiles = `find $PASTIX_DOC/files -regex ".*-[ch]\.html"`;
    for my $file (@listfiles)
    {
	if (system "grep displayequation $file > /dev/null" != 0)
	{
	    system "sed -i -e \"s/<head>/<head><script>function displayequation(latex){document.write(\\\"<img src='http:\\/\\/www.codecogs.com\\/eq.latex?\\\" + latex + \\\"' alt='\\\" + latex + \\\"' \\/>\\\");}<\\/script>/g\" $file";

	    system "perl -pe 's/\\\$\\\$\([^\\\$]*\)\\\$\\\$/<script>displayequation\(\"\$1\"\)<\\\/script>/g' -i $file";
	}
	    
    }
	

}


sub diffFiles {
    my ($srcdir, $destdir, $file) = @_;
	
    if (0 != system "diff $destdir/$file $srcdir/$file > /dev/null")
    {
	print "Fichiers $file différents : \n";
	system ("diff -d -u -tBbwE $destdir/$file $srcdir/$file");
	print "1) Edit files manualy\n";
	print "2) cp $srcdir/$file to $destdir/$file \n";	
	if ($autoyes == 1)
	{
	    $answer = "2";
	}
	else
	{
	    $answer = <STDIN>;
	    chomp $answer;
	}
	if (defined($answer) && ($answer=~/1/ || $answer=~/2/))
	{
	    if ($answer=~/1/)
	    {
		system "emacs $destdir/$file $srcdir/$file";
	    }
	    else
	    {
		system "cp $srcdir/$file $destdir/$file";
	    }
			    
	}
	else
	{
	    print "Undefined answer : $answer\n";
	    exit;
	}
    }
    else
    {
#	print  "Fichiers $file identiques\n";
    }
}


sub Usage {

    print "GenNaturalDocs.pl\n";
    print "   -h   Affiche cette aide\n";
    print "   -y   Auto answer questions\n";
    print "Génère la documentation au format html.\n";
}


getopts("hy",\%opts);

if ( defined $opts{h} ){
    Usage();
    exit;
}

if ( defined $opts{y} ){
    $autoyes = 1;
}

GenDoc();
