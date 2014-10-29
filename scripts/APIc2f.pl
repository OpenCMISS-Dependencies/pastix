#!/usr/bin/perl

use strict;
use Getopt::Std;

my $increment;
my %opts;
my $file;
my $format;

sub Convert {

    open (APIc, $file); 
    
    my $iparmkey  = "^[ \t]*([ID]PARM_[0-9A-Z_]*)[ \t]*= ([0-9]*),?";
    my $iparmkey2 = "^[ \t]*([ID]PARM_SIZE)[ \t]*= ([0-9]*),?";
    my $comment = "";
    my $chaine; 
    foreach my $line (<APIc> )
    {
	$chaine = $line;
	if (!($comment eq "! ") &&
	    $line =~ /\/\*(.*)/)
	{
	    $comment = "! ";
	    $chaine  = $1;
	}
	if ($comment eq "! " &&
	    $line =~ /(.*)\*\//)
	{
	    $comment = "";
	    $chaine  = $1;
	}
	if ($comment eq "! ")
	{
	    print $comment.$chaine;
	}
	if ($line =~ /$iparmkey/)
	{
	    if ($line =~ /$iparmkey2/)
	    {
		$chaine = sprintf($format, $1, ($2));
	    }
	    else
	    {
		$chaine = sprintf($format, $1, ($2+$increment));
	    }
	    print $comment.$chaine;
	}
	else
	{
	    if ($line =~ /^[ \t]*([0-9A-Z_]*)[ \t]*=[ \t]*([0-9]*),?/)
	    {
		$chaine = sprintf($format, $1, $2);
		print $comment.$chaine
	    }
	    else
	    {
		if ($line =~ /[ ]*#define[ ]*([0-9A-Z_]*)[ \t]*([0-9]*)/)
		{
		    if ( (!($1 =~ /MTX\_IS/)) && (!($1 =~ /API\_H/)) )
		    {
			$chaine = sprintf($format, $1, $2);
			
			print $comment.$chaine;
		    }
		}
		else 
		{
		   if (($line =~ /\_POS\_/) || ($line =~ /Matrix\ type/))
		   {
		       #print "!- ".$1."\n" if (/^\/\**(.*)\*\//);
		       print "\n" if (/^[ ]*$/);
		   }
		}
	    }

	}
	
    }
    close APIc;
}

sub Usage {

    print "APIc2f.pl -f api.h\n";
    print "  convertit le fichier api.h en un header pour le fortran.\n";
    print "   -h   Affiche cette aide\n";
    print "   -m   If using with murge\n";
}


getopts("hmf:",\%opts);

if ( defined $opts{h} || !(defined $opts{f})){
    Usage();
    exit;
}
if (defined $opts{m}){
    $increment = 1;
    $format    = "INTEGER(KIND=MURGE_INTS_KIND), PARAMETER :: %-30s = %d\n";
}
else 
{
    $increment = 1;
    $format    = "#define %-30s INT(%d , KIND=PASTIX_INT_KIND)\n";
}
$file = $opts{f};

Convert();
