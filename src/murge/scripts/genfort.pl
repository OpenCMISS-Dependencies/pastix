#!/usr/bin/env perl

###############################################################################
#
#   Script: genfort.pl
#
###############################################################################
#
#   A scrip that converts <murge.h> file into a Fortran include file.
#
#   Usage:
#    > ./genfort.pl -f murge.h
#    >    converts the file murge.h into a fortran header.
#    >     -h          Shows this help
#    >     -c <0,1>    Defines if COEF ar complexes (1) or reals (0)
#    >     -r <size>   Defines COEF/REAL kind (4 or 8 usually)
#    >     -s <size>   Defines INTS kind (4 or 8 usually)
#    >     -l <size>   Defines INTL kind (4 or 8 usually)
#
#   Authors:
#     Xavier Lacoste - lacoste@labri.fr
#
###############################################################################
use POSIX;
use strict;
use Getopt::Std;

###############################################################################
# Group: Variables

#
#   integer: real
#     Kind of the reals/complex
#   
#   integer: ints
#     Kind of the INTS
#
#   integer: intl
#     Kind of the INTL
#
#   integer: complex
#     Indicate if COEFs are complex or real
#
#   string: fichier
#     Path to the <murge.h> file
#
#   string: format
#     Format to print PARAMETERs.
#
#   hash: opts
#     Options given to the script
my $real     = 0;
my $ints     = 0;
my $intl     = 0;
my $complex  = 0; 
my $fichier;
my $format = "INTS, PARAMETER :: %-30s = %d\n";
my %opts;

###############################################################################
# Group: Functions

#
# Function: Usage
#
# Prints usage.
#
sub Usage {

    print "./genfort.pl -f murge.h\n";
    print "  converts the file murge.h into a fortran header.\n";
    print "   -h          Shows this help\n";
    print "   -c <0,1>    Defines if COEF ar complexes (1) or reals (0)\n";
    print "   -r <size>   Defines COEF/REAL kind (4 or 8 usually)\n";
    print "   -s <size>   Defines INTS kind (4 or 8 usually)\n";
    print "   -l <size>   Defines INTL kind (4 or 8 usually)\n";
}
#
# Function: printTab
# 
# Print *chaine* with *tabcount* indentations.
#
# If *comm* is 0, it will also replace INTS, INTL, REAL and COEF by the 
# Correct value.
#
# Parameters:
#   chaine   - String to print
#   tabcount - Number of indentations to add.
#   comm     - Indicate if we are in a comments section.
#
sub printTab # ($chaine, $tabcount, $comm)
  {
    my ($chaine, $tabcount, $comm) = @_;
    for (my $i = 0; $i < $tabcount; $i++)
    {
	$chaine = sprintf("   %s",$chaine);
    }
    if ($comm == 0)
    {
	if ($ints != 0)
	{
	    $chaine =~ s/INTS,/INTEGER(KIND=$ints),/g;
	}
	else
	{
	    $chaine =~ s/INTS,/INTEGER,/g;
	}
	if ($intl != 0)
	{
	    $chaine =~ s/INTL,/INTEGER(KIND=$intl),/g;
	}
	else
	{
	    $chaine =~ s/INTL,/INTEGER,/g;
	}
	if ($complex)
	{
	    if ($real != 0)
	    {
		$chaine =~ s/COEF,/COMPLEX(KIND=$real),/g;
		$chaine =~ s/REAL,/REAL(KIND=$real),   /g;
	    }
	    else
	    {
		$chaine =~ s/COEF,/COMPLEX,/g;
		$chaine =~ s/REAL,/REAL,   /g;
	    }
	}
	else
	{
	    if ($real != 0)
	    {
		$chaine =~ s/COEF,/REAL(KIND=$real),   /g;
		$chaine =~ s/REAL,/REAL(KIND=$real),   /g;
	    }
	    else
	    {
		$chaine =~ s/COEF,/REAL,   /g;
		$chaine =~ s/REAL,/REAL,   /g;
	    }
	}
	$chaine =~ s/MPI_COMM,/MPI_COMM,           /g; 
	$chaine =~ s/CHARACTER\(len=\*\),/CHARACTER(len=*),           /g; 
    }
    print $chaine;

}
#
# Function: Convert 
#
# Main function.
#
# Converts the header <murge.h> file into a Fortran include.
#
sub Convert {

    my $startcom  = 0;
    my $startenum = 0;
    my $countenum;
    my $chaine;
    my $tabcount = 0;
    my $interfaceprinted = 0;
    open (APIc, $fichier); 
    
    foreach my $line (<APIc> )
    {
	if ($startcom == 0)
	{
	    if ($startenum == 0)
	    {
		if ($line =~ /\/\*(.*)/)
		{
		    if ($line =~ /^[^#]/)
		    {
			$startcom = 1;
			$chaine = sprintf("! %s\n", $1);
			if ($line =~ /(.*)\*\//) {
			    $chaine =~ s/\*\///;
			}
			printTab( $chaine, $tabcount, 1);
		    }
		}
		elsif($line =~ /enum/)
		{
		    $startenum = 1;
		    $countenum = 0;		    
		    #$countenum = $countenum + 1 if (/PARAM/);    
		}
	    }
	    else 
	    {
		if ($line =~ /}/)
		{
		    $startenum = 0;
		    $countenum = 0;
		}
		elsif($line =~ /[ ]*([^ ]*[IR]PARAM[^ ]*)[ ]*=[ ]*([0-9]*),?/)
		{
		    $countenum = $2+1;
		    my $key       = $1;
		    #$countenum = $countenum + 1 if ($key =~ /PARAM/);
		    $chaine = sprintf($format, $key, $countenum);
		    printTab($chaine,$tabcount, 0);
		    $countenum++;
		}
		elsif($line =~ /[ ]*([^ ]*)[ ]*=[ ]*([0-9]*),?/)
		{
		    $countenum = $2;
		    my $key       = $1;
		    #$countenum = $countenum + 1 if ($key =~ /PARAM/);
		    $chaine = sprintf($format, $key, $countenum);
		    printTab($chaine,$tabcount, 0);
		    $countenum++;
		}
		elsif($line =~ /[ ]*([^ |^,]*)[ ]*,?/)
		{
		    my $key = $1;
		    chomp $key;
		    $chaine = sprintf($format, $key, $countenum);
		    printTab($chaine,$tabcount, 0);
		    $countenum++;
		}
		
		    
	    }
	}
	else
	{ 
	    if ($line =~ /^[ ]*> (.*)/)
	    {
		if ($interfaceprinted == 0)
		{
		    $chaine = "INTERFACE\n";
		    printTab($chaine, $tabcount);
		    $tabcount = 1;
		    $interfaceprinted = 1;
		}
		$chaine = sprintf("%s\n", $1);
		printTab($chaine, $tabcount, 0);
	    }
	    elsif ($line =~ /(.*)\*\//)
	    {
		$startcom = 0;
		$chaine = sprintf("! %s\n", $1);
		printTab($chaine, $tabcount, 1);
	    }
	    elsif($line =~ /(.*)/)
	    {

		
		$chaine = sprintf("! %s\n", $1);
		if ($line =~ /Murge's constants/)
		{
		    my $chaine2 = "END INTERFACE\n\n";
		    $chaine2 .= "!\n";
		    $chaine2 .= "!   Enum: KINDS\n";
		    $chaine2 .= "!\n";
		    $chaine2 .= "!   Solvers type kinds\n";
		    $chaine2 .= "!\n";
		    $chaine2 .= "!   Contains:\n";
		    $chaine2 .= "!     MURGE_COEF_KIND - Kind to use for COEF\n";
		    $chaine2 .= "!     MURGE_REAL_KIND - Kind to use for REAL\n";
		    $chaine2 .= "!     MURGE_INTS_KIND - Kind to use for INTS\n";
		    $chaine2 .= "!     MURGE_INTL_KIND - Kind to use for INTL\n";
		    $chaine2 .= "!\n";
		    if ($real != 0)
		    {
			$chaine2 .= "INTEGER, PARAMETER :: MURGE_COEF_KIND                = $real\n";
		    }	
		    if ($real != 0)
		    {
			$chaine2 .= "INTEGER, PARAMETER :: MURGE_REAL_KIND                = $real\n";
		    }	
		    if ($ints != 0)
		    {
			$chaine2 .= "INTEGER, PARAMETER :: MURGE_INTS_KIND                = $ints\n";
		    }	
		    if ($intl != 0)
		    {
			$chaine2 .= "INTEGER, PARAMETER :: MURGE_INTL_KIND                = $intl\n";
		    }
		    $tabcount --;
		    printTab($chaine2, $tabcount, 0);
		}

		printTab($chaine, $tabcount, 1);
	    }
	}
      
	
    }
    close APIc;
}


getopts("hf:c:r:s:l:",\%opts);

if ( defined $opts{c} ){
    $complex = $opts{c};
}
if ( defined $opts{r} ){
    $real = $opts{r};
}
if ( defined $opts{s} ){
    $ints = $opts{s};
}
if ( defined $opts{l} ){
    $intl = $opts{l};
}

if ( defined $opts{f} ){
    $fichier = $opts{f};
}
else {
  Usage();
  exit;
}

if ( defined $opts{h} ){
    Usage();
    exit;
}

Convert();
