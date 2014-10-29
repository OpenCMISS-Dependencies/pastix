#!/usr/bin/perl

use POSIX;
use strict;

sub concat {
    my ($filename) = @_;
    my $deep = 0;
    my $line = "";
    my $out  = "";
        
    open(FILE,  "$filename") or die "Couldn’t open file \"$filename\"\n";

    # Parcours du fichier de log
    while ($line = <FILE>)
    {
	chomp $line;
	$line =~ s/ //g;
	$line =~ s/\n//g;
	$out  = "$out$line"
    }
    close FILE;
    print "$out\n"
}

my $ifile = $ARGV[0];
&concat($ifile);

