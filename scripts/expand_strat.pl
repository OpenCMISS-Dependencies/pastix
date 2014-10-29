#!/usr/bin/perl

use POSIX;
use strict;

sub expand {
    my ($filename) = @_;
    my $deep = 0;
    my $line = "";
        
    open(FILE,  "$filename") or die "Couldn’t open file \"$filename\"\n";
    open(TFILE, ">temp") or die "Couldn’t create temp file\n";
    # Parcours du fichier de log
    while ($line = <FILE>)
    {
	chomp $line;
	$line =~ s/{/{\n/g;
	$line =~ s/}/\n}\n/g;
	$line =~ s/,/,\n/g;
	$line =~ s/}\n,/},/g;
	$line =~ s/\|/\|\n/g;
	print TFILE "$line\n";	
    }
    close FILE;
    close TFILE;
    close FILE;
    open(TFILE, "<temp");
    while($line = <TFILE>) {
	$deep = $deep - 1 if ($line =~ /.*}.*/);
	for(my $i=0;$i<$deep;$i++) { $line="  $line";}
	$deep = $deep + 1 if ($line =~ /.*{.*/);
	print $line;
    }
    close TFILE;
}

my $ifile = $ARGV[0];
&expand($ifile);

