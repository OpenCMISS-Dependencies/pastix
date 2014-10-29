#!/usr/bin/perl

#use strict;
use Getopt::Std;
use Modules::APIparser;

sub byindex {
    my $a1;
    my $b1;
    $a1 = $$a[1];
    $b1 = $$b[1];
    $a1 <=> $b1;
}

sub Convert {
    my $keyword;    
    my $index;
    my $definition;
    my $default;
    my $IO;
    my $chaine;
    my @tableau  = ();
    my $i;
    my $lastindex;
    my $format   = "%-10s#  %-60s index : %s %s\n";

    
    my @iparm    = ();
    my @dparm    = ();
    my @api      = ();
    my $filename = $ARGV[0];

    my ($refiparm, $refdparm, $refapi) = API_parser($filename);
    
    @iparm = @$refiparm;
    @dparm = @$refdparm;
    @api   = @$refapi; 

    #produce iparm.txt
    my $outf = "iparm.txt";
    printIparm($outf, @iparm);

    #produce dparm.txt
    my $outf = "dparm.txt";
    printDparm($outf, @dparm);
}


sub Usage {

    print "APIc2ltex.pl api.h";
    print "  convertit le fichier api.h en fichiers pour la refcard. ";
    print "   -h   Affiche cette aide\n";
}


getopts("h",\%opts);

if ( defined $opts{h} ){
    Usage();
    exit;
}

Convert();
