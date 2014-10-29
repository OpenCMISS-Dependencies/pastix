#!/usr/bin/perl

use Getopt::Std;

sub Convert {

    $nbligne = $ARGV[0];
    open (M, $ARGV[1]); @m = <M>;  close M;
    
    open (I, $ARGV[2]);   
    open (J, $ARGV[3]);   
    open (V, $ARGV[4]);   
    open (S, ">".$ARGV[5]);   

    if ( $o_format =~ /ijv/ )
    {
	print S $m[0];
    }
    if ( $o_format == "mm" )
    {
	print S $m[0];
	print S $m[1];
    }

    for ($k = 0; $k<$nbligne; $k++){
	
	$i = <I>;
	$j = <J>;
	$v = <V>;
	
	chomp $i;
	chomp $j;
	chomp $v;

	print S $i." ".$j." ".$v."\n";
    }
    
    close I;
    close J;
    close V;
    close S;
}

sub Usage {

    print "convert.pl [ -o format ] NNZ header i j v output \n";
    print "  convertit 3 fichiers contenant chacun respectivement \n";
    print "  i, j et v en un seul au format ijv\n\n";
    print "   -h   Affiche cette aide\n";
    print "   -o   gives output format (default ijv)\n";
    print "        formats : MatrixMarket \"mm\"\n";
    print "                  ijValues     \"ijv\"\n";
}


getopts("ho:",\%opts);

if ( defined $opts{h} ){
    Usage();
    exit;
}

$o_format = "ijv";
if (defined $opts{o}){
    $o_format = $opts{o};
}

Convert();
