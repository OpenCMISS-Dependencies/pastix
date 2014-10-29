#!/usr/bin/perl

# use strict;
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
    my $i;
    my $title;

    my @iparm    = ();
    my @dparm    = ();
    my @api      = ();
    my $filename = $ARGV[0];

    my ($refiparm, $refdparm, $refapi) = API_parser($filename);
    

    @iparm = @$refiparm;
    @dparm = @$refdparm;
    @api   = @$refapi;

    @iparm = latexifyparam(@iparm);
    @dparm = latexifyparam(@dparm);


    #produce iparm.txt
    my $outf = "iparm.tex";
    open OUTF, ">$outf" or die "ERROR: could not open file $outf - $!\n";
    printf(OUTF "  \\begin{tabParameters}{Integer parameters and outputs.}{Integer parameters and outputs.}{iparm_array} \n");

    for $i (0..$#iparm)
    {
	$keyword    = $iparm[$i][0];
	$index      = $iparm[$i][1];
	$definition = $iparm[$i][2];
	$default    = $iparm[$i][3];
	$IO         = $iparm[$i][5];
	if (!($definition =~ /.*IGNORE.*/))
	{
	    $chaine = sprintf("    %-40s &  %-2d   & %-60s & %-16s & %-6s \\\\\n", 
			      $keyword, $index, $definition, $default, $IO);
	    printf(OUTF $chaine);
	}
    }
    printf(OUTF "  \\end{tabParameters}\n");

    close OUTF;

    #produce dparm.txt
    my $outf = "dparm.tex";
    open OUTF, ">$outf" or die "ERROR: could not open file $outf - $!\n";
    printf(OUTF "  \\begin{tabParameters}{Floating point parameters and outputs.}{Floating point parameters and outputs.}{dparm_array} \n");

	
    for $i (0..$#dparm)
    {
	$keyword    = $dparm[$i][0];
	$index      = $dparm[$i][1];
	$definition = $dparm[$i][2];
	$default    = $dparm[$i][3];
	$IO         = $dparm[$i][4];
	if (!($definition =~ /.*IGNORE.*/))
	{
	    $chaine = sprintf("    %-40s &  %-2d   & %-60s & %-16s & %-6s \\\\\n", 
			  $keyword, $index, $definition, $default, $IO);
	    printf(OUTF $chaine);
	}
    }
    printf(OUTF "  \\end{tabParameters}\n");

    close OUTF;


    my $outf = "api.tex";
    open OUTF, ">$outf" or die "ERROR: could not open file $outf - $!\n";
    my $old_api = "abcdef";
    #produce api.tex
    for $i (0..$#api)
    {
	if ($api[$i][1] != -1)
	{
	    if (!($old_api =~ /$api[$i][0]/) )
	    {
		if (!($old_api =~ /abcdef/) )
		{
		    printf(OUTF "    \\hline\n");
		    printf(OUTF "  \\end{tabIbis}\n\n\n");
		}
		printf(OUTF "\\label{$api[$i][3]-modes}\n");
		$title =  latexifyword($api[$i][2]);
		printf(OUTF "  \\begin{tabIbis}{$title}\n");
		printf(OUTF "    \\hline\n");
	    }
	    $keyword    = latexifyword($api[$i][4]);
	    $index      = $api[$i][5];
	    $definition = latexifyword($api[$i][6]);

	    $chaine = sprintf("    %-40s &  %-2d   & %-60s\\\\\n", 
			      $keyword, $index, $definition);
	    printf(OUTF $chaine);
	    
	    $old_api = $api[$i][0];
	}
    }
    printf(OUTF "    \\hline\n");
    printf(OUTF "  \\end{tabIbis}\n\n\n");
    close OUTF;
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
