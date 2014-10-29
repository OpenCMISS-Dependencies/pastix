##################################################
#
#   Package: Modules::APIparser
#
#       Module to parse api.h file
#################################################

use strict;
use warnings;
use Getopt::Std;

#
# Function: byindex
#
# Compares the second element of to tabulars
#
sub byindex {
    my $a1;
    my $b1;
    $a1 = $$a[1];
    $b1 = $$b[1];
    $a1 <=> $b1;
}

my $formatparam   = "%-10s#  %-60s index : %s %s\n";
#
# Function: byindex
#
# Compares the second element of to tabulars, 
# if equals compare the fourth
#
sub byposindex {
    my $a1;
    my $b1;
    $a1 = $$a[1];
    $b1 = $$b[1];
    if ($a1 == $b1)
    {
	$a1 = $$a[5];
	$b1 = $$b[5];
    }
    $a1 <=> $b1;
}

#
# Function: API_Parser
#
# Construct *iparms*, *dparms* and *api* tabulars
# from api.h file.
#
# *iparms* contains 6 columns and is sorted 
# on value :
#
# > KEYWORD VALUE DEFINITION DEFAULT DEFAULTVALUE IN/OUT
#
# *dparms* contains 5 columns and is sorted
# on value :
#
# > KEYWORD VALUE DEFINITION DEFAULT IN/OUT
#
# *api* contains 7 columns and is sorted iwith 
# to keys : position and value
#
# > ENUM_ID POSITION TITLE LABEL KEYWORD VALUE DEFINITION 
#
#
sub API_parser {
    my ($filename) = @_;

    my @iparms               = ();
    my @dparms               = ();
    my @api                  = ();

    my %iparms_values        = ();
    my %iparms_definition    = ();
    my %iparms_default       = ();
    my %iparms_IO            = ();

    my %dparms_values        = ();
    my %dparms_definition    = ();
    my %dparms_default       = ();
    my %dparms_IO            = ();

    my %api_definition   = ();
    my %api_index        = ();
    my %api_title        = ();
    my %api_label        = ();
    my %api_pos          = ();

    my $api_id;
    my $position;
    my $keyword;    
    my $index;
    my $definition;
    my $default;
    my $defaultval;
    my $IO;
    my $title;
    my $label;
    my $chaine;
    my @tableau = ();
    my $currentapi;    

    my $iparmkey      = "^[ \t]*(IPARM_[0-9A-Z_]*)[ \t]*=[ \t]*([0-9]*),?";
    my $iparmdefkey   = "^[ \t]*(IPARM_[0-9A-Z_]*)[ \t]*-[ \t]*(.*)Default:(.*)(IN|OUT)";
    my $dparmkey      = "^[ \t]*(DPARM_[0-9A-Z_]*)[ \t]*=[ \t]*([0-9]*),?";
    my $dparmdefkey   = "^[ \t]*(DPARM_[0-9A-Z_]*)[ \t]*-[ \t]*(.*)Default:(.*)(IN|OUT)";
    my $enumkey       = "^[ \t]*Enum:[ \t]*(API_[0-9A-Z_]*)[ \t]*.*";
    my $apikey        = "^[ \t]*(API_[0-9A-Z_]*)[ \t]*=[ \t]*([0-9]*),?";
    my $apidefkey     = "^[ \t]*(API_[0-9A-Z_]*)[ \t]*-[ \t]*(.*)";
    my $apimodekey    = "^[ \t]*(.*modes.*)"; 
    my $apiposkey     = "\/\* \_POS\_ (-?[0-9]*) .\/";
    
    open (APIc, $filename); 

    #Parse api.h
    foreach (<APIc> )
    {
	if (/$iparmkey/)
	{
	    $iparms_values{$1} = $2;
	}
	elsif (/$iparmdefkey/)
	{
	    $iparms_definition{$1} = $2;
	    $iparms_default{$1}    = $3;
	    $iparms_IO{$1}         = $4;
	}
	elsif (/$dparmkey/)
	{
	    $dparms_values{$1} = $2;
	}
	elsif (/$dparmdefkey/)
	{
	    $dparms_definition{$1} = $2;
	    $dparms_default{$1}    = $3;
	    $dparms_IO{$1}         = $4;
	}
	elsif (/$enumkey/)
	{
	    $currentapi = $1;
	    $api_index{$currentapi}      = ();
	    $api_definition{$currentapi} = ();
	}
	elsif (/$apikey/)
	{
	    $api_index{$currentapi}{$1} = $2;
	}
	elsif (/$apidefkey/)
	{
	    my $tmp  = $2;
	    my $tmp1 = $1;
	    $api_definition{$currentapi}{$tmp1} = $tmp;
	}
	elsif (/$apimodekey/)
	{
	    my $tmp = $1;
	    $api_title{$currentapi} = $tmp ; 
	    $tmp =~ s/(.*) modes.*/$1/;
	    $tmp =~ s/\ /-/g;
	    $api_label{$currentapi} = $tmp
	}
	elsif (/$apiposkey/)
	{
	    $api_pos{$currentapi} = $1;
	}
    }
    close APIc;

    # Build iparms data set
    foreach my $k (keys(%iparms_values))
    {
	$keyword    = $k;
	$index      = $iparms_values{$k};
	$definition = $iparms_definition{$k};
	$definition =~ s/^[ \t]+//;
	$definition =~ s/[ \t]+$//;
	$default    = $iparms_default{$k};	
	$default    =~ s/^[ \t]+//;
	$default    =~ s/[ \t]+$//;
	$defaultval = $default;
	foreach my $k (keys(%api_index))
	{
	    foreach my $k2 (keys(%{$api_index{$k}}))
	    {
		if ($k2 =~ /$defaultval/)
		{
		    $defaultval = $api_index{$k}{$k2};
		}
	    }
	}
	$IO         = $iparms_IO{$k};
	$iparms[$#iparms+1] = [$keyword, $index, $definition, $default, $defaultval, $IO];
    }

    @iparms = sort byindex @iparms;


    # Build dparms data set
    foreach my $k (keys(%dparms_values))
    {
	$keyword    = $k;
	$index      = $dparms_values{$k};
	$definition = $dparms_definition{$k};
	$default    = $dparms_default{$k};
	$default    =~ s/^\s+//;
	$default    =~ s/\s+$//;
	$defaultval = $default;
	$IO         = $dparms_IO{$k};
	$dparms[$#dparms+1] = [$keyword, $index, $definition, $default, $IO];
    }

    @dparms = sort byindex @dparms;


    foreach $api_id (keys(%api_index))
    {
    	$position   = $api_pos{$api_id};
	$position   =~ s/^\s+//;
	$position   =~ s/\s+$//;
	$label      = $api_label{$api_id};
	$title      = $api_title{$api_id};

	# Build api data set
	foreach $keyword (keys(%{$api_index{$api_id}}))
	{
	    $index      = $api_index{$api_id}{$keyword};
	    $definition = "";
	    $definition = $api_definition{$api_id}{$keyword} if (defined($api_definition{$api_id}{$keyword}));
	    $definition =~ s/^\s+//;
	    $definition =~ s/\s+$//;
	    
	    $api[$#api+1] = [$api_id, $position, $title, $label, $keyword, $index, $definition];
	}
    }
    
    @api = sort byposindex @api;

    return (\@iparms, \@dparms, \@api);
}

sub latexifyparam
{
    my (@params) = @_;
    my $tmp;

    for my $i (0..$#params)
    {	
	$params[$i][0] = latexifyword($params[$i][0]);
	$params[$i][2] = latexifyword($params[$i][2]);
	$params[$i][3] = latexifyword($params[$i][3]);
    }
    return @params;
	
}

sub latexifyword
{
    my (@arg) = @_;
    my $word;
    $word = $arg[0];
    $word =~ s/\_/\\\_/g if (!($word =~ /.*\$.*/ ));
    $word = "\$$word\$" if ($word =~ /.*e\^.*/);
    $word =~ s/([ID]PARM[^ ]*)/\\texttt{$1}/g;
    $word =~ s/(API[^ ]*)/\\texttt{$1}/g;
    $word =~ s/(PaStiX[^ ]*)/\\pastix{}/g;
    $word =~ s/(Scotch[^ ]*)/\\scotch{}/g;
    
    return $word;
}

sub printIparm
{
    my ($outf, @iparm) = @_;

    my $keyword;    
    my $index;
    my $definition;
    my $default;
    my $IO;
    my $chaine;
    my $lastindex;
    my $i;
    my $iparm_size = 64;

    open OUTF, ">$outf" or die "ERROR: could not open file $outf - $!\n";
    $lastindex = 0;
    for $i (0..$#iparm)
    {
	$keyword    = $iparm[$i][0];
	$index      = $iparm[$i][1];
	$definition = $iparm[$i][2];
	$default    = $iparm[$i][4];
	$default    = 0 if ($default eq "-" );
	$default    = 0 if ($default eq "");
	$IO         = $iparm[$i][5];

	if ( $keyword =~ /IPARM_SIZE/ )
	{
	    $iparm_size = $index;
	    last;
	}
	while ($lastindex < $index)
	{
	    
	    $chaine = sprintf($formatparam, 0, "unused", $lastindex, " ");
	    printf (OUTF $chaine);
	    $lastindex++;
	}
	$chaine = sprintf($formatparam, $default, $definition, $keyword, "($IO)");
	printf (OUTF $chaine);
	$lastindex++; 
    }
    while ($lastindex < $iparm_size)
    {
	$chaine = sprintf($formatparam, 
			  0,"unused", $lastindex, "");
	printf (OUTF $chaine);
	$lastindex++;
    }
    close OUTF;
}

sub printDparm
{
    my ($outf, @dparm) = @_;

    my $keyword;    
    my $index;
    my $definition;
    my $default;
    my $IO;
    my $chaine;
    my $lastindex;
    my $i;
    my $dparm_size = 64;

    open OUTF, ">$outf" or die "ERROR: could not open file $outf - $!\n";

    $lastindex = 0;
    for $i (0..$#dparm)
    {
	$keyword    = $dparm[$i][0];
	$index      = $dparm[$i][1];
	$definition = $dparm[$i][2];
	$default    = $dparm[$i][3];
	$default    = 0 if ($default eq "-" );
	$default    = 0 if ($default eq "");
	$default    =~ s/[\^\{\}]//g;
	$IO         = $dparm[$i][4];

	if ( $keyword =~ /DPARM_SIZE/ )
	{
	    $dparm_size = $index;
	    last;
	}
	while ($lastindex < $index)
	{
	    
	    $chaine = sprintf($formatparam, 0, "unused", $lastindex, " ");
	    printf (OUTF $chaine);
	    $lastindex++;
	}
	$chaine = sprintf($formatparam, 
			  $default,$definition,$index, "($IO)");
	printf(OUTF $chaine);  
	$lastindex++;
    }
    while ($lastindex < $dparm_size)
    {
	
	$chaine = sprintf($formatparam, 0, "unused", $lastindex, " ");
	printf (OUTF $chaine);
	$lastindex++;
    }
    close OUTF;
}

sub setiparm
{
    my ($keyword, $value, @params) = @_;
    my $i;
    for $i (0..$#params)
    {
	$params[$i][4] = $value if ($params[$i][0] eq $keyword);
    }
    return @params;
}
sub setdparm
{
    my ($keyword, $value, @params) = @_;
    my $i;
    for $i (0..$#params)
    {
	if ($params[$i][0] eq $keyword)
	{
	    print "$params[$i][0]: $params[$i][3] -> $value\n";
	}
	$params[$i][3] = $value if ($params[$i][0] eq $keyword);
    }
    return @params;
}

sub getiparm
{
    my ($keyword, @params) = @_;
    my $i;
    for $i (0..$#params)
    {
      return $params[$i][4] if ($params[$i][0] eq $keyword);
    }
    return -1;
}1;
