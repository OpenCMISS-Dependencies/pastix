#!/usr/bin/env perl

###############################################################################
#
#   Script: geninterface.pl
#
###############################################################################
#
#   A Script that convert <murge.h> into a C to Fortran interface.
#
#   Usage:
#    > ./geninterface.pl -f murge.h
#    >    converts the file murge.h into a fortran header.
#    >     -h          Shows this help
#    >     -n          Use fortran numbering
#    >                 (decrement index in murge_setoption[int|real])
#
#   Authors:
#     Pascal Jacq    - .
#     Xavier Lacoste - lacoste@labri.fr
#
###############################################################################
use POSIX;
use strict;
use warnings;
use Data::Dumper;
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
#   hash: opts
#     Options given to the script
my $real     = 0;
my $ints     = 0;
my $intl     = 0;
my $complex  = 0;
my $fortran_numbering = 0;
my $fichier;
my %opts;

###############################################################################
# Group: Functions

#
# Function: Usage
#
# Prints usage.
sub Usage {

    print "./geninterface.pl -f murge.h\n";
    print "  converts the file murge.h into a fortran header.\n";
    print "   -h          Shows this help\n";
    print "   -n          Use fortran numbering\n";
    print "               (decrement index in murge_setoption[int|real]\n";
}

#
# Function: getTypePtrNames
#
# Retrieve from the arguments list, the list of arguments names, with added ierror,
# separated with comas.
#
# Also build the list of arguments, as pointers, the list of arguments names preceded by a star
# and the name of the string if there is one.
#
# Parameters:
#   chaine - String containing the argument list.
#
# Returns:
#   The four builded strings : arguments names, pointers, arguments star prefixed, string name.
#
sub getTypePtrNames # ($chaine)
  {
    my ($chaine) = @_;

    my $types = "";
    my $ptrs = "";
    my $names = "";
    my $type;
    my $arg;
    my @args;
    my $isptr;
    my $hasChaine = "";
    my $hasComm = "";

    $chaine =~ s/\r//g;
    @args = split(',' , $chaine);
    if ($#args+1 > 0)
    {
	foreach( @args )
	{
	    $isptr = 0;
	    $arg = $_;
	    if ($arg =~ /\*/)
	    {
		$types .= $arg .", ";
		$arg =~ s/.*\*(.*)/$1, /;
		$ptrs .= $arg;
		$names .= $arg;
	    }
	    else
	    {
		$arg =~ s/.* (.*)/$1/;
		$names .= "$arg, ";

		$type = $_;
		$type =~ s/(.*) .*/$1/;
		if ($type =~ /MPI_Comm/)
		{
		    $type = "MPI_Fint";

		    $hasComm = $arg;
		    $ptrs .= "tmpcomm, ";
		}
		else
		{
		    $ptrs .= "*$arg, ";
		}
		$types .= $type . " *" . $arg .", ";
	    }

	    if ($_ =~ /char[ ]*\*/)
	    {
		$hasChaine = $_;
		$hasChaine =~ s/char[ ]*\*(.*)/$1/;
		$hasChaine =~ s/^[ \t]+//;

		$types .= "INTS *stringlength, ";
		$names .= "stringlength, ";
		$ptrs =~ s/$hasChaine/tmp/;
	    }
	}
    }
    $names, $types, $ptrs, $hasChaine, $hasComm;
}

#
# Function: Convert
#
# Main function.
#
# Converts the header <murge.h> file into a C to Fortran interface.
#
sub Convert {

    my $startcom  = 0;
    my $startenum = 0;
    my $countenum;
    my $chaine;
    my $hasChaine;
    my $hasComm;
    my $tabcount = 0;
    my $fct_name;
    my $fct_name_lc;
    my $fct_name_uc;
    my $fct_args;
    my $fct_args_name;
    my $fct_args_type;
    my $fct_args_ptr;


    print "#ifdef _WIN32\n";
    print "#  define MURGE_DLL_EXPORT\n";
    print "#endif\n";

    print "#include \"murge.h\"\n\n";

    print "/** Macro from Scotch **/\n\n";
    print "#define FORTRAN_NAME(nu,nl,pl,pc) \\\n";
    print "void nu pl;                       \\\n";
    print "void nl pl                        \\\n";
    print "{ nu pc; }                        \\\n";
    print "void nl##_ pl                     \\\n";
    print "{ nu pc; }                        \\\n";
    print "void nl##__ pl                    \\\n";
    print "{ nu pc; }                        \\\n";
    print "void nu pl\n\n";



    open (APIc, $fichier);

    foreach (<APIc> )
    {
	if ($startcom == 0)
	{
	    # If we are not in comments area
	    if ($startenum == 0)
	    {
		# If we are not in an Enum
		if (/DECLSPEC/)
		{
		    #ignore DECLSPEC
		}
		elsif (/\/\*(.*)/)
		{
		    # If we start a commentary
		    $startcom = 1;
		}
		elsif(/enum/)
		{
		    # If we start an enum
		    $startenum = 1;
		    $countenum  = 0;
		    #$countenum = $countenum + 1 if (/PARAM/);
		}
		elsif(/^#/)
		{
		    # Preprocessor command
		}
		elsif(/(.*)/)
		{
		    # Else, function declaration
		    my $temp = $1;
		    $temp =~ s/^[ \t]+//;
		    $temp =~ s/[ \t]+$//;
		    $chaine .= sprintf("%s", $temp, 1);
		    if ($chaine =~ /;/)
		    {
			# If we end the function
			$fct_name = $chaine;
			$fct_name =~ s/INTS (.*)\(.*/$1/;
			$fct_args = $chaine;
			$fct_args =~ s/INTS .*\((.*)\).*/$1/;

			chomp $fct_name;
			chomp $fct_args;

			($fct_args_name, $fct_args_type, $fct_args_ptr, $hasChaine, $hasComm) = getTypePtrNames($fct_args);

			chop $fct_args_ptr;
			chop $fct_args_ptr;
			if ($fct_name ne "")
			{
			    $fct_name_lc = lc $fct_name;
			    $fct_name_uc = uc $fct_name;
			    if ($fct_args_type =~ /MPI/)
			    {
				print "#ifdef MPI_VERSION\n";
			    }
			    print "FORTRAN_NAME($fct_name_uc,\n";
			    print "             $fct_name_lc,\n";
			    print "             ($fct_args_type INTS *ierror),\n";
			    print "             ($fct_args_name ierror))\n";
			    print "{\n";
			    if ($hasChaine ne "")
			    {
				print "  char * tmp = NULL;\n";
				print "  tmp = (char *) malloc ((*stringlength+1)*sizeof(char));\n";
				print "  strncpy(tmp, $hasChaine, *stringlength);\n";
				print "  tmp[*stringlength] = '\\0';\n";

			    }
			    if ($hasComm ne "")
			    {
				print "  MPI_Comm tmpcomm;\n";
				print "  tmpcomm = MPI_Comm_f2c(*$hasComm);\n"
			    }
			    if ( $fct_name_lc =~ "murge_setoptionint" ||
				 $fct_name_lc =~ "murge_setoptionreal")
			    {
				if ($fortran_numbering)
				{
				    print "  int my_number = *number - 1;\n";
				}
				else
				{
				    print "  int my_number = *number;\n";

				}
				$fct_args_ptr =~ s/\*number/my_number/g;
			    }
			    print "  *ierror = $fct_name($fct_args_ptr);\n";
			    print "}\n";
			    if ($fct_args_type =~ /MPI/)
			    {
				print "#endif /* MPI_VERSION */\n";
			    }
			    print "\n";
			}
			$chaine = "";
		    }
		    else
		    {
			chomp $chaine;
		    }
		}
	    }
	}
	else
	{
	  if (/(.*)\*\//)
	    {
		$startcom = 0;
	    }
	}

    }
    close APIc;
}

getopts("hnf:",\%opts);
if ( defined $opts{h} ) {
  Usage();
  exit;
}

if ( defined $opts{n} ){
    $fortran_numbering = 1;
}

if ( defined $opts{f} ){
    $fichier = $opts{f};
    Convert();
}
else {
    Usage();
    exit;
}

