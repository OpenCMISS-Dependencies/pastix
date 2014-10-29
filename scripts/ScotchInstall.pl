#!/usr/bin/perl

$machinename=`uname -n`;


$toto  = "prefix  = _PREFIX_\n";
$toto .= "EXE     =               \n";
$toto .= "LIB     = .a		  \n";
$toto .= "OBJ     = .o		  \n";
$toto .= "\n";
$toto .= "AR      = ar		  \n";

if ( $machinename =~ "decrypthon" )
{
    $toto .= "MAKE    = gmake	  \n";
    $toto .= "CCS     = xlc_r		  \n";
    $toto .= "CCP     = mpcc_r            \n";
    $toto .= "CCD     = xlc_r -I/usr/lpp/ppe.poe/include\n";
    $toto .= "ARFLAGS = -X32_64 -ruv	  \n";
    $toto .= "CFLAGS  = -ma -q64 -qarch=auto -O3 -qstrict -qtune=auto -qlanglvl=extc99 -s -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_PTHREAD -DSCOTCH_RENAME _TYPE_ -Drestrict= -DSCOTCH_COLLECTIVE\n";
    $toto .= "LDFLAGS = -bmaxdata:0x80000000 -lm\n";
    $toto .= "LEX     = lex	  \n";
    $toto .= "YACC    = yacc          \n";
}
elsif ( $machinename =~  /vargas/ )
{
    $toto .= "MAKE    = gmake	  \n";
    $toto .= "CCS     = xlc_r		  \n";
    $toto .= "CCP     = mpcc_r            \n";
    $toto .= "CCD     = xlc_r -I/usr/lpp/ppe.poe/include\n";
    $toto .= "ARFLAGS = -X32_64 -ruv	  \n";
    $toto .= "CFLAGS  = -ma -q64 -qarch=auto -O3 -qstrict -qtune=pwr6 -qlanglvl=extc99 -s -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_PTHREAD -DSCOTCH_RENAME _TYPE_ -Drestrict= -DSCOTCH_COLLECTIVE\n";
    $toto .= "LDFLAGS = -lpthread -lm\n";
    $toto .= "LEX     = lex	  \n";
    $toto .= "YACC    = yacc          \n";
}
else
{
    $toto .= "MAKE    = make	  \n";
    $toto .= "CCS     = gcc		  \n";
    $toto .= "CCP     = mpicc         \n";
    $toto .= "CCD     = mpicc         \n";
    $toto .= "ARFLAGS = -ruv	  \n";
    $toto .= "CFLAGS  = -O3 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_PTHREAD -DSCOTCH_RENAME _TYPE_ -Drestrict= -DSCOTCH_COLLECTIVE \n";
    $toto .= "LDFLAGS = -lz -lm -lrt  \n";
    $toto .= "LEX		= flex -Pscotchyy -olex.yy.c\n";
    $toto .= "YACC    = bison -pscotchyy -y -b y\n";
}


$toto .= "CAT     = cat		  \n";
$toto .= "CP      = cp		  \n";
$toto .= "LN      = ln	          \n";
$toto .= "MKDIR   = mkdir	  \n";
$toto .= "MV      = mv		  \n";
$toto .= "RANLIB  = ranlib	  \n";

$tutu  = "";
$tutu .= "LIB=scotch\n";
$tutu .= "\n";
$tutu .= "export SCOTCH_HOME=_PREFIX_\n";
$tutu .= "\n";
$tutu .= "for i in PATH LD_LIBRARY_PATH LD_RUN_PATH INCLUDE_PATH\n";
$tutu .= "do\n";
$tutu .= "\n";
$tutu .= '  cmd1="echo \$$i | sed -r \'s/(\(.*:\)|)[^:]*${LIB}[^:]*(|\(:.*\))/\1\2/\'"'."\n";
$tutu .= '  cmd2="echo \$temp | sed \'s/::/:/\' | sed \'s/^://\' | sed \'s/:$//\' "'."\n";
$tutu .= "\n";
$tutu .= '  temp=`eval $cmd1`;'."\n";
$tutu .= '  temp=`eval $cmd2`;'."\n";
$tutu .= '  eval "$i=$temp";'."\n";
$tutu .= "done\n";
$tutu .= "\n";
$tutu .= 'export PATH=$PATH:$SCOTCH_HOME/bin'."\n";
$tutu .= 'export LD_RUN_PATH=$LD_RUN_PATH:$SCOTCH_HOME/lib'."\n";
$tutu .= 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SCOTCH_HOME/lib'."\n";
$tutu .= 'export INCLUDE_PATH=$INCLUDE_PATH:$SCOTCH_HOME/include'."\n";

my %liste = (
    "int"   => "",
    "int32" => "-DINTSIZE32",
    "int64" => "-DINTSIZE64",
    "long"  => "-DLONG"
    );

if ($#ARGV < 0)
{
    print "Preciser le repertoire d'installation\n";
    exit;
}
$prefix = $ARGV[0];
foreach $lib ("scotch", "ptscotch")
{
    foreach $version (keys %liste)
    {
	$Makefile = $toto;
	$env      = $tutu;

	$rep = $prefix."/".$lib."/".$version;
    
	$Makefile =~ s/_PREFIX_/$rep/;
	$Makefile =~ s/_TYPE_/$liste{$version}/;
	$env      =~ s/_PREFIX_/$rep/;
   
	print $Makefile."\n";
	open(M, ">Makefile.inc");
	print M $Makefile;
    close M;

	system("make realclean") ;
	system("make $lib") ;
	if ( -w $prefix )
	{
	    system("mkdir -p $rep");
	    system("make install") ;
	}
	else
	{
	    system("mkdir -p $rep");
	    system("make install") ;
	}
 
	open(M, ">$rep/env.sh");
	print M $env;
	close M;
    }
}
