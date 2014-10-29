#!/usr/bin/perl

###############################################################################
#
#   Script: PastixDatabase.pm
#
###############################################################################
#
#   Perl module implementing functions to fill in pastix database.
#
#   Authors:
#     Xavier Lacoste - lacoste@labri.fr
#
###############################################################################

package Modules::PastixDatabase;
use Term::ANSIColor;
use strict;
use warnings;
use vars qw(@ISA @EXPORT);
use Exporter;
use Modules::Common;
use DBI;
use Modules::LogParser;

###############################################################################
# Group: Variables

#
#
#   Array: ISA
#     Inheritance array
#
#   Array: EXPORT
#     List of symbols to export

our @ISA    = qw( Exporter );
our @EXPORT = qw( checkMatrixPresence  checkMachinePresence fillCompiler 
                  fillExecutable nextExecutionEntry nextProgrammeEntry 
                  addExecution lastProgrammeEntry);
###############################################################################
# Group: Functions

#
# Function: checkMatrixPresence
#
# Checks that the matrix already exists in Database, 
# otherwise adds the matrix.
#
# Parameters:
#   matrix - Name of the matrix
sub checkMatrixPresence # ($matrix)
{
    my ($matrix) = @_;
    my $dbh = DBI->connect($connectionInfo,$userid,$passwd);
    my $query = "SELECT * from MATRICE where NOM = '$matrix'";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    if ($sth->fetch)
    {
	$sth->finish();
    }
    else
    {
	$sth->finish();
	$query = "INSERT INTO `MATRICE` ( `NOM`) VALUES ( '$matrix');";
	$dbh->do($query);
	#print "$query\n";
    }
    $dbh->disconnect();
}

#
# Function: checkMachinePresence
#
# Checks that the machine already exists in Database, 
# otherwise adds the machine.
#
# Parameters:
#   machinename - Name of the machine
#
sub checkMachinePresence # ($machinename)
{
    my ($machinename) = @_;
    my $dbh = DBI->connect($connectionInfo,$userid,$passwd);
    my $query     = "SELECT * from MACHINE where NOM = '$machinename'";
    my $sth       = $dbh->prepare($query);
    $sth->execute();

    if ($sth->fetch)
    {
	$sth->finish();
#	    print "$machine déjà présente.\n";
    }
    else
    {
	$sth->finish();
	$query = "INSERT INTO `MACHINE` ( `NOM` ) VALUES ( '$machinename');";
	#print "$query\n";
	$dbh->do($query);
    }
    $dbh->disconnect();
}

#
# Function: fillCompiler
#
# Fills the COMPILER database if the entry does not exist. 
#
# Parameters:
#   nomprogramme  - Name of the programme
#   CompilationOK - Boolean indicating if the compilation worked 
#                   ( OK = 0, ERROR = 1)
#   nberr         - Number of error in log file
#   nbwarn        - Number of warning in log file
#   log           - Log file path
#   mydate        - Date of the compilation YYYYMMDD-HHMM
#   user          - User name
#
sub fillCompiler # ($nomprogramme, $CompilationOK, $nberr, $nbwarn, $log, $mydate, $user)
{
    my ($nomprogramme, $CompilationOK, $nberr, $nbwarn, $log, $mydate, $user) = @_;
    my $dbh = DBI->connect($connectionInfo,$userid,$passwd);
    my $query = "";
    my $sth;
    my $MACH;
    my $VERS;
    my $PROG;
    my $OK;
    my $annee = substr $mydate,0,4;
    my $mois = substr $mydate,4,2;
    my $jour = substr $mydate,6,2;
    my $heure= substr $mydate,9,2;
    my $min= substr $mydate,11,2;
    my $day = "$annee-$mois-$jour";
    chomp $day;
    my $hour = "$heure:$min";
    chomp $hour;
    my $logdest;
    my $programme = lastProgrammeEntry($nomprogramme);

    $logdest = "log-compile/$machine-$programme-$date" ;
    `cp $log $SITE_PATH/$logdest`;

    # Check if the entry does not exist
    $query = "SELECT * from COMPILER where MACHINE = '$machine'"
	." and PROGRAMME = '$programme' and `UTILISATEUR` = '$user' and ".
	"`DATE` = '$day' and `HEURE` = '$hour'";
    $sth = $dbh->prepare($query);
    $sth->execute();

    if ($sth->fetch)
    {
	#Already added
	$sth->finish();
    }
    else
    {
	$sth->finish();
	$query = "INSERT INTO `COMPILER` ( `MACHINE` , `PROGRAMME` , "
	    . "`UTILISATEUR` , `RESULTAT` , `ERRORS`, `WARNINGS`, "
	    . " `DATE` , `HEURE`)"
	    . " VALUES ( '$machine', '$programme' , '$user' , "
	    . "'$CompilationOK' , '$nberr', '$nbwarn', "
	    . " '$day' , '$hour');";
	#print "$query\n";
	$dbh->do($query);
	
	$dbh->disconnect();
	return $logdest;
    }
}

#
# Function: fillExecutable
#
# Fills the PROGRAMME database. 
#
# Parameters:
#   nomprogramme  - Name of the programme
#   branche       - SVN branche of the programme
#   patch         - Path to the file corresponding to the svn diff
#   log           - Log file path
#   logerr        - Error log file path
#   config        - config.in file path
#
sub fillExecutable # ($nomprogramme, $branche, $patch, $log, $logerr, $config)
{
    my ($nomprogramme, $branche, $patch, $log, $logerr, $config) = @_;
    my $dbh = DBI->connect($connectionInfo,$userid,$passwd);

    my %metriques = LogGetMetriques($log, $logerr);
    my $version  = $metriques{'vers'};
    my $patchmd5 = `md5sum $patch`; 
    $patchmd5 =~ s/([^ ]*) .*/$1/;
    chomp $patchmd5;
    my $configmd5 = `md5sum $config`; 
    $configmd5 =~ s/([^ ]*) .*/$1/;
    chomp $configmd5;
    my $nextid = nextProgrammeEntry();
    `cp $patch $SITE_PATH/patch/$nextid`;
    #remplit la table programme
    my $query = "SELECT * from PROGRAMME where NOM = '$nomprogramme'"
	." and BRANCHE = '$branche' and `patch.md5` = '$patchmd5' and ".
	"`config.md5` = '$configmd5' and `VERSION` = '$version'";

    my $sth = $dbh->prepare($query);
    $sth->execute();


    if ($sth->fetch)
    {
	$sth->finish();
	print "executable déjà présent.\n";
    }
    else
    {
	$sth->finish();
	my $query = "INSERT INTO `PROGRAMME` ( `NOM`, `BRANCHE`, `patch.md5`,"
	    ." `config.md5`";
	foreach my $opt (@optlist)
	{
	    $query .= ",\n `".$metriqueconf{$opt}{'db'}."`";
	}
	$query .= ") VALUES ( '$nomprogramme' , '$branche', '$patchmd5', "
	    ."'$configmd5'";
	foreach my $opt (@optlist)
	{
	    $query .= ",\n '".$metriques{$opt}."'";
	}	
	$query .= ");";
	#print "$query\n";	
	$dbh->do($query);
    
    }
    $dbh->disconnect();
}

#
# Function: nextExecutionEntry
#
# Computes ID of the next execution entry. 
#
sub nextExecutionEntry
{
    my $dbh = DBI->connect($connectionInfo,$userid,$passwd);
    my $query = "SELECT `ID` FROM `EXECUTION` ORDER BY `ID` DESC LIMIT 0 , 1";
    my $sth = $dbh->prepare($query);
    my $numexec = 0;
    $sth->execute();

    $sth->bind_columns(\$numexec);
    $sth->fetch;
    $numexec++;
    $sth->finish;
    $dbh->disconnect();

    return $numexec;
}

#
# Function: nextProgrammeEntry
#
# Computes ID of the next programme entry. 
#
sub nextProgrammeEntry
{
    my $dbh = DBI->connect($connectionInfo,$userid,$passwd);
    my $query = "SELECT `ID` FROM `PROGRAMME` ORDER BY `ID` DESC LIMIT 0 , 1";
    my $sth = $dbh->prepare($query);
    my $numprog = 0;
    $sth->execute();

    $sth->bind_columns(\$numprog);
    $sth->fetch;
    $numprog++;
    $sth->finish;
    $dbh->disconnect();

    return $numprog;
}

#
# Function: lastProgrammeEntry
#
# Computes higher ID of the PROGRAMME database. 
#
sub lastProgrammeEntry
{
    my ($nomprogramme) = @_;
    my $dbh = DBI->connect($connectionInfo,$userid,$passwd);
    my $query = "SELECT `ID` FROM `PROGRAMME` WHERE `NOM` = '$nomprogramme' ORDER BY `ID` DESC LIMIT 0 , 1";
    my $sth = $dbh->prepare($query);
    my $numprog = 0;
    $sth->execute();

    $sth->bind_columns(\$numprog);
    $sth->fetch;
    $sth->finish;
    $dbh->disconnect();

    return $numprog;
}

#
# Function: addExecution
#
# Adds an entry to the Execution database. 
#
# Parameters:
#   log           - Log file path
#   logerr        - Error log file path
#   iparmlog      - iparm.txt file path
#   dparmlog      - dparm.txt file path
#   nomprogramme  - Name of the programme
#   branche       - SVN branche of the programme
#   patchmd5      - md5sum corresponding to the svn diff
#   configmd5     - config.in file md5sum
#   user          - User name
#   version       - Revision
#
sub addExecution # ($log, $logerr, $iparmlog, $dparmlog, $nomprogramme, $branche, $patchmd5, $configmd5, $user , $version)
{
    my ($log, $logerr, $iparmlog, $dparmlog, $nomprogramme, $branche, $patchmd5, $configmd5, $user, $version) = @_;
    my $dbh = DBI->connect($connectionInfo,$userid,$passwd);
    my @st;
    my $i;
    my $iparm;
    my $iparmkey;
    my $dparm;
    my $dparmkey;

    my $numexec   = nextExecutionEntry();
    my %metriques = LogGetMetriques($log, $logerr);
    my %params = ();
    my $progid = -1;
    my $query = " SELECT `ID`FROM `PROGRAMME`"
	." WHERE `NOM`      LIKE '$nomprogramme' "
	." AND `VERSION`    LIKE '$version' "
	." AND `BRANCHE`    LIKE '$branche'"
	." AND `patch.md5`  LIKE '$patchmd5'"
	." AND `config.md5` LIKE '$configmd5'"
	." ORDER BY `PROGRAMME`.`ID` DESC;";


    my $sth = $dbh->prepare($query);
    $sth->execute();
    $sth->bind_columns(\$progid);
    $sth->fetch;
    $sth->finish;

    if ($progid < 1)
    {
	print $query."\n";
	exit;
    }
    # Ouverture du fichier de log
    open (M, "$iparmlog"); 
    @st = <M>;  
    close M;
    
    # Parcours du fichier de log
    $i = 0;
    
    foreach (@st)
    {
	if ($i < $MAX_IPARM)
	{
	    $iparmkey = "";
	    $iparm = $1 if (/[ ]*([0-9]*)[ ].*index : IPARM_.*/);
	    $iparmkey = $1 if (/.*index : (IPARM_[^ ]*)/);
	    if ($iparmkey ne "")
	    {
		$params{$iparmkey} = $iparm;
	    }
	    $i++;
	}
    }

    # Ouverture du fichier de dparm
    open (M, "$dparmlog"); 
    @st = <M>;  
    close M;
    $i = 0;    
    # Parcours du fichier de dparm
    foreach (@st)
    {
	if ($i < $MAX_DPARM)
	{
	    $dparmkey = "";
	    $dparm = $1 if (/[ ]*([0-9]*)[ ].*index : DPARM_.*/);
	    $dparmkey = $1 if (/.*index : (DPARM_[^ ]*)/);
	    if ($dparmkey ne "")
	    {
		params{$dparmkey} = $dparm;
	    }
	    $i++;
	}
    }

    `cp $log $SITE_PATH/logs/$numexec` if ( -f $log );
    `cat $logerr >> $SITE_PATH/logs/$numexec` if ( -f $logerr );
    $query = "INSERT INTO `EXECUTION` ( `ID`, `PROGRAMME`, `UTILISATEUR`";
    foreach my $k (keys(%params))
    {
	$query .= ",\n `$k`";
    }
    foreach my $metr (@execlist)
    {
	$query .= ",\n `".$metriqueconf{$metr}{'db'}."` "
    }
    $query .= ") VALUES ( 'NULL', '$progid' , '$user'";
    foreach my $v (values(%params))
    {
	$query .= ",\n '$v'";
    }
    foreach my $metr (@execlist)
    {
	$query .= ",\n '".$metriques{$metr}."'";
    }
    $query .= ");";
    #print "$query\n";
    $dbh->do($query);
    $dbh->disconnect();
}
1;
