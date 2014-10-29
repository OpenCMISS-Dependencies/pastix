#!/usr/bin/perl

package Modules::LogParser;
use Term::ANSIColor;
use strict;
use warnings;
use vars qw(@ISA @EXPORT);
use Exporter;
use Modules::Common;

our @ISA    = qw( Exporter );
our @EXPORT = qw( LogGetMetriques 
                  LogPrintMetriques 
                  LogPrintMetriquesHeader 
                  LogPrintMetriquesHtml 
                  LogPrintMetriquesHeaderHtml 
                  LogPrintMetriquesHelp );

my $prefixdist = "d?";
####################################################
#             Compilation                          #
####################################################

## Parametre les options de compilation en fonction du cas
sub LogGetMetriques {
    
    my ($log, $logerr)     = @_ ;
    my %metriques = ();

    #Recuperation de la date
    $metriques{'date'} = $log;
    $metriques{'date'} =~ s/.*$metriqueconf{'date' }{'key'}.*/$1/;
    $metriques{'day'} = $log;
    $metriques{'day'}  =~ s/.*$metriqueconf{'day' }{'key'}.*/$1-$2-$3/;
    $metriques{'hour'} = $log;
    $metriques{'hour'} =~ s/.*$metriqueconf{'hour' }{'key'}.*/$1:$2:00/;
    $metriques{'exemple'} = $log;
    $metriques{'exemple'} =~ s/$metriqueconf{'exemple' }{'key'}/$1/;;

#	$data{$executable}{$mat}{$cas}{$filedate} = ();
	
	
    $metriques{'fact'} = -1.0;
    $metriques{'solv'} = -1.0;
    $metriques{'refin'}= -1.0;
    $metriques{'mem'}  = -1;
    $metriques{'mend'} = -1;
    $metriques{'ovh'}  = -1;
    $metriques{'ovhm'} = -1;
    $metriques{'ovhl'} = 950;
    $metriques{'fr2'}  = -1;
    $metriques{'iter'} = -1;
    $metriques{'norm'} = -1.0;
    $metriques{'warn'} = 0;
    $metriques{'err'}  = 0;
    $metriques{'vers'} = "";
    foreach my $opt (@optlist)
    {
	$metriques{$opt} = "";
    }
    # Si ce n'est pas un fichier on le passe
    return (%metriques) if ( ! -f $log );
    # Ouverture du fichier de log
    open (M, "$log");
	
    # Parcours du fichier de log
    foreach (<M>)
    {
	$metriques{'fact'}  = $1  if (/$metriqueconf{'fact' }{'key'}/ && $1 > $metriques{'fact'} );
	$metriques{'solv'}  = $1  if (/$metriqueconf{'solv' }{'key'}/ && $1 > $metriques{'solv'} );
	$metriques{'refin'} = $1  if (/$metriqueconf{'refin'}{'key'}/ && $1 > $metriques{'refin'});
	$metriques{'mem'}   = $1  if (/$metriqueconf{'mem'  }{'key'}/ );
	$metriques{'mend'}  = $1  if (/$metriqueconf{'mend' }{'key'}/ );
	$metriques{'ovh'}   = $1  if (/$metriqueconf{'ovh'  }{'key'}/ && $1 > $metriques{'ovh'}  );
	$metriques{'ovhm'}  = $1  if (/$metriqueconf{'ovhm' }{'key'}/ && $1 > $metriques{'ovhm'} );
	$metriques{'ovhl'}  = $1  if (/$metriqueconf{'ovhl' }{'key'}/ && $1 < $metriques{'ovhl'} && $1 > 0.5);
	$metriques{'fr2'}   = $1  if (/$metriqueconf{'fr2'  }{'key'}/ && $1 > $metriques{'fr2'}  && $1 > 0);
	$metriques{'iter'}  = $1  if (/$metriqueconf{'iter' }{'key'}/ );
	$metriques{'norm'}  = $1  if (/$metriqueconf{'norm' }{'key'}/ );
	$metriques{'norm'}  = $1  if (/$metriqueconf{'norm' }{'key2'}/ );
	$metriques{'warn'} += 1   if (/$metriqueconf{'warn' }{'key'}/ );
	$metriques{'err'}  += 1   if (/$metriqueconf{'err'  }{'key'}/ );
	$metriques{'sym'}   = 'Y' if (/$metriqueconf{'sym'  }{'key'}/ && ($1 =~ /LLt/));
	$metriques{'sym'}   = 'N' if (/$metriqueconf{'sym'  }{'key'}/ && ($1 =~ /LU/ ));
	foreach my $opt (@optlist)
	{
	    $metriques{$opt}   = $1 if (/$metriqueconf{$opt}{'key'}/ );
	}
	foreach my $opt (@paramlist2)
	{
	    $metriques{$opt}   = $1 if (/$metriqueconf{$opt}{'key'}/ );
	}
    }
    close M;
    if ( -f $logerr )
    {
	# Ouverture du fichier de logerr
	open (M, "$logerr");
	
	# Parcours du fichier de logerr
	foreach (<M>)
	{
	    $metriques{'fact'}  = $1  if (/$metriqueconf{'fact' }{'key'}/ && $1 > $metriques{'fact'} );
	    $metriques{'solv'}  = $1  if (/$metriqueconf{'solv' }{'key'}/ && $1 > $metriques{'solv'} );
	    $metriques{'refin'} = $1  if (/$metriqueconf{'refin'}{'key'}/ && $1 > $metriques{'refin'});
	    $metriques{'mem'}   = $1  if (/$metriqueconf{'mem'  }{'key'}/ );
	    $metriques{'mend'}  = $1  if (/$metriqueconf{'mend' }{'key'}/ );
	    $metriques{'ovh'}   = $1  if (/$metriqueconf{'ovh'  }{'key'}/ && $1 > $metriques{'ovh'}  );
	    $metriques{'ovhm'}  = $1  if (/$metriqueconf{'ovhm' }{'key'}/ && $1 > $metriques{'ovhm'} );
	    $metriques{'ovhl'}  = $1  if (/$metriqueconf{'ovhl' }{'key'}/ && $1 < $metriques{'ovhl'} && $1 > 0.5);
	    $metriques{'fr2'}   = $1  if (/$metriqueconf{'fr2'  }{'key'}/ && $1 > $metriques{'fr2'}  && $1 > 0);
	    $metriques{'iter'}  = $1  if (/$metriqueconf{'iter' }{'key'}/ );
	    $metriques{'norm'}  = $1  if (/$metriqueconf{'norm' }{'key'}/ );
	    $metriques{'norm'}  = $1  if (/$metriqueconf{'norm' }{'key2'}/ );
	    $metriques{'warn'} += 1   if (/$metriqueconf{'warn' }{'key'}/ );
	    $metriques{'err'}  += 1   if (/$metriqueconf{'err'  }{'key'}/ );
	    $metriques{'sym'}   = 'Y' if (/$metriqueconf{'sym'  }{'key'}/ && ($1 =~ /LLt/));
	    $metriques{'sym'}   = 'N' if (/$metriqueconf{'sym'  }{'key'}/ && ($1 =~ /LU/ ));
	    foreach my $opt (@optlist)
	    {
		$metriques{$opt}   = $1 if (/$metriqueconf{$opt}{'key'}/ );
	    }
	    foreach my $opt (@paramlist2)
	    {
		$metriques{$opt}   = $1 if (/$metriqueconf{$opt}{'key'}/ );
	    }
	}
	close M;
    }
    if (defined($metriques{'nnza'}) && ($metriques{'sym'} =~ /N/ ))
    {
	$metriques{'nnza'} = ($metriques{'nnza'} + $metriques{'size'}) / 2;
    }
    if (defined($metriques{'nnzl'}) && ($metriques{'sym'} =~ /N/ ))
    {
	$metriques{'nnzl'} = ($metriques{'nnzl'} + $metriques{'size'}) / 2;
    }

    $metriques{'vers'}  =~ s/\ //g;
    # On remplace Defined Not defined int64_t... par les valeurs correspondantes.
    foreach my $opt (@optlist)
    {
	$metriques{$opt} =~ s/[ ]*([A-Za-z0-9].*)/$1/;
	$metriques{$opt} =~ s/(.*[A-Za-z0-9])[ ]*/$1/;
	foreach my $k (keys(%keywords))
	{
	    my $reftmp = $keywords{$k};
	    my %tmp = %$reftmp;
	    foreach my $k2 (keys(%tmp))
	    {
		$metriques{$opt} = $tmp{$k2} if ($metriques{$opt}  eq $k2);
	    }
	}

	if (!($metriques{$opt} =~ /[0-9]/))
	{
	    $metriques{$opt} = -1;
	}
    }
    $metriques{'fact'}  = sprintf("%.2f", $metriques{'fact'});
    $metriques{'solv'}  = sprintf("%.2f", $metriques{'solv'});
    $metriques{'refin'} = sprintf("%.2f", $metriques{'refin'});
    $metriques{'ovh'}   = sprintf("%.2f", $metriques{'ovh'});
    $metriques{'ovhm'}  = sprintf("%.2f", $metriques{'ovhm'});
    $metriques{'ovhl'}  = sprintf("%.2f", $metriques{'ovhl'});
    $metriques{'fr2'}   = sprintf("%.2f", $metriques{'fr2'});


    #if ($RESULT_PATH eq $FILL_DB_PATH)
    #{
    #	$metriques{'mach'} = `uname -n`;
    #	chop $metriques{'mach'};
    #
    #	$metriques{'exec'}  = $log;
    #	#Suppression du répertoire de stockage
    #	$metriques{'exec'}  =~ s/$RESULT_PATH\///;
    #	# Tout ce qu'il y a avant le premier slash
    #	$metriques{'exec'}  =~ s/([_a-zA-Z0-9]*)\/.*/$1/;
    #		
    #    # Récupération du nom de la matrice 
    #	$metriques{'matr'}  = $log;
    #	$metriques{'matr'}  =~ s/$RESULT_PATH\///;
    #	$metriques{'matr'}  =~ s/$metriques{'exec'}\///;
    #	# Tout ce qu'il y a avant le premier slash
    #	$metriques{'matr'}  =~ s/([_a-zA-Z0-9]*)\/.*/$1/;
    #	
    #	# Récupération du cas XXP_XXT_... (suppression de tout ce qu'il y a avant le dernier /)
    #	$metriques{'cas'}   = $log;	    
    #	$metriques{'cas'}  =~ s/$RESULT_PATH\///;
    #	$metriques{'cas'}  =~ s/$metriques{'exec'}\///;
    #	$metriques{'cas'}  =~ s/$metriques{'matr'}\///;
    #	# Tout ce qu'il y a avant le premier slash
    #	$metriques{'cas'}  =~ s/([_a-zA-Z0-9]*)\/.*/$1/;
    #}
    #else
    #{
	
    $metriques{'user'} = $log;
    $metriques{'user'}  =~ s/$RESULT_PATH\///;
    # Tout ce qu'il y a avant le premier slash
    $metriques{'user'}  =~ s/([_a-zA-Z0-9]*)\/.*/$1/;

    $metriques{'mach'} = $log;
    $metriques{'mach'}  =~ s/$RESULT_PATH\///;
    $metriques{'mach'}  =~ s/$metriques{'user'}\///;
    # Tout ce qu'il y a avant le premier slash
    $metriques{'mach'}  =~ s/([_a-zA-Z0-9]*)\/.*/$1/;
    
    $metriques{'exec'}  = $log;
    #Suppression du répertoire de stockage
    $metriques{'exec'}  =~ s/$RESULT_PATH\///;
    $metriques{'exec'}  =~ s/$metriques{'user'}\///;
    $metriques{'exec'}  =~ s/$metriques{'mach'}\///;
    # Tout ce qu'il y a avant le premier slash
    $metriques{'exec'}  =~ s/([_a-zA-Z0-9]*)\/.*/$1/;
    
    
    # Récupération du nom de la matrice 
    $metriques{'matr'}  = $log;
    $metriques{'matr'}  =~ s/$RESULT_PATH\///;
    $metriques{'matr'}  =~ s/$metriques{'user'}\///;
    $metriques{'matr'}  =~ s/$metriques{'mach'}\///;
    $metriques{'matr'}  =~ s/$metriques{'exec'}\///;
    # Tout ce qu'il y a avant le premier slash
    $metriques{'matr'}  =~ s/([_a-zA-Z0-9]*)\/.*/$1/;
    
    # Récupération du cas XXP_XXT_... (suppression de tout ce qu'il y a avant le dernier /)
    $metriques{'cas'}   = $log;	    
    $metriques{'cas'}   =~ s/$RESULT_PATH\///;
    $metriques{'cas'}   =~ s/$metriques{'user'}\///;
    $metriques{'cas'}   =~ s/$metriques{'mach'}\///;
    $metriques{'cas'}   =~ s/$metriques{'exec'}\///;
    $metriques{'cas'}   =~ s/$metriques{'matr'}\///;
    # Tout ce qu'il y a avant le premier slash
    $metriques{'cas'}  =~ s/([_a-zA-Z0-9]*)\/.*/$1/;



    $metriques{'nproc'}  = $log;
    $metriques{'nproc'}  =~ s/.*\/([0-9][0-9])P_.*/$1/;
    $metriques{'nthrd'}  = $log;
    $metriques{'nthrd'}  =~ s/.*P_([0-9][0-9])T_.*/$1/;
    if ($metriques{'nproc'} =~ /[a-zA-Z]+/)
    {
	$metriques{'nproc'} = 0;
    }    
    if ($metriques{'nthrd'} =~ /[a-zA-Z]+/)
    {
	$metriques{'nthrd'} = 0;
    }

    $metriques{'nbtask'} = $metriques{'nproc'} * $metriques{'nthrd'};
    
    return (%metriques);
}

sub LogPrintMetriques {
    
    my ($log, $logerr, @listemetriques) = @_ ;
    my $str = "";
    my (%metriques) = LogGetMetriques($log, $logerr);

    foreach my $metr (@listemetriques)
    {
	#print $metr." ".$metriqueconf{$metr}{'fmt'}."\n";
	#print $metr." ".$metriques{$metr}."\n";
	$str = sprintf("%s $metriqueconf{$metr}{'fmt'}", $str, $metriques{$metr});
    }
    return ($str, %metriques);
}
    
sub LogPrintMetriquesHeader {
    
    my (@listemetriques) = @_ ;
    my $str = "";

    foreach my $metr (@listemetriques)
    {
	$str = sprintf("%s $metriqueconf{$metr}{'fmt'}", $str, $metriqueconf{$metr}{'nom'});
    }
    return ($str);
}
  
sub LogPrintMetriquesHtml {
    
    my ($log, $logerr, @listemetriques) = @_ ;
    my $str = "";
    my (%metriques) = LogGetMetriques($log, $logerr);

    foreach my $metr (@listemetriques)
    {
	#print $metr." ".$metriqueconf{$metr}{'fmt'}."\n";
	#print $metr." ".$metriques{$metr}."\n";
	$str = sprintf("%s <td>$metriqueconf{$metr}{'fmt'}</td>", $str, $metriques{$metr});
    }
    return ($str, %metriques);
}

sub LogPrintMetriquesHeaderHtml {
    
    my (@listemetriques) = @_ ;
    my $str = "";

    foreach my $metr (@listemetriques)
    {
	$str = sprintf("%s <td>$metriqueconf{$metr}{'fmt'}</td>", $str, $metriqueconf{$metr}{'nom'});
    }
    return ($str);
}

sub LogPrintMetriquesHelp {
    
    foreach my $opt ( @paramlist )
    {
	my $str = "";
	$str = sprintf("                       %-6s : %s\n", $opt, $metriqueconf{$opt}{'desc'});
	print $str;
    }
    foreach my $opt ( @paramlist2 )
    {
	my $str = "";
	$str = sprintf("                       %-6s : %s\n", $opt, $metriqueconf{$opt}{'desc'});
	print $str;
    }
}
1;
