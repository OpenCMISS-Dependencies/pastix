#!/usr/bin/perl

use strict;
use FindBin qw($Bin);
use lib "$Bin";
use Getopt::Std;
use Modules::Common;
my %opts;

# Variables globales 
my $tracefilefmt = "trace*.[0-9][0-9]";
my $t0 = 0.;
my @t0byproc;
my $max_task = 0;
my @tps_by_proc;
my @old_state;
my @old_time;
my @print_state;
my @print_time;
my @nbthrdbyproc;
my @indice_thread;
my $nbproc = 0;
my @traces;
my @memproc;

my %ind = (
    "time" => 0,
    "proc" => 1,
    "thrd" => 2,
    "lvl"  => 3,
    "type" => 4,
    "stte" => 5,
    "id"   => 6,
    "size" => 7,
    "idbbl"=> 7,
    "fcand"=> 8,
    "lcand"=> 9,
    "cand" => 10,
    "tskid"=> 11,
    "local"=> 12,
    );
my $username=`id -un`;
chop $username;

# On doit spécifier le nom du fichier de sortie
####################################################
#                Usage                             #
####################################################

sub Usage {

    print "MakeTrace.pl [ -h ] [ -b ] [ -d | -c ] [ -v ] output\n";
    print "   -h   Affiche cette aide\n"; 
    print "   -b   Calcule la trace de blend a partir de traceBlend.trf\n";
    print "   -v   Remplace les . par des , pour les problèmes de langue dans Paje\n";
    print "  Par defaut, la trace est coloriée en fonction du nombre de candidat\n";
    print "   -d   Colorie en fonction de l'indice de la bulle propriétaire\n";
    print "   -c   Colorie en fonction du thread candidat dans la simulation statique\n";
    print "   -l   Colorie en fonction de l'emplacement des donnees\n";
    print "   -L   Colorie en rouge/vert suivant que la donnee est locale ou non\n";
    print "   -s   Colorie en fonction du type de tache\n";
}

sub MaZtimestamp  
{
    my ($filein, $fileout) = @_;
    my $first;

    open(M, $filein);
    open(N, ">".$fileout);

    foreach(<M>)
    {
	my ($stamp, $a) = split /\s+/, $_, 2;
	if (defined $first) {
	    $stamp = $stamp - $first;
	} else {
	    $first = $stamp;
	    $stamp = 0;
	}
	
	printf N "%09.9f %s", $stamp, $a;
    }

    close(M);
    close(N);
}

# Tri des fichiers de trace par date 
sub Trace_Tri {
    
    my ($tracefile) = @_;
    
    # Suppression du premier timestamp dans chacun des fichier de trace 
    if (0) # Inutile
    {
	foreach my $trace (@traces)
	{
	    chomp $trace;
	    #Tri des fichiers de traces dans un seul fichier
	    MaZtimestamp("$trace", "/tmp/tmp.$username..$trace");
	}    
    }

    # Tri de tous les fichiers dans un fichier temporaire 
    # en fonction du timsestamp de chaque evenement
    system("sort -n $tracefilefmt -s -k 1 > $tracefile") == 0 or die("ERROR : Pb de tri des traces");

    # Version ou on supprime le premier timestamp global
    # a tout le monde
    if ( 0 ) 
    {
	MaZtimestamp("$tracefile", "/tmp/tracetmp2.$username");
	system("mv /tmp/tracetmp2.$username $tracefile");
    }
}

sub Trace_GetArchi
{
    my ($tracefile) = @_;
    my $i = 0;
    open (M, $tracefile)  or die("ERROR : Pb d'ouverture du fichier de trace (GetArchi)"); 

    my $time;
    my $proc;
    my $thrd;
    my $lvl;
    my $type;
    my $stte;
    my $id;

    foreach (<M>)
    {
	my @ligne = split ;

	$time = $ligne[$ind{"time"}];
	$proc = $ligne[$ind{"proc"}];
	$thrd = $ligne[$ind{"thrd"}];
	$lvl  = $ligne[$ind{"lvl"}];
	$type = $ligne[$ind{"type"}];
	$stte = $ligne[$ind{"stte"}];
	$id   = $ligne[$ind{"id"}];

	#next if ( ( ($type == 0) || ($type == 1)) && ($lvl != 1) );
	next if ($type != -1);
	
	$nbproc              = $proc if ($proc > $nbproc);
	$nbthrdbyproc[$proc] = 0     if (!defined($nbthrdbyproc[$proc]));
	$nbthrdbyproc[$proc] = $thrd if ($thrd > $nbthrdbyproc[$proc]);
	$max_task            = $id   if ($id   > $max_task);
	
	# Récupération du temps initial 
	$t0              = $time if ($i == 0);
	$t0byproc[$proc] = $time if (!defined($t0byproc[$proc]));
	$i++;
    }
    $nbproc++;
    close M;
    print $nbproc.' '.$nbthrdbyproc[0]."\n";
}

sub Trace_Gen
{

    my ($output, $tracefile) = @_;
    my $time; # Timestamp
    my $proc; # Numero du Proc
    my $thrd; # Numero du thread
    my $lvl;  # Niveau de trace
    my $type; # Type de trace
    my $stte; # Etat de la trace
    my $id;   # id de la trace
    my $gid;  # Indice global du thread
    my $state;
    my $nbprocend = 0;

    open(OUT, ">".$output) or die("ERROR : Pb ouverture du fichier de sortie de traces");

    # Recopie de l'entete commune a chaque trace
    open(N,   "$SCRIPT_PATH/head.trace") or die("ERROR : Pb d'ouvertuure du fichier d'entete");
    foreach (<N>)
    {
	print OUT;
    }
    close N;

    # Creation de la structure de la machine
    $indice_thread[0] = 0;
    for(my $i=0; $i < $nbproc; $i++)
    {
	# Nouveau noeud (Processus MPI)
	$indice_thread[$i+1] = $indice_thread[$i];
	print OUT '7 0.000000 C_N'.$i.' CT_Node C_Net0 "Node'.$i.'"'."\n";
	for(my $j=0; $j < $nbthrdbyproc[$i]+1; $j++)
	{
	    # Nouveau Processeur (Thread)
	    my $k = $indice_thread[$i+1];
	    print OUT '7 0.000000 C_P'.$k.' CT_Proc C_N'.$i.' "Proc'.$k.'"'."\n";
	    $indice_thread[$i+1]++;
	}
    }

    # Conversion au format Paje du fichier de trace
    open (M, $tracefile) or die("ERROR : Pb d'ouverture du fichier de trace (Gen)");
    foreach (<M>)
    {
	my @ligne   = split ;

	$time = $ligne[$ind{"time"}];
	$proc = $ligne[$ind{"proc"}];
	$thrd = $ligne[$ind{"thrd"}];
	$lvl  = $ligne[$ind{"lvl"}];
	$type = $ligne[$ind{"type"}];
	$stte = $ligne[$ind{"stte"}];
	$id   = $ligne[$ind{"id"}];
	$gid  = $indice_thread[$proc] + $thrd;
	
	if (($type == 1) && ($stte == 8))
	{
	    printf OUT "8 %f C_P%d CT_Proc\n", $time, $gid;
	    $tps_by_proc[$gid] = $time;
	    next;
	}		

	#$nbprocend++ if ( ($type == 1) && ($stte == 5) );
	#last if ($nbprocend == $nbproc);

	# Etat 
	if ( $type == -1 ) # Fin de traces
	{
	    #print $time.' '.$proc.' '.$thrd.' '.$lvl.' '.$type.' '.$stte.' '.$id."\n";
	    
	    next if ( $stte == 0 );
	    printf OUT "8 %f C_P%d CT_Proc\n", $time, $gid  if ( $thrd != -1 );
	    printf OUT "8 %f C_N%d CT_Node\n", $time, $proc if ( $thrd == -1 );
	    next;
	}

	if ($type == 0) # Debut de tache
	{
		# Etat du noeud 
	    if ($lvl == 0) 
	    {
		printf OUT "10 %f ST_NodeState C_N%d SN_%d\n", $time, $proc, $stte;
	    }
	    elsif ($lvl == 1)
	    {
		$state = $stte;
		if ( !(defined($opts{s})) && ( $stte > 11 && $stte < 20 ))  # Si on est dans l'etat de calcul on change l'etat
		{
		    $state =(($ligne[$ind{"lcand"}]-$ligne[$ind{"fcand"}]+1)%30 + 20);
		    $state = ($ligne[$ind{"idbbl"}]%30 + 20) if ( defined($opts{d}) );
		    $state = ($ligne[$ind{"cand" }]%30 + 20) if ( defined($opts{c}) );
		    $state = ($ligne[$ind{"local"}]%30 + 20) if ( defined($opts{l}) );
		    if ( defined($opts{L}) ) {
			if ($ligne[$ind{"local"}] != $gid) {
			    $state = 21; 
			} else {
			    $state = 20;
			}
		    }
		}
		
		printf OUT "10 %f ST_ProcState C_P%d S_%d\n", $time, $gid, $state;
		$tps_by_proc[$gid] = $time;
	    } 
	    else
	    {
		# On passe la ligne si il s'agit d'une tache de niveau superieur a 1
		next;
	    }
	}
	# Fin de tache
	elsif ($type == 1)
	{
	    if ($lvl == 1 && ($stte > 13 && $stte < 18))
	    {
# 		$state = $stte;
# 		$state = ($id%30 + 20) if ( $stte == 12); # Si on est dans l'etat de calcul on change l'etat
# 		$state = ($id%30 + 20) if ( $stte > 13 && $stte < 18); # Trace Blend
		
		printf OUT "10 %f ST_ProcState C_P%d S_0\n", $time, $gid;
		$tps_by_proc[$gid] = $time;
	    } 
	    next;
	}
	# Envoi
	elsif ($type == 2)
	{
	    my $size =  $ligne[$ind{"size"}];
	    printf OUT "42 %f L_%d C_Net0 C_P%d %d %s\n", $time, $stte, $gid, $id, $proc."_".$lvl."_".$id."_".$size;
	}
	    # Reception
	elsif ($type == 3)
	{
	    my $size =  $ligne[$ind{"size"}];
	    printf OUT "43 %f L_%d C_Net0 C_P%d %d %s\n", $time, $stte, $gid, $id, $lvl."_".$proc."_".$id."_".$size;
	}
	# Memoire
	elsif ($type == 4)
	{
	    $id = int($id);
	    if ( !defined($memproc[$proc])
		 || ($id != $memproc[$proc]))
	    {
		$memproc[$proc] = $id;
		printf OUT "51 %f V_Mem C_N%d %f\n", $time, $proc, $id;
	    }
	}
	else
	{
	    print STDERR "ERROR(0) : ";
	    for (my $i = 0; $i <= $#ligne; $i++)
	    {
		print STDERR $ligne[$i]." ";
	    }
	    print STDERR "\n";
	    exit;
	}
    }
    close M;
    
}

sub Trace_Pt2Vg
{
    my ($file) = @_;
    system('sed -i \'s/\./,/g\' '.$file);
}

#################################################################################
#
#                   Main
#
#################################################################################

getopts("clLbdsv",\%opts);

# On doit préciser le fichier de sortie
if ($#ARGV < 0 || defined $opts{h})
{
    Usage();
    exit;
}

if (defined $opts{b})
{
    $tracefilefmt = 'traceBlend.trf';
}
@traces = `ls $tracefilefmt`;

if(@traces == 0)
{
    print STDERR "ERROR : Il n'y a pas de fichiers de traces\n";
    exit;
}

Trace_Tri("/tmp/tracetmp.$username");
Trace_GetArchi("/tmp/tracetmp.$username");
Trace_Gen($ARGV[0], "/tmp/tracetmp.$username");
if (defined $Getopt::Std::opts{v})
{
    Trace_Pt2Vg($ARGV[0]);
}


