#!/usr/bin/perl 

$file = "ListeJobs.txt";
$sleeptime = 30;
sub SuppressionPremiereLigne()
{
    `tail -n +2 $file > /tmp/TODO.temp && mv /tmp/TODO.temp $file`;
}

# Ouverture du fichier de log
open (M, "$file"); @st = <M>;  close M;

foreach $cmd (@st)
{

    chop $cmd;
    print $cmd."\n";
    system($cmd);
    while ( $? != 0 )
    {
	print "Submission returned value $?, Trying again in $sleeptime sec...\n";
        sleep($sleeptime);
	system($cmd);
    }

    SuppressionPremiereLigne();

}
