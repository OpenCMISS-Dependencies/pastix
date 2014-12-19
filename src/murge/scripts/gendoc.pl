#!/usr/bin/env perl

###############################################################################
#
#   Script: gendoc.pl
#
###############################################################################
#
#   A script that generate the doc from Murge sources.
#
#   Usage:
#    > ./gendoc.pl -d <NaturalDocsFileDirectory> -n <NaturalDocsExecutable>
#    >    Generate the documentation from Murge sources.
#    >     -h          Shows this help
#    >     -d <dir>    Sources Directory         (Default : ../)
#    >     -i <dir>    Menu/Images Directory     (Default : ./NaturalDocs)
#    >     -n <exe>    NaturalDocs Executable    (Default : NaturalDocs)
#    >     -o <dir>    Defines output directory  (Default : ./doc)
#    >     -p <dir>    Defines project directory (Default : ./pdoc)
#
#   Authors:
#     Xavier Lacoste - lacoste@labri.fr
#
###############################################################################
use POSIX;
use strict;
use Getopt::Std;

###############################################################################
# Group: Variables

#
#   string: naturaldocs
#     Path to root of murge sources.
#
#   string: dirsrc
#     Path to root of murge sources.
#
#   string: dirnd
#     Path to directory containing images.
#
#   string: output
#     Path to output directory.
#
#   string: project
#     Path to project directory.
#
#   hash: opts
#     Options given to the script

my $naturaldocs = "NaturalDocs";
my $dirsrc      = '..';
my $dirnd       = './NaturalDocs';
my $output      = "doc";
my $project     = "pdoc";
my %opts;

###############################################################################
#   Group: Functions

#
#   Function: Usage
#
#   Prints usage.
#
sub Usage {

    print "./gendoc.pl \n";
    print "  Generate the documentation from Murge sources.\n";
    print "   -h          Shows this help\n";
    print "   -d <dir>    Sources Directory         (Default : ../)\n";
    print "   -i <dir>    Menu/Images Directory     (Default : ./NaturalDocs)\n";
    print "   -n <exe>    NaturalDocs Executable    (Default : NaturalDocs) \n";
    print "   -o <dir>    Defines output directory  (Default : ./doc)\n";
    print "   -p <dir>    Defines project directory (Default : ./pdoc)\n";
}
 
#
#   Function: GenDoc
#   
#   Generates Murge documentation directory.
# 
sub GenDoc {
  my $cmd;
  $cmd  = "$naturaldocs -i $dirsrc";
  $cmd .= " -o HTML $output -p $project";
  $cmd .= " -img $dirnd";
    
  #creatind NaturalDocs Directories
  if ( ! -d $output )
    {
      `mkdir -p $output`;
    }
  if ( ! -d $project )
    {
      `mkdir -p $project`;
      #system "cp ".$directory."/scripts/NaturalDocs/Topics.txt ".$project."/Topics.txt ";
    }
  else
    {
	#print "Dossiers présents\n";
    }
    #Running Natural Docs
  system "cp ".$dirnd."/Menu.txt ".$project."/Menu.txt ";
  print $cmd."\n";
  system "$cmd";

}


getopts("hn:d:o:p:i:",\%opts);

if ( defined $opts{o} ){
  $output = $opts{o};
}
if ( defined $opts{p} ){
  $project = $opts{p};
}
if ( defined $opts{d} ){
  $dirsrc = $opts{d};
}
if ( defined $opts{i} ){
  $dirnd = $opts{i};
}
if ( defined $opts{n} ){
  $naturaldocs = $opts{n};
} 
if ( defined $opts{h} ){
    Usage();
    exit;
}

GenDoc();
