#! /usr/local/bin/perl

die "Not enough arguments\n" if @ARGV < 1;

use File::Basename;

my $workdir = "/nfs/farm/babar/AWG36/sacco/analysis-23/workdir";
my $rootcoll = "/nfs/farm/babar/AWG36/sacco/analysis-23/workdir/1237/root/";
my $mode    = "1237";
my $run     = $ARGV[0];

system("rm $workdir/genbnuChain$run");

@contdir=`echo $rootcoll/$run/SP-$mode-BSemiExcl*.root | xargs ls`;
 
open TMP, ">>$workdir/genbnuChain$run";

foreach $contdir(@contdir){
    $file = $contdir;
    print TMP $file;    
}

exit;
