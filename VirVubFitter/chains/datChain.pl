#! /usr/local/bin/perl

use File::Basename;

die "Not enough arguments\n" if @ARGV < 1;

my $workdir = "/nfs/farm/babar/AWG36/sacco/analysis-23/workdir";
my $mode    = "";
my $run     = $ARGV[0];
my $rootcoll = "/nfs/farm/babar/AWG36/sacco/analysis-23/workdir/data/$run/root";

system("rm $workdir/datChain$run");

@contdir=`ls $rootcoll/BSemiExcl*.root`;
 
open TMP, ">>$workdir/datChain$run";

foreach $contdir(@contdir){
    $file = $contdir;
    print TMP $file;    
}

exit;
