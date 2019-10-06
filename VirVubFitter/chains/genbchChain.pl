#! /usr/local/bin/perl

die "Not enough arguments\n" if @ARGV < 1;

use File::Basename;

my $workdir = "/nfs/farm/babar/AWG36/sacco/analysis-23/workdir";
my $rootcoll = "/nfs/farm/babar/AWG36/sacco/analysis-23/workdir/1235/root/";
my $mode    = "1235";
my $run     = $ARGV[0];

system("rm $workdir/genbchChain$run");

@contdir=`echo $workdir/1235/root/$run/SP-$mode-BSemiExcl*.root | xargs ls`;
 
open TMP, ">>$workdir/genbchChain$run";

foreach $contdir(@contdir){
    $file = $contdir;
    print TMP $file;    
#     print TMP "$workdir/$contdir \n"; 
}

exit;
