#! /usr/local/bin/perl
#usage: 

#die "Not enough arguments\n" if @ARGV < 1;

my $workdir = "/nfs/farm/babar/AWG36/sacco/analysis-23/workdir";
my $rootcoll = "/nfs/farm/babar/AWG36/sacco/analysis-23/workdir/2575/root";
#my $run = $ARGV[0];

#print  $ARGV[0]; print "\n";

system("rm $workdir/nreChain$run");

#system("ls $rootcoll/*2575*BSemiExcl*.root");

@contdir=`ls $rootcoll/*2575*BSemiExcl*.root`;

open TMP, ">>$workdir/nreChain";

foreach $contdir(@contdir){
    $file = $contdir;
    print TMP $file;
}


exit;
