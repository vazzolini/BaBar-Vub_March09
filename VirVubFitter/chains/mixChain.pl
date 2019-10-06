#! /usr/local/bin/perl
#usage: mixChain.pl 1 for running over run1

#die "Not enough arguments\n" if @ARGV < 1;

my $workdir = "/nfs/farm/babar/AWG36/sacco/analysis-23/workdir";
my $rootcoll = "/nfs/farm/babar/AWG36/sacco/analysis-23/workdir/3618/root";
my $run = $ARGV[0];

print  $ARGV[0]; print "\n";

system("rm $workdir/mixChain");

system("ls $rootcoll/*3618*BSemiExcl*.root");

@contdir=`ls $rootcoll/*3618*BSemiExcl*.root`;

open TMP, ">>$workdir/mixChain";

foreach $contdir(@contdir){
    $file = $contdir;
    print TMP $file;
}


exit;
