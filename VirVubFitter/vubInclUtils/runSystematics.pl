#! /usr/bin/perl -w

# $Id: runSystematics.pl,v 1.3 2006/07/24 13:54:07 menges Exp $
# first version by Wolfgang Menges
#
# script to run systematics 

require 5.000;

use Getopt::Long;

$help = <<__HELP__; 
This program runs all systematics for the three kinematic variables. 

Options are:

 - running only on a selected variable
 - running on a differnt SF
 - running only part of the systematics (smearing)

See usage for options ($0 --usage).
__HELP__

$usage = <<__USAGE__;
usage: $0 [-hvqna] [OPTIONS] mode run [pattern]

OPTIONS:
            -u, --usage: this message
             -h, --help: help message
          -v, --version: print version to standard output and exit
            -q, --quiet: shut up
           -n, --noexec: do not run an external program

                    -mx: fit Mx (default: all three variables)
                   -mx6: fit Mx with bin edge at 1.625 GeV (default: don't run)
                   -mx7: fit Mx with bin edge at 1.7 GeV (default: don't run)
                 -pplus: fit P+ (default: all three variables)
                  -mxq2: fit Mx/Q2 (default: all three variables)


                     -b: run B systematics (default: all)
                     -d: run D systematics (default: all)
                     -f: run FF systematics (default: all)
                 -other: run all except B, D, FF systematics (default: all)

                  -SF=s: choose SF weight set

Example: $0 [-hvqna] -SF=babar05 
__USAGE__

@options = ('help|h', 'version|v', 'debug|g', 'noexec|n','quiet|q', 'usage|u',
	    'mx', 'mx6', 'mx7', 'pplus', 'mxq2', 'b', 'd', 'f', 'other|o', 'SF=s');

die "$usage" unless &GetOptions(@options);
die "This is $0: v 0.9 2006/07/24 \n" if $opt_version;
die $help if $opt_help;
die $usage if $opt_usage;

# parameters to fit (default all: mx, p+, mx/q2)
my $fitPar = "";
   $fitPar .= "Mx" if defined $opt_mx && $opt_mx;
   $fitPar .= "M6x" if defined $opt_mx6 && $opt_mx6;
   $fitPar .= "M7x" if defined $opt_mx7 && $opt_mx7;
   $fitPar .= "Pp" if $opt_pplus;
   $fitPar .= "MxQ2" if $opt_mxq2;
   $fitPar = "MxPpMxQ2" if $fitPar eq "";

# set up command with common things
my $command = "./subVirFit.pl -cm CM2 -debug -novarfit -que kanga -small";

# set up options
my $runAll = "";
   $runAll .= "O" if defined $opt_other && $opt_other;
   $runAll .= "B" if defined $opt_b && $opt_b;
   $runAll .= "D" if defined $opt_d && $opt_d;
   $runAll .= "F" if defined $opt_f && $opt_f;
   $runAll = "BDFO" if $runAll eq "";

print "Fitting: $runAll\n";

my $rewOptions = "-d -b -FF=1";
my $fitOptions = "-fixshape -tmodel -notunbinmes=3";

my $sf = "hfagcombwin06"; # default
   $sf = "babar05" if defined $opt_SF && $opt_SF eq "babar05";
   $sf = "belle04" if defined $opt_SF && $opt_SF eq "belle04";

my ($SFOptions, $SFTag, $SFPstop, $BFstart);
if ($sf eq "hfagcombwin06") { $SFOptions = $SFTag = "HFAGCombWin06"; $SFPstop = 12; $BFstart = 13; };
if ($sf eq "babar05")       { $SFOptions = "babar05"; $SFTag = "Babar05"; $SFPstop = 11; $BFstart = 12; };
if ($sf eq "belle04")       { $SFOptions = "belle3d_sum04"; $SFTag = "Belle04"; $SFPstop = 8; $BFstart = 11; };

my $hybridOptions = "-rew 12 -wF sysWd/${SFTag}_3d_0.txt";

# --- run default
my $tag = '1111';
wrapper($command, $fitPar, "Su${SFTag}.sys${tag}", $fitOptions, $rewOptions, $hybridOptions) if $runAll =~ /O/;

# --- tracking, neutral, pid sys
foreach $tag (qw/ tracking neutralEff neutralCor pidEkill pidMkill pidKkill pidEmis pidMmis pidKmis klong/ ) {

    my $sysOptions = "-nosmall -dirS=/nfs/farm/babar/AWGsemilep01/menges/summer06/store.3/chains-sys.$tag -que xlong";

    wrapper($command, $fitPar, "Su${SFTag}.sys${tag}", $fitOptions, $rewOptions, $hybridOptions, $sysOptions) if $runAll =~ /O/;
}


# --- form factor reweighting

foreach $seq (100 .. 199) {
    my $sysOptions = "-nod -Sys=${seq}"; my $tag = "B${seq}";

    wrapper($command, $fitPar, "Su${SFTag}.sys${tag}", $fitOptions, $rewOptions, $hybridOptions, $sysOptions) if $runAll =~ /B/;
}
foreach $seq (100 .. 199) {
    my $sysOptions = "-nob -Sys=${seq}"; my $tag = "D${seq}";

    wrapper($command, $fitPar, "Su${SFTag}.sys${tag}", $fitOptions, $rewOptions, $hybridOptions, $sysOptions) if $runAll =~ /D/;
}

foreach $seq (0) {
# foreach $seq (100 .. 199) {
    my $sysOptions = "-FF=${seq}"; my $tag = "FF${seq}";

    wrapper($command, $fitPar, "Su${SFTag}.sys${tag}", $fitOptions, $rewOptions, $hybridOptions, $sysOptions) if $runAll =~ /F/;
}

# -- ssbar popping
{
    my $sysOptions = "-ssbar=1"; my $tag = "SSbar1";
    wrapper($command, $fitPar, "Su${SFTag}.sys${tag}", $fitOptions, $rewOptions, $hybridOptions, $sysOptions) if $runAll =~ /O/;

    $sysOptions = "-ssbar=2"; $tag = "SSbar2";
    wrapper($command, $fitPar, "Su${SFTag}.sys${tag}", $fitOptions, $rewOptions, $hybridOptions, $sysOptions) if $runAll =~ /O/;
}

# --- mes parameters
foreach $seq (1 .. 30) {
    my $sysOptions = "-me $seq"; my $tag = "ME$seq";

    wrapper($command, $fitPar, "Su${SFTag}.sys${tag}", $fitOptions, $rewOptions, $hybridOptions, $sysOptions) if $runAll =~ /O/;
}

# --- SF ellipse
foreach $seq (2 .. $SFPstop) {
    my $hybridOptions = "-rew 12 -wF sysWd/${SFTag}_3d_${seq}.txt"; my $tag = "SFP$seq";

    wrapper($command, $fitPar, "Su${SFTag}.sys${tag}", $fitOptions, $rewOptions, $hybridOptions) if $runAll =~ /O/;
}

# --- SF ellipse
foreach $seq (qw/gau roman/) {
    my $hybridOptions = "-rew 12 -wF sysWd/${SFTag}_3d_0_${seq}.txt"; my $tag = "SFF$tag";

    wrapper($command, $fitPar, "Su${SFTag}.sys${tag}", $fitOptions, $rewOptions, $hybridOptions) if $runAll =~ /O/;
}

# --- BF 
foreach $seq (1 .. 13) {
    my $calc = $seq + $BFstart - 1; if ($sf eq "belle04" && $seq > 5) { $calc += 2; };
    my $hybridOptions = "-rew 12 -wF sysWd/${SFTag}_3d_${calc}.txt"; my $tag = "BF$seq";

    wrapper($command, $fitPar, "Su${SFTag}.sys${tag}", $fitOptions, $rewOptions, $hybridOptions) if $runAll =~ /O/;
}

# --- peaking background
{
    my $fitOptions = "-fixshape -que xlong"; my $tag = "Peak";

    wrapper($command, $fitPar, "Su${SFTag}.sys${tag}", $fitOptions, $rewOptions, $hybridOptions) if $runAll =~ /O/;
}

exit 0;

sub wrapper {
    my $mxOptions   = "-Sun";
    my $mx6Options  = "-Sun -mxbin 1.625";
    my $mx7Options  = "-Sun -mxbin 1.7";
    my $ppOptions   = "-Sun -distfit 2";
    my $mxq2Options = "-Sun -mxbin 1.7 -q2bin 8";

    my $command = $_[0]; shift;
    my $fitPar  = $_[0]; shift;
    my $flag    = $_[0]; shift;
    my @fitOptions = @_;

    my_run("$command $mxOptions   @fitOptions -flag 1DM${flag}")  if $fitPar =~ /Mx/;
    my_run("$command $mx6Options  @fitOptions -flag 1DM6${flag}") if $fitPar =~ /M6x/;
    my_run("$command $mx7Options  @fitOptions -flag 1DM7${flag}") if $fitPar =~ /M7x/;
    my_run("$command $ppOptions   @fitOptions -flag 1DP${flag}")  if $fitPar =~ /Pp/;
    my_run("$command $mxq2Options @fitOptions -flag 2D${flag}")   if $fitPar =~ /Q2/;
}

sub my_run {
    my($cmd) = @_;
    print "start_all: $cmd\n" unless $opt_quiet;
    system $cmd unless $opt_noexec;
}

