#!/usr/local/bin/perl

use Getopt::Long;

# This in order to have a printout of the command line on $tmpsett
my $sring="./submXFit.pl";
for($i=0;$i<@ARGV;$i++){
    $sring=$sring." ".$ARGV[$i];}

die unless &GetOptions( 'debug', \$debug, 'help', \$help, 'int', \$int, 
'norun', \$norun, 'que=s', \$que, 'flag=s', \$flag, 'cm=s', \$cm, 
'Run4', \$Run4, 'Run3', \$Run3, 'Run5', \$Run5, ,'Run6', \$Run6, 
'Run12', \$Run12, 'Run14', \$Run14, 'Run15', \$Run15, 'small!',\$small, 
'novarfit', \$novarfit, 'dirS=s',\$dirS, 'opt=i', \$opt, 'mx1d=s', \$mx1d, 
'mxscan=s', \$mxscan, 'mxbin=s', \$mxbin, 'q2cut=s', \$q2cut, 'b!', \$b, 
'd!', \$d, 'FF=i', \$ff, 'ssbar=i', \$ssbar, 'chb=i', \$chb,
'minintpur=s', \$minintpur, 'rew=i', \$rew, 'wF=s', \$wF, 'FFile=s', \$FFile, 
'unf=s', \$unf, 'unfF=s', \$unfF, 'mhiunf=f', \$mhiunf, 'mx2unf', \$mx2unf, 
'mult', \$mult, 'mpar=s', \$mpar, 'mesfitmodel=i', \$mesfitmodel, 
'depl', \$depl, 'chsys', \$chsys, 'Sys=i', \$Sys, 'SysD=i', \$SysD, 'mu=i', \$mu,
'newbin', \$newbin, 'trulept=f', \$trulept, "brecoqual", \$brecoqual, 
"rel=i", \$rel, "-empm", \$empm, 
"-mesmeancorr", \$mesmeancorr, '-neucorr', \$neucorr, '-nucomb', \$nucomb, 
'-unffit=i', \$unffit, '-comp', \$comp, '-recomputevars', \$recomputevars );

usage(0) if($help);

# warning notice if scratch area is missing
die "Scratch directory in home directory is missing. Please create ~/scra!" unless -e "$ENV{'HOME'}/scra";

if($cm){
    print "Computing Model $cm\n" if($debug);}
else {print "Computing Model CM2\n" if($debug);}

if($rel){
    print "Running on release $rel ntuples\n" if($debug);}
else{print "Running on release 22 ntuples\n" if($debug);}

if($comp){
    print "Running for data-MC comparisons\n";}

if($recomputevars){
    print "Computing vars from track and neutral lists\n";}

if($dirS){
    print "Systematics dir $dirS\n";}

$ENV{MXCUT} =        5.0;
$ENV{MXCUT} = $mxscan if ($mxscan);
print "Cutting Mx at:: $ENV{MXCUT}\n" if($debug);
$ENV{MXBIN} =        1.55;
$ENV{MXBIN} = $mxbin if ($mxbin);
print "Variable Mxbin at:: $ENV{MXBIN}\n" if($debug);
$ENV{Q2CUT} =        0.;
$ENV{Q2CUT} = $q2cut  if ($q2cut);
print "Cutting Q2 at:: $ENV{Q2CUT}\n" if($debug);
#FULL SET ON DATA
$ENV{DOBDECWEIGHT}=      0;
$ENV{DOBDECWEIGHT}=      1 if($b);
$ENV{DODDECWEIGHT}=      0;
$ENV{DODDECWEIGHT}=      1 if($d);
$ENV{DOFFWEIGHT}=      0;
$ENV{DOFFWEIGHT}=      $ff if($ff);
$ENV{DOSSBARWEIGHT} =  0;
$ENV{DOSSBARWEIGHT} =  $ssbar if($ssbar);
$ENV{DOMESMEANCORR} =  1 if($mesmeancorr);
$ENV{BTYPE}=         2;
if($chb eq "0" || $chb eq "1"){
    $ENV{BTYPE} = $chb;}
print "Running against $ENV{BTYPE} Bs\n" if($debug);

$ENV{CHLOW}=      -.5;
$ENV{CHHIGH}=      .5;
$ENV{USECB}= 1;
$ENV{GAUSSFIT}=      0;
$ENV{MNUHIGH}=     0.3;
$ENV{MNULOW}=     -100000.;
$ENV{EMPMLOW}=    -0.3;
$ENV{EMPMHIGH}=    0.3;
$ENV{PRMM2CUT}=      -3.;
$ENV{NEUCUT}= 0.05;
$ENV{LEPTCUT}=       1.;
$ENV{TRULEPTCUT}=    1.;
$ENV{TRULEPTCUT}= $trulept if ($trulept);
$ENV{LEPTTYPE}=      "2"; #0 ==  elec ; 1 == muon ; 2 ==all
$ENV{MININTPUR}=     0.; if ($minintpur > 0. && $minintpur < 1.) { $ENV{MININTPUR} = $minintpur; }
$ENV{MAXINTPUR}=     1000.;
$ENV{RUN}=           0;
$ENV{FITOPTION}=   0; #if ==3 fixes vcb fraction with fittoshape=0
if($opt){
  $ENV{FITOPTION}= $opt;
}
$ENV{DEPL} =         0;
$ENV{DEPL} = 1 if ($depl);

$ENV{REW} = 1;
$ENV{REW} = $rew if($rew);

die  "Please Select mES fit model with flag -mesfitmodel <opt>.\n  <opt> = 0 for Argus and Gauss PDF\n  1 for Argus and Crystal Ball PDF" unless $mesfitmodel;
$ENV{MESFITMODEL} = $mesfitmodel;

$ENV{BRECOQUAL} = 0;
$ENV{BRECOQUAL} = 1 if($brecoqual);

$ENV{NUCUT} = 1;
$ENV{NUCUT} = 2 if($empm);
$ENV{NUCUT} = 3 if($nucomb);

$ENV{UNFFIT} = 0;
if($unf) {
  $ENV{UNFFIT} = $unffit if $unffit;
} else {
    print "This variable is only available in unfolding mode! Use the -unf and -mhiunf flags!" if($unffit);}

$ENV{RECOMPUTEVARS} = 0;
$ENV{RECOMPUTEVARS} = 1 if($recomputevars);

$ENV{NEUCORR} = 0;
$ENV{NEUCORR} = 1 if($neucorr);

############ chains here                ##################

##### CB new reduced chains location run subVirFit.pl with -cm CM2 -small and 
#####               -dirS /u/ec/bozzi/work2/test24/VirVubFitter/chains or the chain subdir of VVF in your test release
##### CB separate run periods: use chains in the same subdir as above, cm2 only

my $runflag;

if($cm eq "CM1"){

    $dirP = "/nfs/babar/recoil/Vub_incl/root/anaQA-r00/";
    $dirP = "$dirS" if($dirS);

    $ENV{FILEDATA}=         "${dirP}datChain 1";
    $ENV{FILEVUBTOTALNRE}=  "${dirP}nreChain 1";
    $ENV{FILEVUBTOTALHYB}=  "${dirP}mixChain 1";
    $ENV{FILEVCB}=          "${dirP}genChain 1";
    $ENV{FILEVCB1}=         "${dirP}genbchChain 1";
    $ENV{FILEVCB2}=         "${dirP}genbnuChain 1";
    $ENV{FILEVUBTOTAL}=     "${dirP}mixChain 1";
    $runflag=0;

} else { #CM2

    my $fileTag = "";
    my $fileRun = "14";

    my $dirP;
    if(!$small){
	if($rel eq 14){ 
	    $dirP = "/afs/slac.stanford.edu/u/br/sacco/newsel-23/workdir";}
	else{ if($rel eq 18){ 
	    $dirP = "/afs/slac.stanford.edu/u/br/sacco/vol1/ana-32/workdir";
	    $dirP = "/nfs/farm/babar/AWG63/SemiLept/sacco/rootples/summer07/chains";
	    $dirS = "/nfs/farm/babar/AWG63/SemiLept/sacco/rootples/summer07/chains";}
	      else {#22
#Roberto's chains on which I did my spring tests
		  $dirP = "/afs/slac.stanford.edu/u/br/sacco/vol1/vvf-jan08/VirVubFitter/chains";
#Nicola's July production, just for a few tests
		  $dirP = "/afs/slac.stanford.edu/u/br/gaglio/scra3/Ana42Central/chains";
		  $dirP = "/nfs/farm/babar/AWG37/ISL/kerstin/prod/ana42/sigBcms/chains";
		  if(!$dirS){
		      $dirS = "${dirP}";}}}}
    if($small) {
	if($rel eq 14){ 
	    # default for summer05
	    $dirP="/nfs/farm/babar/AWG38/petrella/reduced/nkinfit-23";
	    # default for summer06
	    $dirP="/nfs/farm/babar/AWGsemilep01/menges/summer06/store.2/chains-1111";
	    
	    $fileTag = Red;
	    print "Running on reduced CM2 files\n" if($debug);}
	else{ if($rel eq 18){
	    $dirP="/nfs/farm/babar/AWG37/ISL/kerstin/reducedH07/def";
	    $dirS="/nfs/farm/babar/AWG37/ISL/kerstin/reducedH07/def";
	    $fileTag = Red;
	    print "Running on release 18 ntuples\n" if($debug);
	}else {#22
	    $dirP="/nfs/farm/babar/AWG37/ISL/kerstin/prod/ana42/sigBcms/reduced";
	    if(!$dirS){
		$dirS = "${dirP}";}
#	    $dirP="/nfs/farm/babar/AWG37/ISL/kerstin/prod/ana42/splitoffs/reduced";
#	    $dirS="/nfs/farm/babar/AWG37/ISL/kerstin/prod/ana42/splitoffs/reduced";
	    $fileTag = Red;
	    print "Running on release 22 ntuples\n" if($debug);
	}}
    }
#    $dirP = "$dirS" if($dirS);

    if($rel eq 14){
	$runflag=14;
    }else{ if($rel eq 18){
	$runflag=15;
    }else{
	$runflag=16;
    }}
    $runflag=4  if($Run4);
    $runflag=3  if($Run3);
    $runflag=5  if($Run5);
    $runflag=6  if($Run6);
    $runflag=12 if($Run12);
    $runflag=14 if($Run14);
    $runflag=15 if($Run15);

    my %generic; $generic{12} = "1-2"; $generic{3} = "3";   $generic{4} = "4"; $generic{14} = "1-4"; $generic{15} = "1-5"; $generic{5} = "5"; $generic{16} = "1-6"; $generic{6} = "6";
    my %signal;  
    if($rel eq 14){
	$signal{12}  = "1-3"; $signal{3}  = "1-3"; $signal{4}  = "4"; $signal{14}  = "1-4";}
    else{
	$signal{12}  = "1-2"; $signal{3}  = "3"; $signal{4}  = "4"; $signal{14}  = "1-4"; $signal{15} = "1-5"; $signal{5} = "5";$signal{16} = "1-6"; $signal{6} = "6";}


    print "Running on Run$runflag !\n" if($debug);
    print "Systematics dir $dirS\n";

    $ENV{FILEVUBTOTALNRE}=  "${dirS}/nre${fileTag}ChainRun$signal{$runflag} 1";
    $ENV{FILEVUBTOTALHYB}=  "${dirS}/mix${fileTag}ChainRun$signal{$runflag} 1";

    $ENV{FILEDATA}=         "${dirP}/dat${fileTag}ChainRun$generic{$runflag} 1";
    $ENV{FILEVCB}=          "${dirS}/gen${fileTag}ChainRun$generic{$runflag} 1";
    if($rel eq 14){
	$ENV{FILEVCB1}=         "${dirP}/gen${fileTag}bchChain 1";
	$ENV{FILEVCB2}=         "${dirP}/gen${fileTag}bnuChain 1";
    }else{
	$ENV{FILEVCB1}=         "${dirP}/gen${fileTag}bchChainRun$generic{$runflag} 1";
	$ENV{FILEVCB2}=         "${dirP}/gen${fileTag}bnuChainRun$generic{$runflag} 1";
    }
    $ENV{FILEVUBTOTAL}=     "${dirP}/mix${fileTag}Chain 1";
    $ENV{FILEVUBTRUTHHYB}=  "/nfs/farm/babar/AWG37/ISL/kerstin/reducedH07/def/VubTruth_res_july_merged.root";
    $ENV{FILEVUBTRUTHNRE}=  "/nfs/farm/babar/AWG37/ISL/kerstin/reducedH07/def/VubTruth_nonres_july_merged.root";

}

if($chsys){  
  my $tmpname=$flag;
  my %chneu=(0,'neu',1,'ch');
  for($i=0;$i<2;$i++)
    {   
      $ENV{BTYPE}=$i;
      if($flag) {
	print "$chneu{$i}:\n";
	$flag=$tmpname."_$chneu{$i}";
	$ENV{TEST}= "$flag";
      } else {$ENV{TEST}=   "MainTest_$chneu{$i}";}
      
      &scriptfit();
    }
} else {
  $ENV{TEST}=   "MainTest";
  if($flag) {$ENV{TEST}= "$flag"};
  &scriptfit(); 
}

sub scriptfit {
  
  &envfit(); # set up some defaults

  # test, if directories exist
  if (-e $ENV{DIRSETTING} && ! -d $ENV{DIRSETTING}) { die "$ENV{DIRSETTING} should be a directory!";}
  mkdir $ENV{DIRSETTING} unless -d $ENV{DIRSETTING};
  if (-e $ENV{DIRFILES} && ! -d $ENV{DIRFILES}) { die "$ENV{DIRFILES} should be a directory!";}
  mkdir $ENV{DIRFILES} unless -d $ENV{DIRFILES};

  my $tmpsett = "$ENV{DIRSETTING}/mysettings.dat_$ENV{TEST}";
  my $tmpfile = "$ENV{DIRFILES}/myfiles.dat_$ENV{TEST}";
  open OUTFILE,">$tmpfile";
  open OUTSETUP,">$tmpsett";
  system("touch  $tmpsett");system("\rm  $tmpsett");system("touch  $tmpsett");
  system("touch  $tmpfile");system("\rm  $tmpfile");system("touch  $tmpfile");

  print OUTSETUP " ndata             $ENV{NDATA} \n";
  print OUTSETUP " nvcb              $ENV{NVCB} \n";
  print OUTSETUP " nvub              $ENV{NVUB} \n";
  print OUTSETUP " nvubdata          $ENV{NVUBDATA} \n";
  print OUTSETUP " nvcbdata          $ENV{NVCBDATA} \n";
  print OUTSETUP " mxCut             $ENV{MXCUT} \n";
  print OUTSETUP " mxBin             $ENV{MXBIN} \n";
  print OUTSETUP " gauss             $ENV{GAUSSFIT} \n";
  print OUTSETUP " useCB             $ENV{USECB} \n";
  print OUTSETUP " gauss             $ENV{GAUSSFIT}\n";
  print OUTSETUP " fixMeanValue      $ENV{FIXMEAN} \n";
  print OUTSETUP " fixSigma          $ENV{FIXSIGMA} \n";
  print OUTSETUP " fixArgus1         $ENV{FIXARG1} \n";
  print OUTSETUP " fixArgus2         $ENV{FIXARG} \n";
  print OUTSETUP " fixCB1            $ENV{FIXCB1} \n";
  print OUTSETUP " fixCB2            $ENV{FIXCB2} \n";
  print OUTSETUP " leptonPCut        $ENV{LEPTCUT} \n";
  print OUTSETUP " TrueleptonPCut    $ENV{TRULEPTCUT} \n";
  print OUTSETUP " prmm2Cut          $ENV{PRMM2CUT} \n";
  print OUTSETUP " neuCut            $ENV{NEUCUT} \n";
  print OUTSETUP " mnuSqLow          $ENV{MNULOW} \n";
  print OUTSETUP " mnuSqHigh         $ENV{MNUHIGH} \n";
  print OUTSETUP " EmpmLow           $ENV{EMPMLOW} \n";
  print OUTSETUP " EmpmHigh          $ENV{EMPMHIGH} \n";
  print OUTSETUP " chLow             $ENV{CHLOW} \n";
  print OUTSETUP " chHigh            $ENV{CHHIGH} \n";
  print OUTSETUP " depl              $ENV{DEPL} \n";
  print OUTSETUP " Btype             $ENV{BTYPE} \n";
  print OUTSETUP " q2Cut             $ENV{Q2CUT} \n";
  print OUTSETUP " lepttype          $ENV{LEPTTYPE} \n";
  print OUTSETUP " mixcorr           $ENV{MIXCORR} \n";
  print OUTSETUP " doBdecreweight    $ENV{DOBDECWEIGHT} \n";
  print OUTSETUP " doFFreweight      $ENV{DOFFWEIGHT} \n";
  print OUTSETUP " doDdecreweight    $ENV{DODDECWEIGHT} \n";
  print OUTSETUP " dossbarweight   $ENV{DOSSBARWEIGHT} \n";
  print OUTSETUP " domesmeancorr     $ENV{DOMESMEANCORR} \n";
  print OUTSETUP " minintpur         $ENV{MININTPUR} \n";
  print OUTSETUP " maxintpur         $ENV{MAXINTPUR} \n";
  print OUTSETUP " run               $ENV{RUN} \n";
  print OUTSETUP " rew               $ENV{REW} \n";
  print OUTSETUP " mesfitmodel       $ENV{MESFITMODEL} \n";
  print OUTSETUP " fitOption         $ENV{FITOPTION} \n";
  print OUTSETUP " brecoqual         $ENV{BRECOQUAL} \n";
  print OUTSETUP " nucut             $ENV{NUCUT} \n";
  print OUTSETUP " neucorr           $ENV{NEUCORR} \n";
  print OUTSETUP " unffit            $ENV{UNFFIT} \n";
  print OUTSETUP " recomputevars     $ENV{RECOMPUTEVARS} \n";

  if($mxscan) {
    print OUTSETUP "  userbinning    $ENV{USERBINNING}\n";
    print OUTSETUP "  bin1           $ENV{BIN1}\n";
    print OUTSETUP "  bin2           $ENV{BIN2}\n";
    print OUTSETUP "  bin3           $ENV{BIN3}\n";
    print OUTSETUP "  bin4           $ENV{BIN4}\n";
    print OUTSETUP "  bin5           $ENV{BIN5}\n";
    print OUTSETUP "  bin6           $ENV{BIN6}\n";
    print OUTSETUP "  bin7           $ENV{BIN7}\n";
    print OUTSETUP "  bin8           $ENV{BIN8}\n";
    print OUTSETUP "  bin9           $ENV{BIN9}\n";
    print OUTSETUP "  bin10          $ENV{BIN10}\n";
   }

  close OUTSETUP;

  #Saving setup of files

  print OUTFILE " fileVubTotalnres  $ENV{FILEVUBTOTALNRE} \n";
  print OUTFILE " fileVubTotalres   $ENV{FILEVUBTOTALHYB} \n";
  print OUTFILE " fileVubTotal $ENV{FILEVUBTOTAL} \n";
  print OUTFILE " fileVcb2 $ENV{FILEVCB2} \n";
  print OUTFILE " fileVcb1 $ENV{FILEVCB1} \n";
  print OUTFILE " fileVcb  $ENV{FILEVCB} \n";
  print OUTFILE " fileData $ENV{FILEDATA} \n";
  print OUTFILE " fileVubTruthnres $ENV{FILEVUBTRUTHNRE} \n";
  print OUTFILE " fileVubTruthres $ENV{FILEVUBTRUTHHYB} \n";
  print OUTFILE "  \n";
  close OUTFILE;

  system("rm -rf $ENV{DIR}");
  system("mkdir $ENV{DIR}");

  # --- start building options for VirFit call
  my $fitflag;

  # --- get mes input parameter file
  if (!$mpar) {                              # don't override user defined filename

      if (!$cm || $cm eq "CM2") {
          $mpar = "mesparsetting_CM2.dat";   # CM2
      } elsif ($cm eq "CM1") {
	  $mpar = "mesparsetting_CM1.dat";   # CM1
      }

      $mpar = "mesparsetting.dat";           # WARNING: override CM1/CM2 choice 
  }
  $fitFlag .= " -Mes $mpar";
  

  if($unf) {
      $fitFlag = $fitFlag." -unf $unf";
      $fitFlag = $fitFlag." -U $unfF" if $unfF;
  } 
  if($mhiunf)  {$fitFlag = $fitFlag." -hiunf $mhiunf";}
  if($mx2unf) {$fitFlag = $fitFlag." -mx2unf ";}
  if($mu)  {$fitFlag = $fitFlag." -mu $mu";}
  if($mult)  {$fitFlag = $fitFlag." -mult";}
  if($Sys) {$fitFlag = $fitFlag." -Sys $Sys";}
  if($SysD) {$fitFlag = $fitFlag." -SysD $SysD";}
  if($newbin) {$fitFlag = $fitFlag." -newbin ";}
  
  if($mx1d){
    $fitFlag = $fitFlag." -mx1d $mx1d";
  }else{
    if ($newbin) {
      $fitFlag = $fitFlag." -mx1d spa/mxmulti1d.dat";
    } else {
      $fitFlag = $fitFlag." -mx1d spa/mxdef1d.dat";
    }
  }

  if($wF) {
    $fitFlag = $fitFlag." -W $wF";
  }else{
    $fitFlag = $fitFlag." -W sysWd/3dweights";
  }

  #####################
  if($cm) { 
      $fitFlag = $fitFlag. " -cm $cm";
  } else {  
      $fitFlag = $fitFlag. " -cm CM2";
  }
  if((!$cm)||($cm eq "CM2")){
    if($runflag){ 
      $fitFlag=$fitFlag." -RunFl $runflag";}
    else{
      $fitFlag=$fitFlag." -RunFl 12";}
  }else{
    $fitFlag=$fitFlag." -RunFl 0";
  }
    
  if((!$cm)||($cm eq "CM2")){
      if($rel) {
	  $fitFlag = $fitFlag. " -rel $rel";
      }else{#22
	  $fitFlag = $fitFlag. " -rel 22";
      }  
  }else{
      $fitFlag = $fitFlag. " -rel 14";
  }

  
  if($novarfit){ 
    $fitFlag = $fitFlag. " -varfit 0"; ## We are using q2, mxhad
  } else {
    $fitFlag = $fitFlag. " -varfit 1"; ## We are using q2fit, mxhadfit
  }

  if($comp){
    $fitFlag = $fitFlag. " -comp";} ## data-MC comparison mode

  $FFile = "wfermifile.dat" unless $FFile; # default
#  $FFile = "" if $ENV{REW} == 12;          # not needed for this reweighting
  $fitFlag = $fitFlag. " -FFile $FFile" if $FFile;

  my $queue = "medium";
  $queue = $que if($que);

  my $command = "../bin/$ENV{BFARCH}/mXFit -F $tmpfile -C $tmpsett -D $ENV{DIR} -P $ENV{TEST} $fitFlag";
  if (!$int) {
      print  "bsub -q $queue -o $ENV{DIR}/vubfit.out $command\n" if $debug;
      system("bsub -q $queue -o $ENV{DIR}/vubfit.out $command") unless $norun;
  } else {
      print "$command\n" if $debug;
      system("$command") unless $norun;
  }
}

sub envfit {

  # init some variables

  $ENV{DIR} = "~/scra/Ibur22_sigBcms/Ibu$ENV{TEST}/";
  $ENV{DIRSETTING} = "theset";
  $ENV{DIRFILES} = "thefiles";

  $ENV{NDATA} =     15000000;
  $ENV{NVUB} =       4000000;
  $ENV{NVCB} =      15000000;

#  $ENV{NDATA} =     15000; #just for test
#  $ENV{NVUB} =       4000; #just for test
#  $ENV{NVCB} =      15000; #just for test

#  $ENV{NDATA} =     20;
#  $ENV{NVUB} =       100;
#  $ENV{NVCB} =      100;

  $ENV{FIXMEAN} =      1;
  $ENV{FIXSIGMA} =     1;
  $ENV{FIXARG1} =      0;
  $ENV{FIXARG} =       1;
  $ENV{FIXCB1} =       1;
  $ENV{FIXCB2} =       1;
  $ENV{MIXCORR} =      2;
######  $ENV{MIXCORR} =      0;
  if($mxscan) {
    $ENV{USERBINNING} = 1;
    $ENV{BIN1} = $mxscan;
    $ENV{BIN2} = 1.9;
    $ENV{BIN3} = 2.2;
    $ENV{BIN4} = 2.5;
    $ENV{BIN5} = 2.8;
    $ENV{BIN6} = 3.1;
    $ENV{BIN7} = 3.4;
    $ENV{BIN8} = 3.7;
    $ENV{BIN9} = 4.2;
    $ENV{BIN10} = 5.;
  }
}

sub usage {
  my($exit, $message) = @_;
  
  print STDERR $message if defined $message;
  print STDERR <<INLINE_LITERAL_TEXT; #'

Usage: $0 <options> 
    Script that submits jobs for Fits.

Options:
  -help    :       print this message.
  -debug   :       Additional output for job building in this script

  -int     :       job is run interactively
  -norun   :       job is not started, but directories/input files created
  -que s   :       job is run in queue s (default xlong)
  -flag    :       Flagging output directory and files

  -cm      :       Select computing model CM1 or CM2
  -rel i   :       Selects release i=14 or i=18 or i=22 ntuples for CM2
  -dirS=s  :       Take chains from directory s
  -small   :       Run on reduced files (CM2 only; take chains from chain/ subdir)
  -Run12   :       Run only on Run1-2 data (CM2 only; take chains from chain/ subdir)
  -Run3    :       Run only on Run3 data (CM2 only; take chains from chain/ subdir)
  -Run4    :       Run only on Run4 data (CM2 only; take chains from chain/ subdir)

  -mesfitmodel  i:  Select a model to fit mES distributions: 
                       i=0 (Argus+Gaussian); 
		       i=1 (Argus+Crystal Ball); 
  -distfit :       1D fit with mx (0), q2 (1) or pplus (2) distribution
  -opt 1   :       do a 3-parameter (vub,vcb,other) 1D fit. Default is 2-parameter. 
  -mx1d    :       choose mx favourite binning (1d only)

  -b       :       B reweighting
  -d       :       D reweighting
  -ssbar i :       ssbar reweighting (non-res vub: 1=-30% , 2=+30%)
  -mxscan  :       Mx upper cut
  -mxbin   :       Mx bin scan
  -q2cut   :       Q2 upper cut
  -q2bin   :       Q2 bin scan
  -novarfit:       use measured variables instead of fitted ones
  -minintpur:      minimum integrated purity (default: 0.)
  -rew     :       Apply vub reweighting 1: 1D, 3:3D, 10 a la dominique, 11 a la dominique with latest mx binning
  -wF      :       Applies Theo model reweighting
  -FFile   :       Use indicated file instead of wfermifile.dat

  -unf     :       Equidistant binning (-unf numbins)
  -unffit i:       Variable to bin in
  -comp    :       Data-MC comparison mode, to be used with unffit etc.
  -mhiunf  :       Endpoint for equidistant binning
  -mx2unf  :       Prepare the spectra for mx2 unfolding

  -mpar=s  :       File with mES starting parameters
  -tmodel  :       use the new pdfs from Thorsten
  -depl    :       Run on depleted sample
  -chsys   :       Used in conjuction with -Sys, submits 2 jobs for charged and neutral Bs with the same
                   random number seed (for B and D systematics)
Examples:
$0 -int -d -b -wF sysWd/3dweights -debug -flag ciappi -cm CM1

The following reproduces the prl result as closely as possible
(binned mES fits, 3-parameters mX fit): 
$0 -int -d -b -binmes -rew 1 -wf sysWd/ -opt 1 -debug -flag prlresult

The following will fit prl-like, but with unbinned mES fits:
$0 -int -d -b -rew 1 -wf sysWd/ -opt 1 -debug -flag unbin_prlresult

The following will fit prl-like, unbinned mES fits, 2-parameter mX fit: 
$0 -int -d -b -rew 1 -wf sysWd/ -debug -flag unbin_prlresult_2par

1d mx fit:
$0 -int -d -b -cm CM2 -debug -novarfit -rew 11 -wF sysWd/belle3d_sum04_0.inc -FFile wfermifile.dat -Sun -flag 1DMR14def

1d pplus fit:
$0 -int -d -b -cm CM2 -debug -novarfit -rew 11 -wF sysWd/belle3d_sum04_0.inc -FFile wfermifile.dat -Sun -distfit 2 -flag 1DPR14def

2d mx/q2 fit:
$0 -int -d -b -cm CM2 -debug -novarfit -rew 11 -wF sysWd/belle3d_sum04_0.inc -FFile wfermifile.dat -comb -Sun -mxbin 1.7 -q2bin 8 -flag 2DR14def

1d mx fit with fixed signal/peakingBG ratios, MC counting, fixed (D**+nonSL)/(D+D*) ratio, new PDFs, small ntuples: 
$0 -int -b -d -cm CM2 -depl -debug -Sun -novarfit -small -tmodel -notunbinmes=3 -rew 12 -wF sysWd/HFAGCombWin06_3d_0.txt -fixSBratio corrratiosigpeakmx_newdatadepl.txt -countmc -dssRatio 0.551496 -flag 1DMR14new


Comments to <azzolini\@slac.stanford.edu>
            <menges\@slac.stanford.edu>
            <petrella\@slac.stanford.edu>
            <sacco\@slac.stanford.edu>.
INLINE_LITERAL_TEXT
#'
  exit($exit) if defined $exit;
}

