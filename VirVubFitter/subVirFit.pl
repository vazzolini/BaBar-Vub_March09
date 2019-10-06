#!/usr/local/bin/perl

use Getopt::Long;

# This in order to have a printout of the command line on $tmpsett
my $sring="./subVirFit.pl";
for($i=0;$i<@ARGV;$i++){
    $sring=$sring." ".$ARGV[$i];}

die unless &GetOptions( 'debug', \$debug, 'help', \$help,
	     'int', \$int, 'norun', \$norun, 'que=s', \$que, 'flag=s', \$flag,
	     'cm=s', \$cm, 'Run4', \$Run4, 'Run3', \$Run3, 'Run5', \$Run5, 'Run6', \$Run6, 'Run12', \$Run12, 'Run14', \$Run14, 'Run15', \$Run15, 'Run16', \$Run16, 'small!',\$small, 'novarfit', \$novarfit,
	     'dirS=s',\$dirS, 'sdDir=s',\$sdDir, 'rdDir=s',\$rdDir,
	     'distfit=i', \$q2fit, 'opt=i', \$opt, 'vcb=s', \$vcb, 'oth=s', \$oth, 'cat', \$cat, 'comb', \$comb, 
	     'mx1d=s', \$mx1d, 'pplus1d=s', \$pplus1d, 'q2V=s', \$q2vec, 'mxV=s', \$mxvec, 
	     'mxscan=s', \$mxscan, 'mxbin=s', \$mxbin, 'q2cut=s', \$q2cut, 'q2bin=s', \$q2bin, 'ppluscut=s', \$pplusbin, 'ppluscut=s', \$pplusbin, 'mnscan=s', \$mnscan,
	     'b!', \$b, 'd!', \$d, 'exclFF=i', \$exclff, 'FF=i', \$ff, 'ssbar=i', \$ssbar, 'me=i', \$me, 'chb=i', \$chb,
	     'nore', \$nore, 're', \$re, 
	     'ewpwlo=s', \$ewpwlo, 'ewpwhi=s', \$ewpwhi, 'csilo=s', \$csilo, 'csihi=s', \$csihi, 
	     'minintpur=s', \$minintpur, 'xlo=s', \$xlo, 'xhi=s', \$xhi, 'wlo=s', \$wlo, 'whi=s', \$whi, 'qlo=s', \$qlo, 'qhi=s', \$qhi,
	     'rew=i', \$rew, 'wF=s', \$wF, 'FFile=s', \$FFile, 'Sun', \$Sun, 'ckBauer', \$ckBauer, 'mulf=s', \$mulf, 'fermi=s', \$fermi,
	     'unf=s', \$unf, 'unfF=s', \$unfF, 'mhiunf=f', \$mhiunf, 'mx2unf', \$mx2unf, 'mult', \$mult,
	     'mpar=s', \$mpar, 'mesfitmodel=i', \$mesfitmodel, 'binmes', \$binmes, 'notunbinmes=i', \$notunbinmes, 'fixshape', \$fixsh, 'depl', \$depl, 'fitdss=i', \$ftdss, 
             'chsys', \$chsys,'pref=s', \$pref, 'Sys=i', \$Sys, 'mu=i', \$mu, 'wisys=s', \$wisys, 'bincut=i', \$bincut, 'newbin', \$newbin, 'trulept=f', \$trulept, 'lptyp=i', \$lptyp, 
	     'fixSBratio=s', \$fixSBratio, 'countmc', \$countmc, 'dssRatio=f', \$dssRatio, 'fitmc', \$fitmc, 'SPseed=i', \$SPseed, 'dssFile=s', \$dssFile, 
	     "brecoqual", \$brecoqual, "smallstatcorr", \$smallstatcorr, "rel=i", \$rel, "-empm", \$empm, "-mesmeancorr", \$mesmeancorr, 'subtractpkgbkg' ,\$subtractpkgbkg,
			'fixcorrratio',\$fixcorrratio, 'subtractpkgbkgmx', \$subtractpkgbkgmx, 'fitallmesrange', \$fitallmesrange, 'service=i', \$service, 
			'cascade=f', \$cascade, 'effdfn=s', \$effdfn, 'bsys=i' ,\$bsys, 'thecomparison', \$thecomparison, 'savepdftree', \$savepdftree,
			'readpdftree=s', \$readpdftree, 'toyhistogrames', \$toyhistogrames);

usage(0) if($help);

# warning notice if scratch area is missing
die "Scratch directory in home directory is missing. Please create ~/scra!" unless -e "$ENV{'HOME'}/scra";

if($cm){
    print "Computing Model $cm\n" if($debug);}
else {print "Computing Model CM2\n" if($debug);}

if($rel){
    print "Running on release $rel ntuples\n" if($debug);}

$ENV{MXCUT} =        5.0;
$ENV{MXCUT} = $mxscan if ($mxscan);
print "Cutting Mx at:: $ENV{MXCUT}\n" if($debug);
$ENV{MXBIN} =        1.55;
$ENV{MXBIN} = $mxbin if ($mxbin);
print "Variable Mxbin at:: $ENV{MXBIN}\n" if($debug);
$ENV{Q2CUT} =        0.;
$ENV{Q2CUT} = $q2cut  if ($q2cut);
print "Cutting Q2 at:: $ENV{Q2CUT}\n" if($debug);
$ENV{Q2BIN} =       8.;
$ENV{Q2BIN} = $q2bin if ($q2bin);
print "Variable Q2bin at:: $ENV{Q2BIN}\n" if($debug);
$ENV{PPLUSCUT} =     5.28;
$ENV{PPLUSCUT} = $ppluscut  if ($ppluscut);
print "Cutting P+ at:: $ENV{PPLUSCUT}\n" if($debug);
$ENV{PPLUSBIN} =     0.66;
$ENV{PPLUSBIN} = $pplusbin if ($pplusbin);
print "Variable P+bin at:: $ENV{PPLUSBIN}\n" if($debug);

#FULL SET ON DATA
$ENV{DOVARSTU}=        "1";
$ENV{NORE}=         1 if($nore);
$ENV{RE}=           1 if($re);
$ENV{DOBDECWEIGHT}=      0;
$ENV{DOBDECWEIGHT}=      1 if($b);
$ENV{DODDECWEIGHT}=      0;
$ENV{DODDECWEIGHT}=      1 if($d);
$ENV{DOFFWEIGHT}=      0;
$ENV{DOEXCLFFWEIGHT}=      0;
$ENV{DOCASCADEWEIGHT} = 0;
$ENV{DOCASCADEWEIGHT} = 1 if($cascade);
$ENV{DOFFWEIGHT}=      $ff if($ff);
$ENV{DOEXCLFFWEIGHT}=      $exclff if($exclff);
$ENV{DOSSBARWEIGHT} =  0;
$ENV{DOSSBARWEIGHT} =  $ssbar if($ssbar);
$ENV{DOMESMEANCORR} =  1 if($mesmeancorr);
$ENV{MIXCORR} =      2; #Apply Mixing Correction
$ENV{BTYPE}=         2;
if($chb eq "0" || $chb eq "1"){
    $ENV{BTYPE} = $chb;}
print "Running against $ENV{BTYPE} Bs\n" if($debug);
$ENV{FITTOTSHAPE}=   2;
if($cat) {$ENV{FITTOTSHAPE}=   1;}
$ENV{MULFAC}= "1";
$ENV{DELTAMB}= "0";
$ENV{DELTAA}= "0";
$ENV{FERMIAPP}= "0";
$ENV{EFFDFN} = "-99";
$ENV{EFFDFN} = $effdfn if($effdfn);
$ENV{SAVEPDFTREE} = 0;
$ENV{SAVEPDFTREE} = 1 if($savepdftree);
$ENV{READPDFTREE} = 0;
$ENV{READPDFTREE} = 1 if($readpdftree);
$ENV{TOYHISTOGRAMES} = 0;
$ENV{TOYHISTOGRAMES} = 1 if($toyhistogrames);


print "input::$ENV{MULFAC}; $ENV{DELTAMB}; $ENV{DELTAA}; $ENV{FERMIAPP}\n" if($debug);
if($mulf) {$ENV{MULFAC}= "$mulf";}


if($fermi) {
    print " do fermi $fermi\n";
  open(IN,"$fermi");
  while(<IN>) {
    chomp($_);
    if($_ =~ /deltamb (.*)/) {$ENV{DELTAMB}= "$1"};
    if($_ =~ /deltaa (.*)/) {$ENV{DELTAA}= "$1"};
    $ENV{FERMIAPP}= "1";
  }

}
print "input::$ENV{MULFAC}; $ENV{DELTAMB}; $ENV{DELTAA}; $ENV{FERMIAPP}\n" if($debug);

#$ENV{CHLOW}=      -3.5;
#$ENV{CHHIGH}=      3.5;
$ENV{CHLOW}=      -.5;
$ENV{CHHIGH}=     .5;
$ENV{EWPWLOCUT} = 3;
$ENV{EWPWHICUT} = 5.5;
$ENV{EWPWLOCUT} = $ewpwlo if($ewpwlo);
$ENV{EWPWHICUT} = $ewpwhi if($ewpwhi);
$ENV{Q2LOCUT} = 0;
$ENV{Q2HICUT} = 30;
$ENV{Q2LOCUT} = $qlo if($qlo);
$ENV{Q2HICUT} = $qhi if($qhi);
$ENV{CSILOCUT} = 0;
$ENV{CSIHICUT} = 1;
$ENV{CSILOCUT} = $csilo if($csilo);
$ENV{CSIHICUT} = $csihi if($csihi);
$ENV{XLOCUT} = 0;
$ENV{XHICUT} = 1.5;
$ENV{XLOCUT} = $xlo if($xlo);
$ENV{XHICUT} = $xhi if($xhi);
$ENV{WLOCUT} = 0;
$ENV{WHICUT} = 2;
$ENV{WLOCUT} = $wlo if($wlo);
$ENV{WHICUT} = $whi if($whi);
$ENV{USECB}= 1;
$ENV{USECB}= 0 if($fixsh);
$ENV{GAUSSFIT}=      0;
$ENV{BLIND}=         0;
$ENV{BLINDSIZE}=     0.4;
$ENV{FITMC}= 0;
$ENV{FITMC}= 1 if($fitmc);
$ENV{MNULOW}=     -100000.;
$ENV{EMPMLOW}=  -0.26;
$ENV{EMPMHIGH}=  0.26;
$ENV{PRMM2CUT}=    -3.;
$ENV{LEPTCUT}=       1.;
$ENV{TRULEPTCUT}=    1.;
$ENV{TRULEPTCUT}= $trulept if ($trulept);
$ENV{LEPTTYPE}=      2; #0 ==  elec ; 1 == muon ; 2 ==all
if($lptyp eq "0" || $lptyp eq "1"){
    $ENV{LEPTTYPE} = $lptyp;}
print "Running lepton type $ENV{LEPTTYPE} \n" if($debug);
$ENV{DOBRECOWEIGHT}=  0;
$ENV{DOPSTARWEIGHT}=  0;
$ENV{DOMM2WEIGHT}=  0;
$ENV{MININTPUR}=     0.; if ($minintpur > 0. && $minintpur < 1.) { $ENV{MININTPUR} = $minintpur; }
$ENV{MAXINTPUR}=     1000.;
$ENV{RUN}=           0;
$ENV{FITOPTION}=   0; #if ==3 fixes vcb fraction with fittoshape=0
if($opt){
  $ENV{FITOPTION}= $opt;
  $ENV{VCBCOMP}= 0;
  $ENV{OTHCOMP}= 0;
  if($opt == 5) {
    $ENV{VCBCOMP}= $vcb;
    $ENV{OTHCOMP}= $oth;
  }
}
$ENV{DEPL} =         0;
$ENV{DEPL} = 1 if ($depl);

$ENV{MNUHIGH}=   0.5;
$ENV{MNUHIGH}= $mnscan if($mnscan);
print "Cutting Mnu at:: $ENV{MNUHIGH}\n" if($debug);
if($comb) {
    $ENV{PSTARCOR}=      1.22599;
}else{
    $ENV{PSTARCOR}=      1.09242;
}

$ENV{REW} = 1;
$ENV{REW} = $rew if($rew);

$ENV{BINMES} = 0;
$ENV{BINMES} = 1 if($binmes);

$ENV{NOTUNBINMES} = 0;
$ENV{NOTUNBINMES} = $notunbinmes if($notunbinmes);

die  "Please Select mES fit model with flag -mesfitmodel <opt>.\n  <opt> = 0 for Argus and Gauss PDF\n  1 for Argus and Crystal Ball PDF 
  2 for Argus and Thorsten Signal AKA Frankenstein PDF\n  3 for three PDFs\n" unless $mesfitmodel;
$ENV{MESFITMODEL} = $mesfitmodel;

$ENV{FITDSS} = 0;
$ENV{FITDSS} = $ftdss if($ftdss);

$ENV{BINCUT} = 0.;
$ENV{BINCUT} = $bincut if($bincut);

$ENV{BRECOQUAL} = 0;
$ENV{BRECOQUAL} = 1 if($brecoqual);

$ENV{EMPM} = 0;
$ENV{EMPM} = 1 if($empm);

$ENV{FITALLMESRANGE} = 0;
$ENV{FITALLMESRANGE} = 1 if($fitallmesrange);

$ENV{SMALLSTATCORR} = 0;
$ENV{SMALLSTATCORR} = 1 if($smallstatcorr);

$ENV{SUBTRACTPEAKING} = 0;
$ENV{SUBTRACTPEAKING} = 1 if($subtractpkgbkg);
$ENV{SUBTRACTPEAKING} = 2 if($subtractpkgbkgmx);

$ENV{CASCADE} = 1;
$ENV{CASCADE} = $cascade if($cascade);

if(($subtractpkgbkg || $subtractpkgbkgmx) && !$countmc) {
    print "! ! ! -- W A R N I N G -- ! ! !  No PEAKING SUBTRACTION without MC Truth-matching counting (flag -countmc)\n";
    print " Setting SUBTRACTPEAKING to 0 ! \n";
    $ENV{SUBTRACTPEAKING} = 0;
}


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
    if($rel eq 14){ 
	$dirP = "/afs/slac.stanford.edu/u/br/sacco/newsel-23/workdir";}
    elsif($rel eq 18){ 
	#$dirP = "/afs/slac.stanford.edu/u/br/sacco/vol1/ana-32/workdir";}
	$dirP = "/nfs/farm/babar/AWG63/SemiLept/sacco/rootples/summer07/chains";}
    else{
#	$dirP = "/afs/slac.stanford.edu/u/br/gaglio/scra3/NewKs/chains";} # (Valid until may 2009. Old default)
	$dirP = "/nfs/farm/babar/AWG69/gaglio/SYS2/Standard/chains";} # New default with random numbers for systematics evaluation
#  	$dirP = "/nfs/farm/babar/AWG69/gaglio/SYS2/PidPiz/chains";} # sys
#       $dirP = "/nfs/farm/babar/AWG69/gaglio/SYS2/K0/Rif/chains"; } #default per KLong Sys
#	$dirP = "/nfs/farm/babar/AWG69/gaglio/SYS2/K0/KLP/chains"; } # per K Sys
    	
    if($small) {
	if($rel eq 14){ 
	    # default for summer05
	    $dirP="/nfs/farm/babar/AWG38/petrella/reduced/nkinfit-23";
	    # default for summer06
	    $dirP="/nfs/farm/babar/AWGsemilep01/menges/summer06/store.2/chains-1111";
	    
	    $fileTag = Red;
	    print "Running on reduced CM2 files\n" if($debug);}
	else{
	    #$dirP="/nfs/farm/babar/AWGsemilep01/petrella/Rootples/1111-R18/reduced"; #track killing prod
	    $dirP="/nfs/farm/babar/AWG63/SemiLept/petrella/summer07/reduced"; 
#	    $dirP="/afs/slac.stanford.edu/u/br/petrella/scra/"; 
	    $fileTag = Red;
	    print "Running on release 18 ntuples\n" if($debug);
	}
    }
    $dirP = "$dirS" if($dirS);

    if($rel eq 14){
	$runflag=14;
    }else{
	$runflag=15;
    }
    $runflag=4  if($Run4);
    $runflag=3  if($Run3);
    $runflag=5  if($Run5);
    $runflag=6  if($Run6);
    $runflag=12 if($Run12);
    $runflag=13 if($Run13);
    $runflag=14 if($Run14);
    $runflag=15 if($Run15);
    $runflag=16 if($Run16);
    
    my $extension="";
    if($small){ $extension = ".root"};
    
    my %generic; $generic{12} = "1-2"; $generic{3} = "3";   $generic{4} = "4"; $generic{14} = "1-4"; $generic{15} = "1-5"; $generic{5} = "5"; $generic{13} = "1-3"; $generic{6} = "6"; $generic{16} = "1-6";
    my %signal;  
    if($rel eq 14){
	$signal{12}  = "1-3"; $signal{3}  = "1-3"; $signal{4}  = "4"; $signal{14}  = "1-4";}
    else{
	$signal{12}  = ""; $signal{3}  = ""; $signal{4}  = ""; $signal{14}  = ""; $signal{15}=""; $signal{13}=""; $signal{16}=""; $signal{6}="";}

    print "Running on Run$runflag only!\n" if($debug);

    if($rel eq 14){
	$ENV{FILEVUBTOTALNRE}=  "${dirP}/nre${fileTag}ChainRun$signal{$runflag} 1";
	$ENV{FILEVUBTOTALHYB}=  "${dirP}/mix${fileTag}ChainRun$signal{$runflag} 1";
    }else{
	$ENV{FILEVUBTOTALNRE}=  "${dirP}/nre${fileTag}ChainRun$generic{$runflag}$extension 1";
	$ENV{FILEVUBTOTALHYB}=  "${dirP}/mix${fileTag}ChainRun$generic{$runflag}$extension 1";
    }
    $ENV{FILEDATA}=         "${dirP}/dat${fileTag}ChainRun$generic{$runflag}$extension 1";
    $ENV{FILEVCB}=          "${dirP}/gen${fileTag}ChainRun$generic{$runflag}$extension 1";
    if($rel eq 14){
	$ENV{FILEVCB1}=         "${dirP}/gen${fileTag}bchChain 1";
	$ENV{FILEVCB2}=         "${dirP}/gen${fileTag}bnuChain 1";
    }else{
	$ENV{FILEVCB1}=         "${dirP}/gen${fileTag}bchChainRun$generic{$runflag} 1";
	$ENV{FILEVCB2}=         "${dirP}/gen${fileTag}bnuChainRun$generic{$runflag} 1";
    }
    $ENV{FILEVUBTOTAL}=     "${dirP}/mix${fileTag}Chain 1";
    $ENV{FILEVUBTRUTHHYB}=  "/nfs/farm/babar/AWG37/ISL/kerstin/ntup/VubTruth_res.root";
    $ENV{FILEVUBTRUTHNRE}=  "/nfs/farm/babar/AWG37/ISL/kerstin/ntup/VubTruth_nonres.root";

    if($fitmc){

      $ENV{FILEVUBTOTALNRE}=  "${dirP}/nre${fileTag}ChainRun$generic{$runflag}$extension 1";
      $ENV{FILEVUBTOTALHYB}=  "${dirP}/mix${fileTag}ChainRun$generic{$runflag}$extension 1";
      $ENV{FILEDATA}=         "${dirP}/gen${fileTag}ChainRun$generic{$runflag}$extension 1";
      $ENV{FILEVCB}=          "${dirP}/gen${fileTag}ChainRun$generic{$runflag}$extension 1";
      $ENV{FILEVCB1}=         "${dirP}/gen${fileTag}bchChain 1";
      $ENV{FILEVCB2}=         "${dirP}/gen${fileTag}bnuChain 1";
      $ENV{FILEVUBTOTAL}=     "${dirP}/gen${fileTag}ChainRun$generic{$runflag}$extension 1";
      $ENV{FILEVUBTRUTHHYB}=  "/nfs/farm/babar/AWG37/ISL/kerstin/ntup/VubTruth_res.root";
      $ENV{FILEVUBTRUTHNRE}=  "/nfs/farm/babar/AWG37/ISL/kerstin/ntup/VubTruth_nonres.root";
      $ENV{DONTWEIGHTVUB}= 0;
    }
}

####################################################################################################

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
  print OUTSETUP " totalStat         $ENV{TOTSTAT} \n";
  print OUTSETUP " totalStatModel    $ENV{TOTSTATMOD} \n";
  print OUTSETUP " BRRatioGenValue   $ENV{BRRATIOGEN} \n";
  print OUTSETUP " BRRatioValueTailU $ENV{BRRATIOTAILU} \n";
  print OUTSETUP " BRRatioValueTailC $ENV{BRRATIOTAILC} \n";
  print OUTSETUP " pstarfact         $ENV{PSTARCOR} \n";
  print OUTSETUP " mxCut             $ENV{MXCUT} \n";
  print OUTSETUP " q2Bin             $ENV{Q2BIN} \n";
  print OUTSETUP " mxBin             $ENV{MXBIN} \n";
  print OUTSETUP " pplusBin          $ENV{PPLUSBIN} \n";
  print OUTSETUP " nore              $ENV{NORE} \n";
  print OUTSETUP " re                $ENV{RE} \n";
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
  print OUTSETUP " mnuSqLow          $ENV{MNULOW} \n";
  print OUTSETUP " mnuSqHigh         $ENV{MNUHIGH} \n";
  print OUTSETUP " EmpmLow           $ENV{EMPMLOW} \n";
  print OUTSETUP " EmpmHigh          $ENV{EMPMHIGH} \n";
  print OUTSETUP " chLow             $ENV{CHLOW} \n";
  print OUTSETUP " chHigh            $ENV{CHHIGH} \n";
  print OUTSETUP " depl              $ENV{DEPL} \n";
  print OUTSETUP " Btype             $ENV{BTYPE} \n";
  print OUTSETUP " q2Cut             $ENV{Q2CUT} \n";
  print OUTSETUP " ewpwloCut         $ENV{EWPWLOCUT} \n";
  print OUTSETUP " ewpwhiCut         $ENV{EWPWHICUT} \n";
  print OUTSETUP " q2loCut           $ENV{Q2LOCUT} \n";
  print OUTSETUP " q2hiCut           $ENV{Q2HICUT} \n";
  print OUTSETUP " csiloCut          $ENV{CSILOCUT} \n";
  print OUTSETUP " csihiCut          $ENV{CSIHICUT} \n";
  print OUTSETUP " xloCut            $ENV{XLOCUT} \n";
  print OUTSETUP " xhiCut            $ENV{XHICUT} \n";
  print OUTSETUP " wloCut            $ENV{WLOCUT} \n";
  print OUTSETUP " whiCut            $ENV{WHICUT} \n";
  print OUTSETUP " lepttype          $ENV{LEPTTYPE} \n";
  print OUTSETUP " fittotshape       $ENV{FITTOTSHAPE} \n";
  print OUTSETUP " mixcorr           $ENV{MIXCORR} \n"; 
  print OUTSETUP " fixcorrratio      $ENV{FIXCORRRATIO} \n"; 
  print OUTSETUP " fitMC             $ENV{FITMC} \n";
  print OUTSETUP " blinding          $ENV{BLIND} \n";
  print OUTSETUP " blindsize         $ENV{BLINDSIZE} \n";
  print OUTSETUP " randomseed        $ENV{RANDOMSEED} \n";
##  print OUTSETUP " seme              $ENV{SEME} \n";
  print OUTSETUP " service           $ENV{SERVICE} \n";
  print OUTSETUP " cascade           $ENV{CASCADE} \n";
  print OUTSETUP " issmearAll        $ENV{ISSMEARALL} \n";
  print OUTSETUP " smearAllMeanValue $ENV{SMEARALLMEANVALUE} \n";
  print OUTSETUP " smearAllSigma     $ENV{SMEARALLSIGMA} \n";
  print OUTSETUP " issmearBkg        $ENV{ISSMEARBKG} \n";
  print OUTSETUP " smearBkgMeanValue $ENV{SMEARBKGMEANVALUE} \n";
  print OUTSETUP " smearBkgSigma     $ENV{SMEARBKGSIGMA} \n";
  print OUTSETUP " dotrkreweight     $ENV{DOTRKWEIGHT} \n";
  print OUTSETUP " doneureweight     $ENV{DONEUWEIGHT} \n";
  print OUTSETUP " doBdecreweight    $ENV{DOBDECWEIGHT} \n";
  print OUTSETUP " doFFreweight      $ENV{DOFFWEIGHT} \n";
  print OUTSETUP " doexclFFreweight  $ENV{DOEXCLFFWEIGHT} \n";
  print OUTSETUP " doCascadeReweight $ENV{DOCASCADEWEIGHT} \n";
  print OUTSETUP " EffDFN            $ENV{EFFDFN} \n";
  print OUTSETUP " savepdftree       $ENV{SAVEPDFTREE} \n";
  print OUTSETUP " readpdftree       $ENV{READPDFTREE} \n";
  print OUTSETUP " toyhistogrames    $ENV{TOYHISTOGRAMES} \n";
  print OUTSETUP " doDdecreweight    $ENV{DODDECWEIGHT} \n";
  print OUTSETUP " dobrecoreweight   $ENV{DOBRECOWEIGHT} \n";
  print OUTSETUP " dossbarreweight   $ENV{DOSSBARWEIGHT} \n";
  print OUTSETUP " domesmeancorr     $ENV{DOMESMEANCORR} \n";
#  print OUTSETUP " dopstarreweight   $ENV{DOPSTARWEIGHT} \n";
#  print OUTSETUP " domm2reweight     $ENV{DOMM2WEIGHT} \n";
  print OUTSETUP " minintpur         $ENV{MININTPUR} \n";
  print OUTSETUP " maxintpur         $ENV{MAXINTPUR} \n";
  print OUTSETUP " run               $ENV{RUN} \n";
  print OUTSETUP " nnpi0             $ENV{CUTNNPI0} \n";
  print OUTSETUP " doleptplot        $ENV{LEPTPLOT} \n";
  print OUTSETUP " fermiapp          $ENV{FERMIAPP} \n";
  print OUTSETUP " deltamb           $ENV{DELTAMB} \n";
  print OUTSETUP " mulfac            $ENV{MULFAC} \n";
  print OUTSETUP " deltaa            $ENV{DELTAA} \n";
  print OUTSETUP " rew               $ENV{REW} \n";
  print OUTSETUP " binmes            $ENV{BINMES} \n";
  print OUTSETUP " notunbinmes       $ENV{NOTUNBINMES} \n";
  print OUTSETUP " mesfitmodel       $ENV{MESFITMODEL} \n";
  print OUTSETUP " fitdss            $ENV{FITDSS} \n";
  print OUTSETUP " bincut            $ENV{BINCUT} \n";
  print OUTSETUP " fitOption         $ENV{FITOPTION} \n";
  print OUTSETUP " vcbcomp           $ENV{VCBCOMP} \n";
  print OUTSETUP " othcomp           $ENV{OTHCOMP} \n";
  print OUTSETUP " dovarstu          $ENV{DOVARSTU} \n";
  print OUTSETUP " toyMC             0 \n";
  print OUTSETUP " dontweightvub     $ENV{DONTWEIGHTVUB} \n";
  print OUTSETUP " brecoqual         $ENV{BRECOQUAL} \n";
  print OUTSETUP " Empmcut           $ENV{EMPM} \n";
  print OUTSETUP " fitallmesrange    $ENV{FITALLMESRANGE}\n";
  print OUTSETUP " smallstatcorr     $ENV{SMALLSTATCORR} \n";
  print OUTSETUP " subtractpeaking   $ENV{SUBTRACTPEAKING} \n";
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
      $mpar = "mesparsetting_thorsten.dat" if ($mesfitmodel==2 || $mesfitmodel == 3);
  }
  $fitFlag .= " -Mes $mpar";
  

  if($thecomparison){ $fitFlag = $figFlag." -thecomparison"; }
  if($q2fit) {$fitFlag = $fitFlag." -q2 $q2fit";}
  if($comb) {$fitFlag = $fitFlag." -comb ";}
  if($Sun) {$fitFlag = $fitFlag." -Sun ";}
  if($ckBauer) {$fitFlag = $fitFlag." -ckBauer ";}
  if($re) {$fitFlag = $fitFlag." -re ";}
  if($nore) {$fitFlag = $fitFlag." -nore ";}
  if($unf) {
      $fitFlag = $fitFlag." -unf $unf";
      $fitFlag = $fitFlag." -U $unfF" if $unfF;
  } 
  if($mhiunf)  {$fitFlag = $fitFlag." -hiunf $mhiunf";}
  if($mx2unf) {$fitFlag = $fitFlag." -mx2unf ";}
  if($me)  {$fitFlag = $fitFlag." -me $me";}
  if($mu)  {$fitFlag = $fitFlag." -mu $mu";}
  if($mult)  {$fitFlag = $fitFlag." -mult";}
  if($Sys) {$fitFlag = $fitFlag." -Sys $Sys";}
  if($endpointCor) {$fitFlag .= " -epCor";}
  if($newbin) {$fitFlag = $fitFlag." -newbin ";}
  if($fixSBratio) {$fitFlag = $fitFlag." -sigpeakcorrmx $fixSBratio";}
  if($countmc) {$fitFlag = $fitFlag." -count";}
  if($dssRatio)  {$fitFlag = $fitFlag." -dssRatio $dssRatio";}
  if($dssFile) {$fitFlag = $fitFlag." -dssFile $dssFile";}
  if($bsys) {$fitFlag = $fitFlag." -bsys $bsys";}

  if($wisys)  {
    $fitFlag = $fitFlag." -wisys $wisys";
  } else {
    $fitFlag = $fitFlag." -wisys 1.";
  }

  if($q2vec){
    $fitFlag = $fitFlag." -q2V $q2vec";
  }else{
    $fitFlag = $fitFlag." -q2V spa/q2default.dat";
  }

  if($mx1d){
    $fitFlag = $fitFlag." -mx1d $mx1d";
  }else{
    if ($newbin) {
      $fitFlag = $fitFlag." -mx1d spa/mxmulti1d.dat";
    } else {
      $fitFlag = $fitFlag." -mx1d spa/mxdef1d.dat";
    }
  }
  if($pplus1d){
    $fitFlag = $fitFlag." -pplus1d $pplus1d";
  }else{
    if ($newbin) {
      $fitFlag = $fitFlag." -pplus1d spa/pplusmulti1d.dat";
    } else {
      $fitFlag = $fitFlag." -pplus1d spa/pplusdef1d.dat";
    }
  }

  if($mxvec){
    $fitFlag = $fitFlag." -mxV $mxvec";
  }else{
    $fitFlag = $fitFlag." -mxV spa/mxdefault.dat";
  }

  if($wF) {
    $fitFlag = $fitFlag." -W $wF";
  }else{
    $fitFlag = $fitFlag." -W sysWd/3dweights";
  }

  if($SPseed) { $fitFlag = $fitFlag." -SPseed $SPseed";}

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
      }else{
	  $fitFlag = $fitFlag. " -rel 18";
      }  
  }else{
      $fitFlag = $fitFlag. " -rel 14";
  }
 

  if($novarfit){ 
    $fitFlag = $fitFlag. " -varfit 0"; ## We are using q2, mxhad
  } else {
    $fitFlag = $fitFlag. " -varfit 1"; ## We are using q2fit, mxhadfit
  }

  $FFile = "wfermifile.dat" unless $FFile; # default
  $FFile = "" if $ENV{REW} == 12;          # not needed for this reweighting
  $fitFlag = $fitFlag. " -FFile $FFile" if $FFile;

###CB add possibility to save datasets (reduces initialization time)
  if($sdDir){
    $fitFlag = $fitFlag. " -sdDir $sdDir";
  }
  if($rdDir){
    $fitFlag = $fitFlag. " -rdDir $rdDir";
  }
  if($readpdftree) {
      $fitFlag = $fitFlag. " -readpdftree $readpdftree";
  }

  open OUT, ">>$tmpsett";
  print OUT "@ $sring\n";
  #####################

  my $queue = "xlong";
  $queue = $que if($que);
  if(!$pref) {$pref = "Dummy";}

  my $command = "../bin/$ENV{BFARCH}/VirFit -F $tmpfile -C $tmpsett -D $ENV{DIR} -P $ENV{TEST} $fitFlag";
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

  $ENV{DIR} = "~/scra/Ibu$ENV{TEST}/";
  $ENV{DIRSETTING} = "theset";
  $ENV{DIRFILES} = "thefiles";

  $ENV{NDATA} =     15000000;
  $ENV{NVUB} =      15000000;
  $ENV{NVCB} =      15000000;

#   $ENV{NDATA} =     150000; #just for test
#   $ENV{NVUB} =       40000; #just for test
#   $ENV{NVCB} =      150000; #just for test

  $ENV{TOTSTAT} =      270;
  $ENV{TOTSTATMOD} =   270;
  $ENV{BRRATIOGEN} =   0.0017;
  $ENV{BRRATIOTAILU} =  0.00224; #from published mx result on 80 fb^-1
  $ENV{BRRATIOTAILC} =  0.1061;  #from Incluisve BaBar measurements

  if($fitmc) {$ENV{BRRATIOTAILU} = 0.00215;}

  $ENV{FIXMEAN} =      1;
  $ENV{FIXSIGMA} =     1;
  $ENV{FIXARG1} =      0;
  $ENV{FIXARG} =       1;
  $ENV{FIXCB1} =       1;
  $ENV{FIXCB2} =       1;

  $ENV{FIXCORRRATIO} = 0;
  if($fixcorrratio){
      $ENV{FIXCORRRATIO} = 1;}

  $ENV{RANDOMSEED} =   990717;    # OFFICIAL BLINDING!!!!!
##  $ENV{SEME} = 0;
##  $ENV{SEME} = $SPseed if($SPseed);
  $ENV{SERVICE} = -1;
  $ENV{SERVICE} = $service if($service);
  $ENV{ISSMEARALL} =   0; #Default
  $ENV{SMEARALLMEANVALUE} = 0.1;
  $ENV{SMEARALLSIGMA} = 0.001;
  $ENV{ISSMEARBKG} =   0;
  $ENV{SMEARBKGMEANVALUE} = 0.0;  #Default
  $ENV{SMEARBKGSIGMA} =     0.1;
  $ENV{DOTRKWEIGHT} =      0;
  $ENV{DONEUWEIGHT} =      0;
  #$ENV{INTPUR} =       0.;
  #$ENV{DOREWEIGHT} = 0;
  $ENV{LEPTPLOT} = 0;
  $ENV{CUTNNPI0} =    -1000;
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
  -rel i   :       Selects release i=14 or i=18 ntuples for CM2
  -dirS=s  :       Take chains from directory s
  -sdDir=s :       Save RooDataSets to directory s
  -rdDir=s :       Take RooDataSets to directory s
  -small   :       Run on reduced files (CM2 only; take chains from chain/ subdir)
  -Run12   :       Run only on Run1-2 data (CM2 only; take chains from chain/ subdir)
  -Run3    :       Run only on Run3 data (CM2 only; take chains from chain/ subdir)
  -Run4    :       Run only on Run4 data (CM2 only; take chains from chain/ subdir)

  -mesfitmodel  i:  Select a model to fit mES distributions: 
                       i=0 (Argus+Gaussian); 
		       i=1 (Argus+Crystal Ball); 
                       i=2 (Argus+Thorsten Signal PDF aka Frankenstein);
                       i=3 (3 PDFs);

  -distfit :       1D fit with mx (0), q2 (1) or pplus (2) distribution
  -opt 1   :       do a 3-parameter (vub,vcb,other) 1D fit. Default is 2-parameter. 
  -opt 5   :       xxx
  -cat     :       Fit separately vcb and oth comp
  -comb    :       2d fit with [mx,q2]
  -mx1d    :       choose mx favourite binning (1d only)
  -pplus1d :       choose pplus favourite binning (1d only)
  -q2V     :       choose q2 favourite binning (2d only)
  -mxV     :       choose mx favourite binning (2d only)

  -b       :       B reweighting
  -d       :       D reweighting
  -ssbar i :       ssbar reweighting (non-res vub: 1=-30% , 2=+30%)
  -me=i    :       systematics for floating parameters in PDFs
  -mxscan  :       Mx upper cut
  -mxbin   :       Mx bin scan
  -q2cut   :       Q2 upper cut
  -q2bin   :       Q2 bin scan
  -ppluscut:       P+ upper cut
  -pplusbin:       P+ bin scan
  -mnscan  :       upper cut on missing mass squared (for scans)
  -nore    :       Select nonresonant only 
  -re      :       Select    resonant only
  -novarfit:       use measured variables instead of fitted ones
  -minintpur:      minimum integrated purity (default: 0.)
  -(var)lo :       Low var cut (var can be: csi,x,w,ewpw)
  -(var)hi :       High var cut (var can be: csi,x,w,ewpw)
  -rew     :       Apply vub reweighting 1: 1D, 3:3D, 10 a la dominique, 11 a la dominique with latest mx binning
  -wF      :       Applies Theo model reweighting
  -FFile   :       Use indicated file instead of wfermifile.dat
  -Sun     :       Apply signal unfolding to get efficiencies from Bauer et al
  -ckBauer :       Apply signal unfolding to get efficiencies from MC

  -unf     :       Equidistant mX binning (-unf numbins)
  -mhiunf  :       Endpoint for equidistant binning
  -mx2unf  :       Prepare the spectra for mx2 unfolding

  -mpar=s  :       File with mES starting parameters
  -tmodel  :       use the new pdfs from Thorsten
  -binmes  :       Do binned mES fits (a-la-prl) for old pdfs
  -notunbinmes:    Do binned mES fits
  -fixshape:       Float the shape of the Argus function in mES fits 
  -depl    :       Run on depleted sample
  -thecomparison : Dump dataset to be used with the comparison
  -fitdss  :       (mx 1-D fit only) =1: put D** together with other background
                                     =2: D** is the "other" background (i.e. SL events in vcb)
				   WARNING: fitdss screws up the calculation of Nsl and other things
                                            use only to evaluate the effect on the fit yields 
  -fixSBratio=s :  Use file s to fix the signal/peaking background ratio in mES data fits 
  -countmc :       Take number of truth-matched MC events above mES>5.27 instead of mES fits 
                   when computing efficiencies, corrections, shapes of kinematic variables, etc. 
  -dssRatio=f:     Force the ratio of (D**lnu+nonSL)/(Dlnu+D*lnu) to f in the background (=vcboth) 
                   component when doing signal unfolding fits (option -Sun) on kinematic variables
  -dssFile=s:      Reweight the D** mX spectrum bin-by-bin as specified in file s (example in dssrew.dat)
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

