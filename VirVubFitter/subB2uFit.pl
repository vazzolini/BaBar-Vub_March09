#!/usr/local/bin/perl
use Getopt::Long;
&GetOptions('flag=s',\$flag,'Sys=i',\$Sys,'int',\$int,'q2fit',\$q2fit,'comb',\$comb,'FF',\$ff,'d',\$d,'cat',\$cat,'b',\$b,'que=s',\$que,'debug',\$debug,'chb=i',\$chb,'opt=i',\$opt,'qcut=s',\$qcut,'ewpwlo=s',\$ewpwlo,'ewpwhi=s',\$ewpwhi,'wF=s',\$wF,'csilo=s',\$csilo,'csihi=s',\$csihi,'xlo=s',\$xlo,'xhi=s',\$xhi,'dirO=s',\$dirO,'dirS=s',\$dirS,'wlo=s',\$wlo,'whi=s',\$whi,'qlo=s',\$qlo,'qhi=s',\$qhi,'pref=s',\$pref,'fermi=s',\$fermi,'mxscan=s',\$mxscan,'mnscan=s',\$mnscan,'q2bin=s',\$q2bin,'mxbin=s',\$mxbin,'mulf=s',\$mulf,'vcb=s',\$vcb,'oth=s',\$oth,'help',\$help,'unf=s',\$unf,'mhiunf=f',\$mhiunf,'mult',\$mult,'Sun',\$Sun,'ckBauer',\$ckBauer,'nore',\$nore,'re',\$re,'me=i',\$me,'mu=i',\$mu,'rew=i',\$rew,'fixshape',\$fixsh,'binmes',\$bnms,'q2V=s',\$q2vec,'mxV=s',\$mxvec,'mx1d=s',\$mx1d,'wisys=s',\$wisys,'merged',\$merged,'mxfit',\$mxfit,'newbin',\$newbin,'width=s',\$width,'b2uw=s',\$b2uw);
#'mpar=s',\$mpar

usage(0) if($help);

#FULL SET ON DATA
$ENV{DOVARSTU}=        "1";
$ENV{MXCUT}=        5.0;
$ENV{MXCUT}= $mxscan if($mxscan);
print "Cutting Mx at:: $ENV{MXCUT}\n" if($debug);
$ENV{Q2BIN}=        12.;
$ENV{Q2BIN}= $q2bin if($q2bin);
print "Variable Q2bin at:: $ENV{Q2BIN}\n" if($debug);
$ENV{MXBIN}=        1.67;
$ENV{MXBIN}= $mxbin if($mxbin);
print "Variable Mxbin at:: $ENV{MXBIN}\n" if($debug);
$ENV{NORE}=         1 if($nore);
$ENV{RE}=           1 if($re);
$ENV{DOBDECWEIGHT}=      0;
$ENV{DOBDECWEIGHT}=      1 if($b);
$ENV{DODDECWEIGHT}=      0;
$ENV{DODDECWEIGHT}=      1 if($d);
$ENV{DOFFWEIGHT}=      0;
$ENV{DOFFWEIGHT}=      1 if($ff);
$ENV{BTYPE}=         2;
$ENV{BTYPE}= $chb if($chb);
print "Running against $ENV{BTYPE} Bs\n" if($debug);

$ENV{FITTOTSHAPE}=   2;
if($cat) {$ENV{FITTOTSHAPE}=   1;}
$ENV{MULFAC}= "1";
$ENV{DELTAMB}= "0";
$ENV{DELTAA}= "0";
$ENV{FERMIAPP}= "0";
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
$ENV{CHHIGH}=      .5;
#Added q2 cut
$ENV{Q2CUT} =        0;
$ENV{Q2CUT} =        $qcut  if($qcut);
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
$ENV{FITMC}=         "0";
#######$ENV{FITMC}=         "1";
$ENV{MNULOW}=     -1.;

$ENV{PRMM1CUT}=      -20.;
$ENV{PRMM2CUT}=      -3.;
$ENV{PRMM3CUT}=      -20.;
$ENV{COSTMISSCUTLO}=   -0.95;
$ENV{COSTMISSCUTHI}=   0.95;
$ENV{PCMSTRKLOCUT}=    0.120;
$ENV{PMISSCUTLO}=      0.300;

$ENV{LEPTCUT}=       1.;
$ENV{LEPTTYPE}=      "2"; #0 ==  elec ; 1 == muon ; 2 ==all
$ENV{DOBRECOWEIGHT}=  0;
$ENV{DOPSTARWEIGHT}=  0;
$ENV{DOMM2WEIGHT}=  0;
$ENV{MININTPUR}=     0.;
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

#ENRICHED
#my $dirP = "/u/ec/asarti/Unfol13b/VubAnalysis/ed_fil/";
#$dirP = "$dirS" if($dirS);
###$ENV{DEPL}=         1;
#$ENV{DEPL}=         0; # default
#$ENV{FILEVUBTOTALNRE}= "${dirP}csx-vubnre.root";
##$ENV{FILEVUBTOTALNRE}= "${dirP}csx-cocktail-new.root";
##$ENV{FILEVUBTOTALNRE}=  "/afs/slac.stanford.edu/g/babar/work/k/kerstin/myana1252/fitrt/csx-allBu-mult2-vubnre.root";
#$ENV{FILEVUBTOTALHYB}= "${dirP}csx-vubmix.root";
##$ENV{FILEVUBTOTALHYB}= "${dirP}csx-cocktail-new.root";
##$ENV{FILEVUBTOTALHYB}= "/afs/slac.stanford.edu/g/babar/work/k/kerstin/myana1252/fitrt/csx-allBu-mult-vubmix.root";
#$ENV{FILEDATA}=         "/nfs/farm/babar/AWG6/Recoil/newprod/root/anaQA-prl/csx-data.root";
##$ENV{FILEDATA}=         "${dirP}csx-cocktail-new.root";
##$ENV{FILEDATA}=         "${dirP}csx-data.root";
##$ENV{FILEDATA}=         "/u/ec/asarti/Unfol13b/VubAnalysis/datasec/csx-gen-New.root";
##$ENV{FILEDATA}=         "${dirP}csx-vubnre.root";
#$ENV{FILEVCB1}=         "${dirP}csx-genbch-new.root";
#$ENV{FILEVCB2}=         "${dirP}csx-genbnu-new.root";
#$ENV{FILEVCB}=          "/u/ec/asarti/Unfol13b/VubAnalysis/datasec/csx-gen-New.root";
##$ENV{FILEVCB}=          "${dirP}csx-cocktail-new.root";
##$ENV{FILEVCB}=          "${dirP}csx-vubmix.root";
#$ENV{FILEVUBTOTAL}=     "${dirP}csx-vubmix.root";
##$ENV{FILEVUBTOTAL}=     "${dirP}csx-cocktail-new.root";
##$ENV{FILEVUBTOTAL}=     "/afs/slac.stanford.edu/g/babar/work/k/kerstin/myana1252/fitrt/csx-allBu-mult2-vubmix.root";
############################################################################################
####################################  PRLfiles  ############################################
#my $dirP = "/u/ec/ursl/root/anaQA-prl/";
##$dirP = "$dirS" if($dirS);
#$ENV{DEPL}=         0; # default
#$ENV{FILEVUBTOTALNRE}= "${dirP}csx-vubnre.root";
#$ENV{FILEVUBTOTALHYB}= "${dirP}csx-vubmix.root";
#$ENV{FILEDATA}=         "${dirP}csx-data.root";
#####$ENV{FILEVCB}=          "${dirP}csx-genb-new.root";
####$ENV{FILEVCB1}=         "${dirP}csx-genbch-new.root";
#####$ENV{FILEVCB2}=         "${dirP}csx-genbnu-new.root";
#$ENV{FILEVUBTOTAL}=     "${dirP}csx-vubmix.root";
############################################################################################
####################################  files run by RIC  ############################################
# my $dirP = "/nfs/farm/babar/AWG6/Recoil/newprod/root/anaQA-prl/";
# ##$dirP = "$dirS" if($dirS);
# $ENV{DEPL}=         0; # default
# $ENV{FILEVUBTOTALNRE}= "${dirP}csx-vubnre.root";
# $ENV{FILEVUBTOTALHYB}= "${dirP}csx-vubmix.root";
# $ENV{FILEDATA}=         "${dirP}csx-data.root";
# $ENV{FILEVCB}=          "${dirP}csx-genb-new.root";
# $ENV{FILEVCB1}=         "${dirP}csx-genbch-new.root";
# $ENV{FILEVCB2}=         "${dirP}csx-genbnu-new.root";
# $ENV{FILEVUBTOTAL}=     "${dirP}csx-vubmix.root";
############################################################################################
########### small faster root files, already selected semileptonic events ##################
#my $dirP = "/nfs/farm/babar/AWG6/Recoil/newprod/root/anaQA-prl/reduced/";
#$dirP = "$dirS" if($dirS);
##$ENV{DEPL}=        1;
#$ENV{DEPL}=         0; #default
#$ENV{FILEVUBTOTALNRE}=  "${dirP}csx-vubnre_small.root";
#$ENV{FILEVUBTOTALHYB}=  "${dirP}csx-vubmix_small.root";
#$ENV{FILEDATA}=         "${dirP}csx-data_small.root";
#$ENV{FILEVCB}=          "${dirP}csx-genb-new_small.root";
#$ENV{FILEVCB1}=         "${dirP}csx-genbch-new_small.root";
#$ENV{FILEVCB2}=         "${dirP}csx-genbnu-new_small.root";
#$ENV{FILEVUBTOTAL}=     "${dirP}csx-vubmix_small.root";
############################################################################################
############ chains for systematic studies                                 ##################
# my $dirP = "/nfs/babar/recoil/Vub_incl/root/anaQA-r00/";
# $dirP = "$dirS" if($dirS);
# print "====> $dirP\n"; 
 
# #$ENV{DEPL}=         1;
# $ENV{DEPL}=         0;
# $ENV{FILEVUBTOTALNRE}=  "${dirP}nreChain 1";
# $ENV{FILEVUBTOTALHYB}=  "${dirP}mixChain 1";
# #####$ENV{FILEVUBTOTALNRE}=  "${dirP}nreChain_red 1";
# #####$ENV{FILEVUBTOTALHYB}=  "${dirP}mixChain_red 1";
# $ENV{FILEDATA}=         "/nfs/babar/recoil/Vub_incl/root/anaQA-r00/datChain 1";
# #$ENV{FILEVCB}=          "/nfs/babar/recoil/Vub_incl/root/anaQA-r06/genChain 1 "; ##ff
# $ENV{FILEVCB}=          "${dirP}genChain 1";
# $ENV{FILEVCB1}=         "${dirP}genbchChain 1";
# $ENV{FILEVCB2}=         "${dirP}genbnuChain 1";
# $ENV{FILEVUBTOTAL}=     "${dirP}mixChain 1";
#############################################################################################


if ($merged) {
  my $dirP = "/nfs/babar/recoil/Vub_incl/root/anaQA-r00/";
  $dirP = "$dirS" if($dirS);
  ###$ENV{DEPL}=          1;
  $ENV{DEPL}=             0; #default
  $ENV{FILEVUBTOTALNRE}=  "${dirP}csx-vubnre.root";
  $ENV{FILEVUBTOTALHYB}=  "${dirP}csx-vubmix.root";
  $ENV{FILEDATA}=         "${dirP}csx-data.root";
  $ENV{FILEVCB}=          "${dirP}csx-genb-new.root";
  $ENV{FILEVCB1}=         "${dirP}csx-genbch-new.root";
  $ENV{FILEVCB2}=         "${dirP}csx-genbnu-new.root";
  $ENV{FILEVUBTOTAL}=     "${dirP}csx-vubmix.root";
} else {
  my $dirP = "/nfs/babar/recoil/Vub_incl/root/anaQA-r00/";
  $dirP = "$dirS" if($dirS);
  ###$ENV{DEPL}=        1;
  $ENV{DEPL}=         0; #default
  $ENV{FILEVUBTOTALNRE}=  "${dirP}nreChain 1";
  $ENV{FILEVUBTOTALHYB}=  "${dirP}mixChain 1";
  $ENV{FILEDATA}=         "${dirP}datChain 1";
  $ENV{FILEVCB}=          "${dirP}genChain 1";
  $ENV{FILEVCB1}=         "${dirP}genbchChain 1";
  $ENV{FILEVCB2}=         "${dirP}genbnuChain 1";
  $ENV{FILEVUBTOTAL}=     "${dirP}mixChain 1";
}

$ENV{MNUHIGH}=       .5;
$ENV{MNUHIGH}= $mnscan if($mnscan);
print "Cutting Mnu at:: $ENV{MNUHIGH}\n" if($debug);
$ENV{PSTARCOR}=      1.13152;

$ENV{REW} = 1;
$ENV{REW} = $rew if($rew);


$ENV{BINMES} = 0;
$ENV{BINMES} = 1 if($bnms);

$ENV{TEST}=   "MainTest";
if($flag) {$ENV{TEST}= "$flag"};
&scriptfit();


sub scriptfit {
  
  &envfit();
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
  print OUTSETUP " BRRatioValueTail  $ENV{BRRATIOTAIL} \n";
  print OUTSETUP " pstarfact         $ENV{PSTARCOR} \n";
  print OUTSETUP " mxCut             $ENV{MXCUT} \n";
  print OUTSETUP " q2Bin             $ENV{Q2BIN} \n";
  print OUTSETUP " mxBin             $ENV{MXBIN} \n";
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
  print OUTSETUP " prmm1Cut          $ENV{PRMM1CUT} \n";
  print OUTSETUP " prmm2Cut          $ENV{PRMM2CUT} \n";
  print OUTSETUP " prmm3Cut          $ENV{PRMM3CUT} \n";
  print OUTSETUP " pmissCutlo        $ENV{PMISSCUTLO} \n";
  print OUTSETUP " costmissCutlo     $ENV{COSTMISSCUTLO} \n";
  print OUTSETUP " costmissCuthi     $ENV{COSTMISSCUTHI} \n";
  print OUTSETUP " pcmsTrkloCut      $ENV{PCMSTRKLOCUT} \n";
  print OUTSETUP " mnuSqLow          $ENV{MNULOW} \n";
  print OUTSETUP " mnuSqHigh         $ENV{MNUHIGH} \n";
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
  print OUTSETUP " fitMC             $ENV{FITMC} \n";
  print OUTSETUP " blinding          $ENV{BLIND} \n";
  print OUTSETUP " blindsize         $ENV{BLINDSIZE} \n";
  print OUTSETUP " randomseed        $ENV{RANDOMSEED} \n";
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
  print OUTSETUP " doDdecreweight    $ENV{DODDECWEIGHT} \n";
  print OUTSETUP " dobrecoreweight   $ENV{DOBRECOWEIGHT} \n";
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
  print OUTSETUP " fitOption         $ENV{FITOPTION} \n";
  print OUTSETUP " vcbcomp           $ENV{VCBCOMP} \n";
  print OUTSETUP " othcomp           $ENV{OTHCOMP} \n";
  print OUTSETUP " dovarstu          $ENV{DOVARSTU} \n";
  print OUTSETUP " toyMC             0 \n";
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
  print OUTFILE "  \n";
  close OUTFILE;

  system("rm -rf $ENV{DIR}");
  system("mkdir $ENV{DIR}");


 # if($mpar){
#    $fitFlag = $fitFlag."-Mes $mpar";
#  } else{
#    $fitFlag = $fitFlag."-Mes mesparsetting.dat";
#  }
  my $fitFlag = "-Mes mesparsetting.dat";
  # my $fitFlag = "-Mes mespar_ric.dat";   

  #my $fitFlag = " ";


  if($q2fit) {$fitFlag = $fitFlag." -q2 ";}
  if($comb) {$fitFlag = $fitFlag." -comb ";}
  if($Sun) {$fitFlag = $fitFlag." -Sun ";}
  if($mxfit) {$fitFlag = $fitFlag." -mxfit ";}
  if($newbin) {$fitFlag = $fitFlag." -newbin ";}
  if($ckBauer) {$fitFlag = $fitFlag." -ckBauer ";}
  if($re) {$fitFlag = $fitFlag." -re ";}
  if($nore) {$fitFlag = $fitFlag." -nore ";}
  if($unf) {$fitFlag = $fitFlag." -unf $unf"; $fitFlag = $fitFlag." -U Unfpar.dat";} 
  if($mhiunf)  {$fitFlag = $fitFlag." -hiunf $mhiunf";}
  if($me)  {$fitFlag = $fitFlag." -me $me";}
  if($mu)  {$fitFlag = $fitFlag." -mu $mu";}
  if($mult)  {$fitFlag = $fitFlag." -mult";}
  if($Sys) {$fitFlag = $fitFlag." -Sys $Sys";}
  if($me)  {$fitFlag = $fitFlag." -me $me";}
  if($mu)  {$fitFlag = $fitFlag." -mu $mu";}

  if($wisys)  {
    $fitFlag = $fitFlag." -wisys $wisys";
  } else {
    $fitFlag = $fitFlag." -wisys 1.";
  }

  if($width)  {
    $fitFlag = $fitFlag." -width $width";
  } else {
    $fitFlag = $fitFlag." -width 0.3";
  }

  if($q2vec){
    $fitFlag = $fitFlag." -q2V $q2vec";
  }else{
    $fitFlag = $fitFlag." -q2V spa/q2default.dat";
  }

  if($mxvec){
    $fitFlag = $fitFlag." -mxV $mxvec";
  }else{
    $fitFlag = $fitFlag." -mxV spa/mxdefault.dat";
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

  if($b2uw) {
    $fitFlag = $fitFlag." -b2uw $b2uw";
  }else{
    $fitFlag = $fitFlag." -b2uw b2uweights";
  }

  if($wF) {
    $fitFlag = $fitFlag." -W $wF";
  }else{
    $fitFlag = $fitFlag." -W sysWd/3dweights";
  }
  my $queue = "long";
  $queue = $que if($que);
  if(!$pref) {$pref = "Dummy";}
  if (!$int) {
    if($debug) {print "bsub -G babarISemil -q $queue -o $ENV{DIR}/vubfit.out ../bin/$ENV{BFARCH}/b2uFit -F $tmpfile -C $tmpsett -D $ENV{DIR} -P $pref  $fitFlag\n"};
    system("bsub -G babarISemil -q $queue -o $ENV{DIR}/vubfit.out ../bin/$ENV{BFARCH}/b2uFit -F $tmpfile -C $tmpsett -D $ENV{DIR} -P $ENV{TEST}  $fitFlag");
#    if($debug) {print "bsub -q $queue -o $ENV{DIR}/vubfit.out ../bin/$ENV{BFARCH}/b2uFit -F $tmpfile -C $tmpsett -D $ENV{DIR} -P $pref  $fitFlag\n"};
#    system("bsub -q $queue -o $ENV{DIR}/vubfit.out ../bin/$ENV{BFARCH}/b2uFit -F $tmpfile -C $tmpsett -D $ENV{DIR} -P $ENV{TEST}  $fitFlag");
  } else {
    if($debug) {print "../bin/$ENV{BFARCH}/b2uFit -F $tmpfile -C $tmpsett -D $ENV{DIR} -P $ENV{TEST}  $fitFlag"};
    system("../bin/$ENV{BFARCH}/b2uFit -F $tmpfile -C $tmpsett -D $ENV{DIR} -P $ENV{TEST}  $fitFlag");
  }
  #-Fermi
}

sub envfit {
  #$ENV{TEST data

  $ENV{DIR} = "please/setme/tosomething";
  $ENV{DIR} = "$dirO/b2u-$ENV{TEST}/" if($dirO);
  $ENV{DIRSETTING} = "theset";
  $ENV{DIRFILES} = "thefiles";

  $ENV{NDATA} =     15000000;
  $ENV{NVUB} =       4000000;
  $ENV{NVCB} =      15000000;

  #########################TEMP
  #$ENV{NDATA} =     150000; 
  #$ENV{NVUB} =      150000;
  #$ENV{NVCB} =      150000;
  

  #$ENV{NVUBDATA} =      6400;    ##MC
  #$ENV{NVCBDATA} =    213500;    ##MC
  ############################
 # $ENV{NVUBDATA} =     35857;    ##MC ricomputed for chains
  #$ENV{NVCBDATA} =   1630000;    ##MC ricomputed for chains

  ############################
  # $ENV{NVUBDATA} =      1950;    ## ricompute for reduced file 
  # $ENV{NVCBDATA} =     65000;    ## ricompute for reduced file


  $ENV{TOTSTAT} =      270;
  $ENV{TOTSTATMOD} =   270;
  $ENV{BRRATIOGEN} =   0.0017;
  $ENV{BRRATIOTAIL} =  0.0017;

  $ENV{FIXMEAN} =      1;
  $ENV{FIXSIGMA} =     1;
  $ENV{FIXARG1} =      0;
  $ENV{FIXARG} =       1;
  $ENV{FIXCB1} =       1;
  $ENV{FIXCB2} =       1;

  $ENV{MIXCORR} =      2;
  $ENV{RANDOMSEED} =   990717;    # OFFICIAL BLINDING!!!!!
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
  -b       :       B rew
  -d       :       D rew
  -q2V     :       choose q2 favourite binning (2d only)
  -mxV     :       choose mx favourite binning (2d only)
  -mx1d    :       choose mx favourite binning (1d only)
  -qcut    :       Q2 cut
  -q2fit   :       Fit with Q2 distribution
  -q2bin   :       Bin scan
  -nore    :       Select nonresonant only 
  -re      :       Select    resonant only
  -(var)lo :       Low var cut (var can be: csi,x,w,ewpw)
  -(var)hi :       High var cut (var can be: csi,x,w,ewpw)
  -cat     :       Fti separately vcb and oth comp
  -flag    :       Flagging output
  -int     :       job is run interactively
  -debug   :       Additional output/printout
  -help    :       print this message.
  -wF      :       Applies Theo model reweighting
  -comb    :       Combined cut analysis [mx,q2]
  -Sun     :       Apply signal unfolding to get efficiencies from Bauer et al
  -ckBauer :       Apply signal unfolding to get efficiencies from MC
  -unf     :       Equidistant mX binning (-unf numbins)
  -mhiunf  :       Endpoint for equidistant binning
  -rew     :       Apply vub reweighting 1: 1D, 3:3D, 10 a' la dominique
  -fixshape:       Fix the shape of the Argus function in mES fits 
  -binmes  :       Do binned mES fits (a-la-prl)
  -opt 1   :       do a 3-parameter (vub,vcb,other) mX fit. Default is 2-parameter. 

Examples:
$0 -int -d -b -wF sysWd/3dweights -debug -flag ciappi

The following reproduces the prl result as closely as possible
(binned mES fits, 3-parameters mX fit): 
$0 -int -d -b -binmes -rew 1 -wf sysWd/ -opt 1 -debug -flag prlresult

The following will fit prl-like, but with unbinned mES fits:
$0 -int -d -b -rew 1 -wf sysWd/ -opt 1 -debug -flag unbin_prlresult

The following will fit prl-like, unbinned mES fits, 2-parameter mX fit: 
$0 -int -d -b -rew 1 -wf sysWd/ -debug -flag unbin_prlresult_2par

Comments to <azzolini\@slac.stanford.edu>.
INLINE_LITERAL_TEXT
#'
  exit($exit) if defined $exit;
}

