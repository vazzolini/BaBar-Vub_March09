#!/usr/bin/perl

use Getopt::Long;

die unless &GetOptions( 'data', \$data, 'MC', \$MC, 'help',\$help, 'bkg', \$bkg, 'all',\$all,
			'sig',\$sig, 'oldmodel', \$oldmodel, 'newmodel', \$newmodel, 'allcut', \$allcut, 'lepcut', \$lepcut, 
			'enr', \$enr, 'depl',\$depl, 'dumppar', \$dumppar, 'parfile=s',\$parfile,'mxl=f',\$mxl,'mxh=f',\$mxh,
			'q2l=f',\$q2low, 'q2h=s',\$q2h,'scanmx',\$scanmx,'scanmxq2',\$scanmxq2,'scanpplus',\$scanpplus, 
			'intpur=f', \$intpur, 'chb=i', \$chb,'haspi0=i', \$haspi0, 'ppl=f', \$ppl, 'pph=f', \$pph);

usage(0) if($help);

################ SOME DEFAULT ###########
my $mxsize=1;
my $q2size=1;
my $pplussize=1;
my @mxcut=(-300,300);
my @q2cut=(-999,300);
my @ppluscut=(-100,100);
my $pfile;
my $intp=0;
my $ch=2;
my $haspiz=-1;

############### DATA/ MC  ###########
my $isdata;
$isdata=1 if($data);
$isdata=0 if($MC);
if(!($data  ||  $MC)){ die "\nPlease specify DATA or MC!   See help: ./submit.pl -help\n";}
############## SELECT MODEL (2 or 3 pdfs) ##############
my $model;
$model=0 if($oldmodel);
$model=1 if($newmodel);
if(!($oldmodel || $newmodel)){$model=1; print "\n USING NEW MODEL (Frankestein + Argus)\n";}

############## FIT SIGNAL or BKG, or ALL  ###########
my $pdfcomp;

$pdfcomp=1 if($sig);
$pdfcomp=2 if($bkg);
$pdfcomp=3 if($all);

if (!($all || $sig || $bkg)) {die "\nPlease specity the pdf to be used in the fit!     See help: ./submit.pl -help\n";}
############# ALL CUT/ LEPCUT ###########
my $isallcut;
$isallcut=0 if($lepcut);
$isallcut=1 if($allcut);
if (! ($lepcut || $allcut)) {die "\nPlease specity allcut/lepcut    See help: ./submit.pl -help\n";}
############ DEPLETED/ENRICHED ###########
my $isdepleted;
$isdepleted=0 if($enr);
$isdepleted=1 if($depl);
if (! ($enr || $depl)) {die "\nPlease specity enriched/depleted sample    See help: ./submit.pl -help\n";}
############  DUMP PARAMETERS ###########
my $dumpp=0;
$dumpp=1 if($dumppar);
############  FILE WITH PARAMETERS ###########
if (! ($parfile)) {
    if($data){ 
	print "\nUsing standard parameters file parameters_DATA.txt\n";
	$pfile="parameters_DATA.txt";}
    if($MC) {
	print "\nUsing standard parameters file parameters_MC.txt\n";
	$pfile="parameters_MC.txt";}
}
$pfile=$parfile if($parfile);
########### MX cut ###########

if($scanmx){
   @mxcut=(0,1.55,1.9,2.2,2.5,2.8,3.1,3.4,3.7,4.2,5);
   $mxsize=10;
#   @q2cut=(-999.,300.); forse inutile
#   $q2size=1; forse inutile
}
if($scanmxq2){
    @mxcut=(0,1.77,2.2,2.8,5);
    @q2cut=(0,2,4,6,8,10,12,14,26);
    $mxsize=4;
    $q2size=8;
}
if($mxl || $mxh){
    @mxcut=($mxl,$mxh);
}
if($q2low || $q2high){
    @q2cut=($q2low,$q2high);
}
if($ppl || $pph){
    @ppluscut=($ppl,$pph);
}
if($scanpplus){
    @ppluscut=(0,0.66,1.32,1.98,2.64,3.3,3.96,4.62,5.28);
    $pplussize=8;
   # @q2cut=(-999.,300.); forse inutile
   # $q2size=1; forse inutile
}
    

########### INTEGR. PURITY ###########
$intp=$intpur if($intpur);

########### B Charge ###########
if($chb eq "0" || $chb eq "1") {$ch=$chb}

########### Modes with pi0 in BRECO #########
if($haspi0 eq "0" || $haspi0 eq "1") {$haspiz=$haspi0}

for($ppi=0; $ppi<$pplussize; $ppi++){
    for($q2i=0; $q2i<$q2size; $q2i++){
	for($mxi=0;$mxi<$mxsize; $mxi++){
	    if( -e "tempjob.tmp"){ system("rm tempjob.tmp");}
	    if( -e "autojob.C"){ system("rm autojob.C");}
	    open JOBFILE, ">tempjob.csh";
	    open ROOTMACROFILE, ">autojob.C";
	    
	    print ROOTMACROFILE "void autojob()\n{\n";
	    print ROOTMACROFILE " gROOT->SetStyle(\"Plain\");\n;";
	    print ROOTMACROFILE " gROOT->ProcessLine(\".x setenvironment.C\");\n";
	    print ROOTMACROFILE " gROOT->ProcessLine(\".L fittest.C\");\n";
	    print ROOTMACROFILE " gROOT->ProcessLine(\"fittest t\");\n";
	    print ROOTMACROFILE " gROOT->ProcessLine(\"t.Test($isdata,$model,$pdfcomp,$isallcut,$isdepleted,$dumpp,\\\"$pfile\\\",$intp,$ch,$mxcut[$mxi],$mxcut[$mxi+1],$q2cut[$q2i],$q2cut[$q2i+1],$ppluscut[$ppi],$ppluscut[$ppi+1])\");\n";
	    print ROOTMACROFILE "}\n";
	    
	    print JOBFILE "bbrroot -b -q -l autojob.C;\n";
	    print JOBFILE "echo EXIT;\n";
	    
	    close JOBFILE;
	    close ROOTMACROFILE;
	    
	system("chmod 744 tempjob.csh");
	system("./tempjob.csh");
	}
    }
}
exit;


# PARAMETERS FOR fittest::Test():
#  Test(isdata,model (0=2pdf model, 1=3pdfmodel),MC component (0=AllMC, 1=signal MC 2=bkg MC),isallcut,isdepleted,dumpfitparameters,fitparameters file,intpur,b charge,mx low,mx high,q2 low q2 high)

sub usage {
  my($exit, $message) = @_;
  
  print STDERR $message if defined $message;
  print STDERR <<INLINE_LITERAL_TEXT; #'

Usage: $0 <flags> 
    Script that submits jobs for Fits.

flags:
  -help    :       print this message.
  -data    :       fit on data.
  -MC      :       fit on MC.
  -model   :       Choose PDF model: 0=cb+arg, 1 = frank+argus (default 1)
  -bkg     :       fit only Vbc dataset
  -sig     :       fit only VubIN dataset
  -all     :       fit PStar Dataset
  -allcut  :       Apply AllCuts to dataset.
  -lepcut  :       Apply Lepton Cut to dataset.
  -enr     :       fit on enriched sample.
  -depl    :       fit on depleted sample.
  -dumppar :       Dump the fitted parameters.
  -parfile <file>: Use <file> as input parameters configuration file (default: parameters.txt).
  -intpur  :       integrated purity cut.
  -chb     :       1=charged only; 0=neutral only. (default: charged+neutrals).
  -mxl     :       mX low boundary.
  -mxh     :       mX high boundary.
  -q2l     :       q2 low  boundary.
  -q2h     :       q2 high boundary.
  -scanmx  :       do a 1D mx scan.
  -scanmxq2:       do a 2D mx,q2 scan
  -scanpplus:      do a 1D P+ scan
  -haspi0:         0=select BRECO modes w/o pi0; 1=select BRECO modes with pi0; (default: all BRECO modes)

Comments to <petrella\@slac.stanford.edu>
        
Example:  ./submit.pl -MC -sig -lepcut -enr -dumppar -parfile parameters_MC.txt -intpur 0.3 -chb 1 -mxl 0 -mxh 1.55

INLINE_LITERAL_TEXT
#'
  exit($exit) if defined $exit;
}
