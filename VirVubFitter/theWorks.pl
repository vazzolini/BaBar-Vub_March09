#!/usr/local/bin/perl
# we assume we have all dat files in a single place (actually this has been hacked...)
# 

use Getopt::Long;

Getopt::Long::config('bundling_override');
my  %opts = ();
GetOptions( \%opts, 
            "repo|r=s", 
            "log|l=s",
            "theo|t=s");


###my $LOGREPO = "/u/ec/bozzi/Thefiles/dat";
###my $LOGDIR = "sundfn";

######my $LOGREPO = $opts{'repo'} || "/u/ec/bozzi/Thefiles/dat";
my $LOGREPO = $opts{'repo'} || "/u1/Vubresults/";
my $LOGDIR = $opts{'log'} || "sundfn";

@mxcuts = ("15","17","186");
@q2cuts = ("00","02","04","06","08","10","12","14");

@erresys = ("00","01","02","06","07","08","09","10","11","12","13","14","15","16","17","19");

my %thcleo;
my $LOGSUBDIR = $opts{'theo'} || "theocleo";


if($opts{'theo'} =~ /cleo/) {
  @thcleo = ("0","1","2","3","4","5","6","7","8","9");
  $LOGSUBDIR = "theocleo";
} 
if($opts{'theo'} =~ /belle/) {
  @thcleo = ("0","1","2","3","4","5","6","7","8");
  $LOGSUBDIR = "theobelle";
}

my %pstarcorr;

$pstarcorr{"15"}{"00"} = 1.23715/1.22289;
$pstarcorr{"15"}{"02"} = 1.26811/1.25142;
$pstarcorr{"15"}{"04"} = 1.29261/1.27162;
$pstarcorr{"15"}{"06"} = 1.29559/1.28304;
$pstarcorr{"15"}{"08"} = 1.27996/1.27781;
$pstarcorr{"15"}{"10"} = 1.27101/1.28037;
$pstarcorr{"15"}{"12"} = 1.29936/1.27959;
$pstarcorr{"15"}{"14"} = 1.27557/1.25348;

$pstarcorr{"17"}{"00"} = 1.21172/1.20239;
$pstarcorr{"17"}{"02"} = 1.25389/1.24012;
$pstarcorr{"17"}{"04"} = 1.29565/1.27494;
$pstarcorr{"17"}{"06"} = 1.30756/1.28984;
$pstarcorr{"17"}{"08"} = 1.28801/1.2829;
$pstarcorr{"17"}{"10"} = 1.27938/1.28593;
$pstarcorr{"17"}{"12"} = 1.30085/1.28023;
$pstarcorr{"17"}{"14"} = 1.27493/1.25304;

$pstarcorr{"186"}{"00"} = 1.21587/1.20029;
$pstarcorr{"186"}{"02"} = 1.26130/1.23824;
$pstarcorr{"186"}{"04"} = 1.30127/1.27756;
$pstarcorr{"186"}{"06"} = 1.30805/1.29041;
$pstarcorr{"186"}{"08"} = 1.29076/1.28698;
$pstarcorr{"186"}{"10"} = 1.28011/1.28835;
$pstarcorr{"186"}{"12"} = 1.29925/1.27945;
$pstarcorr{"186"}{"14"} = 1.27493/1.25304;


my $LOGEXT = "dat";

my $LOGJOB = "Ibu\*";
my $brbrline;

my %pbfth;
my %pbftemp;
my %strel;
my %mcrel;

### statistical error and MC statistics
### remember to scale

my %dffp;
my %dffm;

foreach $mxc(@mxcuts){
  foreach $q2c(@q2cuts){
    $dffp{$mxc}{$q2c} = 0;
    $dffm{$mxc}{$q2c} = 0;
  }
}

foreach $esys (@thcleo) {
  foreach $mxc (@mxcuts) {
    foreach $q2c (@q2cuts) {
      $LOGJOB = "Ibu\*$esys\_$q2c$mxc";
      if($opts{'theo'} =~ /belle/) {
	$LOGJOB = "Ibu\*SUNdfn$esys\_$q2c$mxc";
      }
      my $LOGFILE = "$LOGREPO/$LOGDIR/$LOGSUBDIR/$LOGJOB/\*$LOGEXT";
      if($LOGDIR =~ /sun/) {
        $brbrline = `grep  "PartialBRBR =" $LOGFILE | awk 'BEGIN{FS="("}{print \$1,\$2}' | awk '{print \$1,\$3,\$5,\$8}'`;
###        $brbrline = `grep  "PartialBRBR     " $LOGFILE | awk '{print \$2,\$3,\$4}'`;
      } else {
        $brbrline = `grep "BRBR" $LOGFILE`;
      }
      print "r",$esys, " mx<",$mxc, " q2>",$q2c, "   ", $brbrline;
      my ($dum,$pbf,$sterr,$mcerr) = split(" ",$brbrline);
      $strel{$esys}{$mxc}{$q2c} = $sterr/$pbf;
      $mcrel{$esys}{$mxc}{$q2c} = $mcerr/$pbf;
      $pbftemp{$esys}{$mxc}{$q2c} = $pbf*$pstarcorr{$mxc}{$q2c};
      print "th ",$esys," mx<",$mxc," q2>",$q2c," PBF and errors: ", $pbftemp{$esys}{$mxc}{$q2c}," ", $strel{$esys}{$mxc}{$q2c}, " ", $mcrel{$esys}{$mxc}{$q2c};
      print "\n";
####
#### compute theoretical systematics here
####
      my $tmpth=($pbftemp{$esys}{$mxc}{$q2c}-$pbftemp{"0"}{$mxc}{$q2c})/$pbftemp{"0"}{$mxc}{$q2c};
      if($tmpth>$dffp{$mxc}{$q2c}) {
	$dffp{$mxc}{$q2c} = $tmpth;
      }
      if($tmpth<$dffm{$mxc}{$q2c}) {
	$dffm{$mxc}{$q2c} = $tmpth;
      }	

##      my $sqsum = sqrt($strel*$strel+$mcrel*$mcrel);
##      print $sqsum; print "\n";
    }
  }
}

### value for the PBF (pstar corrected)

my %spbfth;
my %pmbfth;

my $slbf=0.1083;

foreach $mxc (@mxcuts) {
  foreach $q2c (@q2cuts) {
    $pbfth{$mxc}{$q2c} = $pbftemp{"0"}{$mxc}{$q2c}*$slbf*1000;
    $spbfth{$mxc}{$q2c} = $strel{"0"}{$mxc}{$q2c};
    $mpbfth{$mxc}{$q2c} = $mcrel{"0"}{$mxc}{$q2c};
    
  }
}

# ### theoretical error with 2 points: half of the sum of absolute differences. Symmetrize here (it's negligible)
# 
# my %dff;
# 
# foreach $mxc (@mxcuts) {
#   foreach $q2c (@q2cuts) {
#     $dff{$mxc}{$q2c}  = 0.5*(abs($pbftemp{"1"}{$mxc}{$q2c}-$pbftemp{"0"}{$mxc}{$q2c})/$pbftemp{"0"}{$mxc}{$q2c}+
# 			     abs($pbftemp{"2"}{$mxc}{$q2c}-$pbftemp{"0"}{$mxc}{$q2c})/$pbftemp{"0"}{$mxc}{$q2c});
#     print " ", $mxc, " ", $q2c, " ", $dff{$mxc}{$q2c}; print "\n";
#   }
# }



### systematic errors here r00-r19
### I will quote relative errors only so no pstarfactor correction is needed here 

my %bfsys;
my %sysrel;
my %systot;

foreach $esys (@erresys) {
  foreach $mxc (@mxcuts) {
    foreach $q2c (@q2cuts) {
      print "$esys,$mxc,$q2c\n";
####TEMP CB
      $LOGREPO = "/u/ec/bozzi/Thefiles/dat/";
      my $LOGSUBDIR = "syst";
      $LOGJOB = "Ibu\*$esys\_$q2c$mxc";
      my $LOGFILE = "$LOGREPO/$LOGDIR/$LOGSUBDIR/$LOGJOB/\*$LOGEXT";
      print "$LOGFILE\n"; 
      my $brbrline = `grep "BRBR  " $LOGFILE`;
      my ($dum,$pbf,$sterr,$mcerr) = split(" ",$brbrline);
      $bfsys{$esys}{$mxc}{$q2c} = $pbf;
      $sysrel{$esys}{$mxc}{$q2c} = abs($pbf-$bfsys{"00"}{$mxc}{$q2c})/$bfsys{"00"}{$mxc}{$q2c};
      $systot{$mxc}{$q2c}+=$sysrel{$esys}{$mxc}{$q2c}*$sysrel{$esys}{$mxc}{$q2c};
      print "r",$esys, " mx<",$mxc, " q2>",$q2c, "   ", $bfsys{$esys}{$mxc}{$q2c}, " ", 
	$bfsys{"00"}{$mxc}{$q2c}, " ",$sysrel{$esys}{$mxc}{$q2c}," ",$systot{$mxc}{$q2c};
      print "\n";
    }
  }
}


### total detector systematics
### add 
print "\n";
print "total detector systematics";
print "\n";
foreach $mxc (@mxcuts) {
  foreach $q2c (@q2cuts) {
    $systot{$mxc}{$q2c} = sqrt($systot{$mxc}{$q2c});
    print $systot{$mxc}{$q2c};
    print "\n";
  }
}


### Breco systematics: mesfits 1.4%, mes shape 3% (BAD540), epsUl/epsSLl 3% (BAD540) = 4.5%

my $brecoSys = 0.045;

### Background systematics: binning 7%, Brew 2.2%, Drew 1% = 7.4%

my $bkgSys = 0.074;

#### B reweighting
#### /u1/Vubresults/sundfn/IbuZ_SUNdfn_????_b_xxx

#### D reweighting
#### /u1/Vubresults/sundfn/IbuZ_SUNdfn_????__xxx

### fare media e deviazione standard

### Signal syst: excl 2.7% (BAD540), incl 1.4% (dominique), had 3% (BAD540), ssbar 3.7% (BAD540) = 5.6%

#####my $sigSys = 0.056;
#####only 3.7% from ssbar popping
my $ssSig = 0.037;

#### exclusive
#### /u1/Vubresults/sundfn/theocleo da 11 a 21
#### /u1/Vubresults/sundfn/theobelle da 11 a 21

if($opts{'theo'} =~ /cleo/) {
  @exbr = ("11","12","13","14","15","16","17","18","19","20","21");
} 
if($opts{'theo'} =~ /belle/) {
  @exbr = ("11","12","13","14","15","18","19","20","21");
}

@pippa = ("17");

my %exclSys; 
my %exclRel;
my %hadSys;
foreach $mxc (@pippa) {
  foreach $q2c (@q2cuts) {
    foreach $esys (@exbr) {
      #####      print "$esys,$mxc,$q2c\n";
      $LOGREPO = "/u1/Vubresults";
      ####my $LOGSUBDIR = "syst";
      $LOGJOB = "Ibu\*$esys\_$q2c$mxc";
      my $LOGFILE = "$LOGREPO/$LOGDIR/$LOGSUBDIR/$LOGJOB/\*$LOGEXT";
      #####      print "$LOGFILE\n"; 
      #####my $brbrline = `grep "BRBR  " $LOGFILE`;
      my $brbrline = `grep  "PartialBRBR =" $LOGFILE | awk 'BEGIN{FS="("}{print \$1,\$2}' | awk '{print \$1,\$3,\$5,\$8}'`;
      my ($dum,$pbf,$sterr,$mcerr) = split(" ",$brbrline);
      $pbftemp{$esys}{$mxc}{$q2c} = $pbf;
      $exclRel{$esys}{$mxc}{$q2c} = abs($pbf-$pbftemp{"0"}{$mxc}{$q2c})/$pbftemp{"0"}{$mxc}{$q2c};
      print "excl",$esys, " mx<",$mxc, " q2>",$q2c, "   ", $pbftemp{$esys}{$mxc}{$q2c}, " ", 
	$pbftemp{"0"}{$mxc}{$q2c}, " ",$exclRel{$mxc}{$q2c};
      print "\n";
    }
    #### 14 & 15: pion (cleo) or pion+rho (Belle)
    if($exclRel{"14"}{$mxc}{$q2c} > $exclRel{"15"}{$mxc}{$q2c}) {
      $exclSys{$mxc}{$q2c}+=$exclRel{"14"}{$mxc}{$q2c}*$exclRel{"14"}{$mxc}{$q2c};
    }	
    if($exclRel{"14"}{$mxc}{$q2c} <= $exclRel{"15"}{$mxc}{$q2c}) {
      $exclSys{$mxc}{$q2c}+=$exclRel{"15"}{$mxc}{$q2c}*$exclRel{"15"}{$mxc}{$q2c};
    }
    ### add 16&17 in case of cleo (rho)
    if($LOGSUBDIR =~ /theocleo/) {
      if($exclRel{"16"}{$mxc}{$q2c} > $exclRel{"17"}{$mxc}{$q2c}) {
	$exclSys{$mxc}{$q2c}+=$exclRel{"16"}{$mxc}{$q2c}*$exclRel{"16"}{$mxc}{$q2c};
      }	
      if($exclRel{"16"}{$mxc}{$q2c} <= $exclRel{"17"}{$mxc}{$q2c}) {
	$exclSys{$mxc}{$q2c}+=$exclRel{"17"}{$mxc}{$q2c}*$exclRel{"17"}{$mxc}{$q2c};
      }	
    }
    #### omega
    if($exclRel{"18"}{$mxc}{$q2c} > $exclRel{"19"}{$mxc}{$q2c}) {
      $exclSys{$mxc}{$q2c}+=$exclRel{"18"}{$mxc}{$q2c}*$exclRel{"18"}{$mxc}{$q2c};
    }	
    if($exclRel{"18"}{$mxc}{$q2c} <= $exclRel{"19"}{$mxc}{$q2c}) {
      $exclSys{$mxc}{$q2c}+=$exclRel{"19"}{$mxc}{$q2c}*$exclRel{"19"}{$mxc}{$q2c};
    }
    #### eta + etaprime
    if($exclRel{"20"}{$mxc}{$q2c} > $exclRel{"21"}{$mxc}{$q2c}) {
      $exclSys{$mxc}{$q2c}+=$exclRel{"20"}{$mxc}{$q2c}*$exclRel{"20"}{$mxc}{$q2c};
    }	
    if($exclRel{"20"}{$mxc}{$q2c} <= $exclRel{"21"}{$mxc}{$q2c}) {
      $exclSys{$mxc}{$q2c}+=$exclRel{"21"}{$mxc}{$q2c}*$exclRel{"21"}{$mxc}{$q2c};
    }
    $exclSys{$mxc}{$q2c} = sqrt($exclSys{$mxc}{$q2c});
    print "the exclusive charmless BR uncertainty for mx<$mxc q2>$q2c is $exclSys{$mxc}{$q2c}\n";
    $hadSys{$mxc}{$q2c} = $exclRel{"13"}{$mxc}{$q2c};
    print "the hadronization  uncertainty for mx<$mxc q2>$q2c is $hadSys{$mxc}{$q2c}\n";
  }
}

#### inclusive
#### /u1/Vubresults/sundfn/theocleo 075 125
#### /u1/Vubresults/sundfn/theobelle 075 125

my %inclRelp;
my %inclRelm;
my %inclRel;

foreach $mxc (@pippa) {
  foreach $q2c (@q2cuts) {
      #####      print "$esys,$mxc,$q2c\n";
      $LOGREPO = "/u1/Vubresults";
      ####my $LOGSUBDIR = "syst";
      $LOGJOB = "Ibu\*075\_$q2c$mxc";
      my $LOGFILE = "$LOGREPO/$LOGDIR/$LOGSUBDIR/$LOGJOB/\*$LOGEXT";
      #####      print "$LOGFILE\n"; 
      #####my $brbrline = `grep "BRBR  " $LOGFILE`;
      my $brbrline = `grep  "PartialBRBR =" $LOGFILE | awk 'BEGIN{FS="("}{print \$1,\$2}' | awk '{print \$1,\$3,\$5,\$8}'`;
      my ($dum,$pbf,$sterr,$mcerr) = split(" ",$brbrline);
      $inclRelp{$mxc}{$q2c} = abs($pbf-$pbftemp{"0"}{$mxc}{$q2c})/$pbftemp{"0"}{$mxc}{$q2c};
      print "inclp mx<",$mxc, " q2>",$q2c, "   ", $pbf, " ", 
	$pbftemp{"0"}{$mxc}{$q2c}, " ",$inclRelp{$mxc}{$q2c};
      print "\n";
    }
}

foreach $mxc (@pippa) {
  foreach $q2c (@q2cuts) {
      #####      print "$esys,$mxc,$q2c\n";
      $LOGREPO = "/u1/Vubresults";
      ####my $LOGSUBDIR = "syst";
      $LOGJOB = "Ibu\*125\_$q2c$mxc";
      my $LOGFILE = "$LOGREPO/$LOGDIR/$LOGSUBDIR/$LOGJOB/\*$LOGEXT";
      #####      print "$LOGFILE\n"; 
      #####my $brbrline = `grep "BRBR  " $LOGFILE`;
      my $brbrline = `grep  "PartialBRBR =" $LOGFILE | awk 'BEGIN{FS="("}{print \$1,\$2}' | awk '{print \$1,\$3,\$5,\$8}'`;
      my ($dum,$pbf,$sterr,$mcerr) = split(" ",$brbrline);
      $inclRelm{$mxc}{$q2c} = abs($pbf-$pbftemp{"0"}{$mxc}{$q2c})/$pbftemp{"0"}{$mxc}{$q2c};
      print "inclm mx<",$mxc, " q2>",$q2c, "   ", $pbf, " ", 
	$pbftemp{"0"}{$mxc}{$q2c}, " ",$inclRelm{$mxc}{$q2c};
      print "\n";
    }
}

foreach $mxc (@pippa) {
  foreach $q2c (@q2cuts) {
    $inclRel{$mxc}{$q2c} = 0.5*($inclRelp{$mxc}{$q2c}+$inclRelm{$mxc}{$q2c});
      print "The inclusive charmless BR uncertainty for mx<$mxc q2>$q2c is $inclRel{$mxc}{$q2c}\n";
    }
}
    
my %sigSys;

foreach $mxc (@pippa) {
  foreach $q2c (@q2cuts) {
    $sigSys{$mxc}{$q2c} = sqrt($ssSig*$ssSig+
			       $exclSys{$mxc}{$q2c}*$exclSys{$mxc}{$q2c}+
			       $hadSys{$mxc}{$q2c}*$hadSys{$mxc}{$q2c}+
			       $inclRel{$mxc}{$q2c}*$inclRel{$mxc}{$q2c});
    print "The total signal uncertainty for mx<$mxc q2>$q2c is $sigSys{$mxc}{$q2c}\n";
  }
}


#### HAD
#### /u1/Vubresults/sundfn/syst/ HAD   this is already included in exclusive point 13
### tot of the above 3: 10.31%

my %totTotp;
my %totTotm;
my %absStat;
my %absDet;
my %absBreco;
my %absTheop;
my %absTheom;
my %absBkg;
my %absSig;
my %absMcStat;
my %absTotp;
my %absTotm;

foreach $mxc (@mxcuts) {
  foreach $q2c (@q2cuts) {
    $totTotp{$mxc}{$q2c} = sqrt($systot{$mxc}{$q2c}*$systot{$mxc}{$q2c} 
			       + $dffp{$mxc}{$q2c}*$dffp{$mxc}{$q2c} 
			       + $spbfth{$mxc}{$q2c}*$spbfth{$mxc}{$q2c} 
			       + $mpbfth{$mxc}{$q2c}*$mpbfth{$mxc}{$q2c} 
			       + $brecoSys*$brecoSys
			       + $bkgSys*$bkgSys
			       + $sigSys{$mxc}{$q2c}*$sigSys{$mxc}{$q2c}
			      );
    $totTotm{$mxc}{$q2c} = sqrt($systot{$mxc}{$q2c}*$systot{$mxc}{$q2c} 
			       + $dffm{$mxc}{$q2c}*$dffm{$mxc}{$q2c} 
			       + $spbfth{$mxc}{$q2c}*$spbfth{$mxc}{$q2c} 
			       + $mpbfth{$mxc}{$q2c}*$mpbfth{$mxc}{$q2c} 
			       + $brecoSys*$brecoSys
			       + $bkgSys*$bkgSys
			       + $sigSys{$mxc}{$q2c}*$sigSys{$mxc}{$q2c}
			      );
    $absStat{$mxc}{$q2c} = $pbfth{$mxc}{$q2c}*$spbfth{$mxc}{$q2c};
    $absDet{$mxc}{$q2c} = $pbfth{$mxc}{$q2c}*$systot{$mxc}{$q2c};
    $absBreco{$mxc}{$q2c} = $pbfth{$mxc}{$q2c}*$brecoSys;
    $absBkg{$mxc}{$q2c} = $pbfth{$mxc}{$q2c}*$bkgSys;
    $absTheop{$mxc}{$q2c} = $pbfth{$mxc}{$q2c}*$dffp{$mxc}{$q2c};
    $absTheom{$mxc}{$q2c} = $pbfth{$mxc}{$q2c}*$dffm{$mxc}{$q2c};
    $absSig{$mxc}{$q2c} = $pbfth{$mxc}{$q2c}*$sigSys{$mxc}{$q2c};
    $absMcStat{$mxc}{$q2c} = $pbfth{$mxc}{$q2c}*$mpbfth{$mxc}{$q2c};
    
    $absTotp{$mxc}{$q2c} = $pbfth{$mxc}{$q2c}*$totTotp{$mxc}{$q2c};
    $absTotm{$mxc}{$q2c} = $pbfth{$mxc}{$q2c}*$totTotm{$mxc}{$q2c};

    my $dettaglio = sprintf "Mx<%f Q2>%f PBF = %.3f +/- %.3f (stat) +/- %.3f (det sys) 
+/- %.3f (Breco) +/- %.3f (Bkg) %.3f +%.3f (theory) +/- %.3f (sig) +/- %.3f (MC stat)", 
  $mxc,$q2c, $pbfth{$mxc}{$q2c}, $absStat{$mxc}{$q2c}, $absDet{$mxc}{$q2c}, $absBreco{$mxc}{$q2c}, $absBkg{$mxc}{$q2c},
    $absTheom{$mxc}{$q2c},$absTheop{$mxc}{$q2c},$absSig{$mxc}{$q2c},$absMcStat{$mxc}{$q2c};

    print $dettaglio; 
    print "\n";
    print "\n";
    my $meanerr = 0.5*($absTotp{$mxc}{$q2c}+$absTotm{$mxc}{$q2c});
    my $risultato = sprintf "...or in other words:  Mx<%f Q2>%f = %.3f +/- %.3f ",$mxc,$q2c,,$pbfth{$mxc}{$q2c},$meanerr;
###$absTotm{$mxc}{$q2c},$absTotp{$mxc}{$q2c};

    print $risultato;
    print "\n";
    print "\n";
    print "\n";
  }
}
			      
			      
