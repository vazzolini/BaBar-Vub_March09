#!/usr/bin/perl

use Getopt::Long;

&GetOptions('int',\$int,'que=s',\$que,'set=s',\$set,'help',\$help);

if($help){
    &usage();
    exit;
}

$ARG{rootfile};
$ARG{dataname};
$ARG{mcname};
$ARG{outdir};
$ARG{isdata};
$ARG{bch};

$que="kanga" if(!$que);
$set="settings.dat" if(!$set);

my $tmpfile;

@pattern=("rootfile","outdir","dataname","mcname","bch","isdata");

&jobplots();

## -------------------------------------------------------- ##
sub jobplots() {
    
    open(IN,"$set"); #open file with settings
    print "-------> Echoing $set\n";
  LINE:  while(<IN>){
      chomp($_);
      my $line=$_;
      
      print "$line\n";
      next LINE if(($_ =~ /^\#+/) || ($_ =!/^./));
      foreach $i(@pattern) {
	  if($line =~ /^$i/){
	      &setvariable($line);
	      next LINE;
	  }
      }
      system("mkdir $ARG{outdir}") if(! -e $ARG{outdir} || ! -d $ARG{outdir}); #create outdir if doesn't exist
      system("rm -f $ARG{outdir}/logfile_$histos[0].out") if(-e "$ARG{outdir}/logfile_$histos[0].out"); #clean old logfile
      
      @histos = split(" ",$line);

      if( $ARG{isdata} ) {
	  $tmpfile = "$histos[0]_3";
      } else  {
	  $tmpfile = "$histos[0]_3MC";
      }
            
      print "\n--------> RUNNING $tmpfile.tmp\n\n";

      &lastcomp();
      
      if($int){ #Now job submitting...
	  system("./$tmpfile.tmp"); 
      } else {
	  system("bsub -q $que -o $ARG{outdir}/logfile_$tmpfile.out $tmpfile.tmp");
      }
  }
    close(IN);
}
exit;
## -------------------------------------------------------- ##
sub lastcomp(){
    
    if( -e "$tmpfile.tmp"){
	system("rm $tmpfile.tmp");
    }
    system("touch $tmpfile.tmp");
    
    $ENV{CWD}=`pwd`;
    chomp($ENV{CWD});
    $ENV{DIRLIB}="$ENV{PWD}/../../shlib/$ENV{BFARCH}:$ENV{PWD}/../../RooFitCore/tmp:$ENV{PWD}/../../RooFitModels/tmp";
    
    if($ENV{ARCHNAME} eq "linux"){
	open OUTFILE,"> $tmpfile.tmp";
	print OUTFILE "#!/bin/tcsh\n setenv LD_LIBRARY_PATH \${LD_LIBRARY_PATH}:$ENV{DIRLIB}\n";
	print OUTFILE "./dataMCcomp $ARG{rootfile} $ARG{outdir} $ARG{bch} $ARG{isdata} $histos[0] $histos[1] $histos[2] $histos[3] $histos[4] $histos[5] $histos[6]\n";
	close OUTFILE;
    }
    system("chmod 744 $tmpfile.tmp");
}
## -------------------------------------------------------- ##
sub setvariable(){
    @line=split("=",$_[0]);
    $line[0]=~ s/\s//;   #cleans spacebars
    $line[1]=~ s/\s//;   #cleans spacebars
    $ARG{$line[0]}=$line[1];
}
## -------------------------------------------------------- ##
sub usage {
  
    print STDERR
"Options:
  -int       :       Run interactive mode
  -que       :       Specify queue type
  -set <filename> :  Specift settings filename
  -help      :       Print this help
 
Example:
   ./compare.pl -que long -set mysetting.dat
 "
}
