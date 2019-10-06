#!/usr/bin/perl

use Getopt::Long;

&GetOptions('int',\$int,'que=s',\$que,'set=s',\$set,'help',\$help);

if($help){
    &usage();
    exit;
}

$ARG{file1};
$ARG{tree1};
$ARG{file2};
$ARG{tree2};
$ARG{cat};
$ARG{outdir};
$ARG{isbch};
$ARG{sys};

$que="medium" if(!$que);
$set="settings.dat" if(!$set);

@pattern=("file1","file2","tree1","tree2","outdir","cat","isbch","sys");

&jobplots();

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
      
      @histos=split(" ",$line);
      print "\n--------> RUNNING @histos[0].tmp\n\n";
      &lastcomp();
      if($int){ #Now job submitting...
	  system("$histos[0].tmp"); 
      } else {
	  system("bsub -q $que -o $ARG{outdir}/logfile_$histos[0].out $histos[0].tmp");
      }
  }
    close(IN);
}
exit;

sub lastcomp(){
    
    if( -e "histos[0].tmp"){
	system("rm $histos[0].tmp");
    }
    system("touch $histos[0].tmp");
    
    $ENV{CWD}=`pwd`;
    chomp($ENV{CWD});
    $ENV{DIRLIB}="$ENV{PWD}/../../shlib/$ENV{BFARCH}";
    
    if($ENV{ARCHNAME} eq "linux"){
	open OUTFILE,"> $histos[0].tmp";
	print OUTFILE "#!/bin/tcsh\n setenv LD_LIBRARY_PATH \${LD_LIBRARY_PATH}:$ENV{DIRLIB}\n";
       	print OUTFILE "$ENV{BFDIST}/releases/$ENV{BFCURRENT}/bin/$ENV{BFARCH}/bbrroot -b -q '$ENV{CWD}/comp.C(\"$ARG{file1}\",\"$ARG{tree1}\",\"$ARG{file2}\",\"$ARG{tree2}\",\"$ARG{outdir}\",\"$histos[1]\",$histos[2],$histos[3],$histos[4],$histos[5],$histos[6],$ARG{isbch},$ARG{cat},$ARG{sys})'\n";
	close OUTFILE;
    }
    system("chmod 744 $histos[0].tmp");
}

sub setvariable(){
    @line=split("=",$_[0]);
    $line[0]=~ s/\s//;   #cleans spacebars
    $line[1]=~ s/\s//;   #cleans spacebars
    $ARG{$line[0]}=$line[1];
}
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






