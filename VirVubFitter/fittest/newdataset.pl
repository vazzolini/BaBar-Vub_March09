#!/usr/bin/perl

my @run=("12","3","4","14","5","15");

#my @run=("5");
#my $j=1; #ONLY MC

for ($j=0;$j<2;$j++) {  # 0=DATA, 1=Montecarlo
    foreach $i(@run){
	my $filename="createdset_$j\_$i\.C";
	my $log="createdsetR18_$j\_$i\.log";
	print "$filename $log\n";
	
	open OUTFILE,"> $filename";
	print OUTFILE "\{\n gROOT->ProcessLine(\".x setenvironment.C\");\n";
	print OUTFILE " gROOT->ProcessLine(\".L fittest.C\");\n";
	print OUTFILE " gROOT->ProcessLine(\"fittest t\");\n";
	print OUTFILE " gROOT->ProcessLine(\"t.WriteDataSet($j,$i)\");\n\}\n";
	close OUTFILE;
#	print ("bsub -q kanga -o ~/scra/$log bbrroot -b -q $filename\n");
	system("bsub -q kanga -o ~/scra/$log bbrroot -b -q $filename\n");

    }
}
