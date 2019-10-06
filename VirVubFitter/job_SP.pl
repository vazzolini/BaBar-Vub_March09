#!/usr/local/bin/perl

my $num_f = 1;
my $fin = 100;
while($num_f <= $fin){
    $num_file = $num_f;
   $name = "_SPright_$num_file";
################### 1d MX
####    $truename="mx155_b_63";

   print "./subVirFit.pl -b -d -que kanga -cm CM2 -debug -Sun -novarfit -small -tmodel -notunbinmes=3 -rew 12 -wF sysWd/HFAGCombWin06_3d_0.txt -fixSBratio corrratiosigpeakmx_right.txt -countmc -dssRatio 0.551496 -SPseed `date +%s` -flag $name \n";
   system("./subVirFit.pl -b -d -que kanga -cm CM2 -debug -Sun -novarfit -small -tmodel -notunbinmes=3 -rew 12 -wF sysWd/HFAGCombWin06_3d_0.txt -fixSBratio corrratiosigpeakmx_right.txt -countmc -dssRatio 0.551496 -SPseed `date +%s` -flag $name \n");

#    system("./subVirFit.pl -b -cm CM2 -debug -novarfit -small -tmodel -notunbinmes=3 -que kanga -rew 12 -wF sysWd/HFAGCombWin06_3d_0.txt -Sun -Sys `date +%s` -flag $truename\n");

    sleep 3;
    $num_f++;

}
