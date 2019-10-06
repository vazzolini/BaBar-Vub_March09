#!/usr/local/bin/perl

my $start = 0;
my $stop =  180;
for($j = $start; $j< $stop; $j+=10){
    my $k = $j+9;
#    print("bsub -q kanga -o ./loggen$j-$k ./toygen -start $j -stop $k -sfd 0 -sfg 0 -sfs 0\n");
    system("bsub -q kanga -o ./loggen$j-$k ./toygen -start $j -stop $k -sfd 0 -sfg 0 -sfs 0\n");
}

exit;
