#!/bin/csh -f

echo ../bin/$BFARCH/anaplot -c $chain1 -D $dir -C $filecut -V $name -min $min -max $max -b $bin

#../bin/$BFARCH/anaplot -c $chain1 -D $dir -C $filecut -V $name -min $min -max $max -b $bin

bsub -C 0 -q kanga -o $dir/$name.out ../bin/$BFARCH/anaplot -c $chain1 -D $dir -C $filecut -V $name -min $min -max $max -b $bin




