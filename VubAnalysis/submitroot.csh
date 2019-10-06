#! /bin/tcsh -f

set input = $1
set output = $2
set queue = $3
set MC = $4
set options = '-d 4 -MC -pk -SN'
if ($MC == 0) then
set options = '-d 4 '
endif


foreach i ( $input )

   set name = `echo $i |awk -F/ '{print $NF}'`
   bsub  -q $queue -o $output/$name ../bin/SunOS58/anaQA -t h1 -c  $i -C b0cuts.dat -D $output $options
# ../bin/Linux2/anaQA -t h1 -c  $i -C b0cuts.dat -D $output $options
   echo bsub  -q $queue -o $output/$name ../bin/SunOS58/anaQA -t h1 -c  $i -C b0cuts.dat -D $output $options

end




