#! /bin/tcsh -f

set dir = $1
set type = $2
set year = $3
set dirstat = $dir/stat
rm -rf  $dirstat
mkdir $dirstat

foreach i ( $dir/* )
    echo $i
    set lines = `wc $i | awk '{print $1}'`
    set runfile =  $i-runs
    set lumifile = $i-lumi
    set eventsfile = $i-totevents

    touch $runfile
    set toteve = 0
    set counter = 0
    while ($counter<$lines) 
	@ counter ++    
	if($type == 'data') then
	    set runnumber = `head -$counter $i | tail -1 | awk -F- '{print $3}' | awk -F. '{print $1}'`
	    set tcldir =  /nfs/farm/babar/AWG12/ISL/tcls_sum02/data/$year
	endif
	if($type == 'cock') then
	    set runnumbertemp = `head -$counter $i | tail -1 | grep b0 | awk -F- '{print $5}' | awk -F. '{print $1}'`
	    echo $runnumbertemp
	    if ($runnumbertemp) then
		set runnumber = $runnumbertemp
		set tcldir =  /nfs/farm/babar/AWG12/ISL/tcls_sum02/b0cocktail
	    endif
	    set runnumbertemp = `head -$counter $i | tail -1 | grep bp | awk -F- '{print $5}' | awk -F. '{print $1}'`
	    if ($runnumbertemp) then
		set runnumber = $runnumbertemp
		set tcldir =  /nfs/farm/babar/AWG12/ISL/tcls_sum02/bpcocktail
	    endif	    
	endif
	if($type == 'gene') then
	    set runnumbertemp = `head -$counter $i | tail -1 | grep new | grep -i b0 | awk -F- '{print $5}' | awk -F. '{print $1}'`
	    if ($runnumbertemp) then
		set runnumber = $runnumbertemp
		set tcldir =  /nfs/farm/babar/AWG12/ISL/tcls_sum02/genSP4/B0/$year/new/genbnu
	    endif
	    set runnumbertemp = `head -$counter $i | tail -1 | grep old | grep -i b0 | awk -F- '{print $5}' | awk -F. '{print $1}'`
	    if ($runnumbertemp) then
		set runnumber = $runnumbertemp
		set tcldir =  /nfs/farm/babar/AWG12/ISL/tcls_sum02/genSP4/B0/$year/old/genbnu
	    endif
	    set runnumbertemp = `head -$counter $i | tail -1 | grep old | grep -i chb | awk -F- '{print $5}' | awk -F. '{print $1}'`
	    if ($runnumbertemp) then
		set runnumber = $runnumbertemp
		set tcldir =  /nfs/farm/babar/AWG12/ISL/tcls_sum02/genSP4/ChB/$year/new/genbch
	    endif
	    set runnumbertemp = `head -$counter $i | tail -1 | grep new | grep -i chb | awk -F- '{print $5}' | awk -F. '{print $1}'`
	    if ($runnumbertemp) then
		set runnumber = $runnumbertemp
		set tcldir =  /nfs/farm/babar/AWG12/ISL/tcls_sum02/genSP4/ChB/$year/old/genbch
	    endif	    
	endif	

	set tclfile = `ls $tcldir/*.tcl |grep $year | grep '\-'$runnumber.tcl`
	set nruns = `grep run $tclfile| wc | awk '{print $1}'`
	set counter2 = 0
	while ($counter2<$nruns)
	    @ counter2 ++
	    set therun = `grep run $tclfile|head -$counter2 | tail -1 | awk -F: '{print $1}' | awk '{print $3}'`
	    set theevents = `grep run $tclfile|head -$counter2 | tail -1 | awk -F: '{print $2}' | awk -F/ '{print $1}'`
	    echo $therun >> $runfile 	  
	    @ toteve = $toteve + $theevents 
	end
    end
    echo $toteve> $eventsfile
    set totlumi = `lumi -f $runfile | grep Lumi|grep Proc| awk '{print $3}'`
    echo $totlumi > $lumifile
    echo total events $toteve
    echo total luminosity $totlumi
    mv $runfile $dirstat
    mv $eventsfile $dirstat
    mv $lumifile $dirstat
end
