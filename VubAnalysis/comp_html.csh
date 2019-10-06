#! /bin/tcsh -f

set dirin = $1
set dirout = $2
set name = $3
set maxvar = 25
set counter = 0
set outfile = $dirout/$name-comp.html

#rm -rf $dirout
#mkdir $dirout
cp $dirin/* $dirout
gzip $dirout/*ps
\rm $outfile
touch $outfile

echo '<table BORDER COLS=4 WIDTH="70%" NOSAVE >' >> $outfile
echo '<tr BGCOLOR="#33CCFF" NOSAVE>' >> $outfile
echo '<td NOSAVE>' $name '</td>' >> $outfile
echo '<td NOSAVE> lepton cuts </td>' >> $outfile
echo '<td NOSAVE> all cuts </td>' >> $outfile
echo '<td NOSAVE> all cuts & breco normalization </td>' >> $outfile
echo '<td NOSAVE> efficiencies (all/lept) </td>' >> $outfile
echo '</tr>' >> $outfile

    
while ($counter < $maxvar)
    
    @ counter ++
    if ( $counter == 1 ) then
	set thevar = nchg
    endif
    if ( $counter == 2 ) then
	set thevar = nneu
    endif
    if ( $counter == 3 ) then
	set thevar = qtot
    endif
    if ( $counter == 4 ) then
	set thevar = nkp
    endif
    if ( $counter == 5 ) then
	set thevar = nks
    endif
    if ( $counter == 6 ) then
	set thevar = pcms
    endif
    if ( $counter == 7 ) then
	set thevar = mm2
    endif
    if ( $counter == 8 ) then
	set thevar = mxhad
    endif
    if ( $counter == 9 ) then
	set thevar = mxhadfit
    endif
    if ( $counter == 10 ) then
	set thevar = nneu80_160
    endif
    if ( $counter == 11 ) then
	set thevar = nneu160_320
    endif
    if ( $counter == 12 ) then
	set thevar = npi0
    endif
    if ( $counter == 13 ) then
	set thevar = kmin
    endif
    if ( $counter == 14 ) then
	set thevar = esneu
    endif
    if ( $counter == 15 ) then
	set thevar = epi0
    endif
    if ( $counter == 16 ) then
	set thevar = efneu
    endif
    if ( $counter == 17 ) then
	set thevar = eneu
    endif
    if ( $counter == 18 ) then
	set thevar = etrk
    endif
    if ( $counter == 19 ) then
	set thevar = intpur
    endif
    if ( $counter == 20 ) then
	set thevar = pur
    endif
    if ( $counter == 21 ) then
	set thevar = modeD0
    endif
    if ( $counter == 22 ) then
	set thevar = modeDc
    endif
    if ( $counter == 23 ) then
	set thevar = modeDstar
    endif
    if ( $counter == 24 ) then
	set thevar = modeDstar0pi0
    endif
    if ( $counter == 25 ) then
	set thevar = modeDstar0gam
    endif
    
    echo '<tr NOSAVE>' >> $outfile
    echo '<td BGCOLOR="#FFCCCC" NOSAVE>' $thevar '</td>' >> $outfile
    echo '<td><a href=comparison7'$thevar'leptoncuts.ps.gz>leptcuts</a></td>' >> $outfile
    echo '<td><a href=comparison7'$thevar'allcuts.ps.gz>allcuts</a></td>' >> $outfile
    echo '<td><a href=comparisonnorm7'$thevar'allcuts.ps.gz>allcuts(breconorm)</a></td>' >> $outfile
    echo '<td><a href=comparisoneff7'$thevar'.ps.gz>ratio eff.</a></td>' >> $outfile
    echo '</tr>' >> $outfile

end 
    echo '</table>' >> $outfile
