#! /bin/tcsh -f

set dirin = $1
set dirout = $2
set name = $3
set maxvar = 15
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
echo '<td NOSAVE> lepton cuts (vubcomp) </td>' >> $outfile
echo '<td NOSAVE> lepton cuts (allcomp) </td>' >> $outfile
echo '<td NOSAVE> signal (leptcuts) </td>' >> $outfile
echo '<td NOSAVE> vub (leptcuts) </td>' >> $outfile
echo '<td NOSAVE> vcb (leptcuts) </td>' >> $outfile
echo '<td NOSAVE> other (leptcuts) </td>' >> $outfile
echo '<td NOSAVE> all cuts (vubcomp) </td>' >> $outfile
echo '<td NOSAVE> all cuts (allcomp) </td>' >> $outfile
echo '<td NOSAVE> signal (allcuts) </td>' >> $outfile
echo '<td NOSAVE> vub (allcuts) </td>' >> $outfile
echo '<td NOSAVE> vcb (allcuts) </td>' >> $outfile
echo '<td NOSAVE> other (allcuts) </td>' >> $outfile
echo '</tr>' >> $outfile

    
while ($counter < $maxvar)
    
    @ counter ++
    if ( $counter == 1 ) then
	set thevar = mes
    endif
    if ( $counter == 2 ) then
	set thevar = nrecoomega
    endif
    if ( $counter == 3 ) then
	set thevar = momega
    endif
    if ( $counter == 4 ) then
	set thevar = mm2omega
    endif
    if ( $counter == 5 ) then
	set thevar = nrecopi0
    endif
    if ( $counter == 6 ) then
	set thevar = mom1piome
    endif
    if ( $counter == 7 ) then
	set thevar = mom2piome
    endif
    if ( $counter == 8 ) then
	set thevar = truemom1piome
    endif
    if ( $counter == 9 ) then
	set thevar = truemom2piome
    endif
    if ( $counter == 10 ) then
	set thevar = costhome
    endif
    if ( $counter == 11 ) then
	set thevar = pcms
    endif
    if ( $counter == 12 ) then
	set thevar = mm2
    endif
    if ( $counter == 13 ) then
	set thevar = nchg
    endif
    if ( $counter == 14 ) then
	set thevar = nneu
    endif
    if ( $counter == 15 ) then
	set thevar = qtot
    endif


    
    echo '<tr NOSAVE>' >> $outfile
    echo '<td BGCOLOR="#FFCCCC" NOSAVE>' $thevar '</td>' >> $outfile
    echo '<td><a href='$thevar'1-vubcomp.eps.gz>leptcuts</a></td>' >> $outfile
    echo '<td><a href='$thevar'1-allcomp.eps.gz>leptcuts</a></td>' >> $outfile
    echo '<td><a href='$thevar'1-signal.eps.gz>plot</a></td>' >> $outfile
    echo '<td><a href='$thevar'1-vub.eps.gz>plot</a></td>' >> $outfile
    echo '<td><a href='$thevar'1-vcb.eps.gz>plot</a></td>' >> $outfile
    echo '<td><a href='$thevar'1-other.eps.gz>plot</a></td>' >> $outfile
    echo '<td><a href='$thevar'0-vubcomp.eps.gz>allcuts</a></td>' >> $outfile
    echo '<td><a href='$thevar'0-allcomp.eps.gz>allcuts</a></td>' >> $outfile
    echo '<td><a href='$thevar'0-signal.eps.gz>plot</a></td>' >> $outfile
    echo '<td><a href='$thevar'0-vub.eps.gz>plot</a></td>' >> $outfile
    echo '<td><a href='$thevar'0-vcb.eps.gz>plot</a></td>' >> $outfile
    echo '<td><a href='$thevar'0-other.eps.gz>plot</a></td>' >> $outfile
    echo '</tr>' >> $outfile

end 
    echo '</table>' >> $outfile
