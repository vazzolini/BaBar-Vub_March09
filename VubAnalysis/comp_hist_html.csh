#! /bin/tcsh -f

set dirin = $1
set dirout = $2
set name = $3
set maxvar = 24
set counter = 0
set outfile = $dirout/$name-comp.html

#rm -rf $dirout
#mkdir $dirout
cp $dirin/* $dirout

\rm $outfile
touch $outfile

echo '<table BORDER COLS=4 WIDTH="70%" NOSAVE >' >> $outfile
echo '<tr BGCOLOR="#33CCFF" NOSAVE>' >> $outfile
echo '<td NOSAVE>' $name '</td>' >> $outfile
echo '<td NOSAVE> lepton cuts </td>' >> $outfile
echo '<td NOSAVE> all cuts </td>' >> $outfile
echo '</tr>' >> $outfile

    
while ($counter < $maxvar)
    
    @ counter ++
    if ( $counter == 1 ) then
	set thevar = pcms
	set theps = a1000.eps
	set thevar2 = pcms_allcuts
	set theps2 = a1020.eps
    endif
    if ( $counter == 2 ) then
	set thevar = pcms_ele
	set theps = e1000.eps
	set thevar2 = pcms_ele_allcuts
	set theps2 = e1020.eps
    endif
     if ( $counter == 3 ) then
	set thevar = pcms_mu
	set theps = m1000.eps
	set thevar2 = pcms_mu_allcuts
	set theps2 = m1020.eps
    endif
     if ( $counter == 4 ) then
	set thevar = theta_ele
	set theps = e1500.eps
	set thevar2 = theta_ele_allcuts
	set theps2 = e1510.eps
    endif
     if ( $counter == 5 ) then
	set thevar = theta_mu
	set theps = m1500.eps
	set thevar2 = theta_mu
	set theps2 = m1510.eps
    endif
    if ( $counter == 6 ) then
	set thevar = trkspectrum
	set theps = a1600.eps
	set thevar2 = trkspectrum_allcuts
	set theps2 = a1610.eps
    endif
    if ( $counter == 7 ) then
	set thevar = trktheta
	set theps = a1700.eps
	set thevar2 = trktheta_allcuts
	set theps2 = a1710.eps
    endif
     if ( $counter == 8 ) then
	set thevar = neuspectrum
	set theps = a1800.eps
	set thevar2 = neuspectrum_allcuts
	set theps2 = a1810.eps
    endif
     if ( $counter == 9 ) then
	set thevar = neutheta
	set theps = a1900.eps
	set thevar2 = neutheta_allcuts
	set theps2 = a1910.eps
    endif
     if ( $counter == 10 ) then
	set thevar = mxhad
	set theps = a2000.eps
	set thevar2 = mxhad_allcuts
	set theps2 = a2020.eps
    endif
    if ( $counter == 11 ) then
	set thevar = mxhadfit
	set theps = a2400.eps
	set thevar2 = mxhadfit_allcuts
	set theps2 = a2420.eps
    endif
    if ( $counter == 12 ) then
	set thevar = mxhadfit_ele_leptoncuts
	set theps = e2400.eps
	set thevar2 = mxhadfit_ele_allcuts
	set theps2 = e2420.eps
    endif
    if ( $counter == 13 ) then
	set thevar = mxhadfit_mu_allcuts
	set theps = m2400.eps
	set thevar2 = mxhadfit_mu_allcuts
	set theps2 = m2420.eps
    endif
    if ( $counter == 14 ) then
	set thevar = mxhadfit2
	set theps = a2800.eps
	set thevar2 = mxhadfit2_allcuts
	set theps2 = a2820.eps
    endif
    if ( $counter == 15 ) then
	set thevar = mxhadfit3
	set theps = a2900.eps
	set thevar2 = mxhadfit3_allcuts
	set theps2 = a2920.eps
    endif
     if ( $counter == 16 ) then
	set thevar = costhetamiss
	set theps = a2600.eps
	set thevar2 = costhetamiss_allcuts
	set theps2 = a2620.eps
    endif
     if ( $counter == 17 ) then
	set thevar = pmiss
	set theps = a2700.eps
	set thevar2 = pmiss_allcuts
	set theps2 = a2720.eps
    endif
    if ( $counter == 18 ) then
	set thevar = mm2
	set theps = a3000.eps
	set thevar2 = mm2_allcuts
	set theps2 = a3020.eps
    endif
     if ( $counter == 19 ) then
	set thevar = ntrk
	set theps = a4000.eps
	set thevar2 = ntrk_allcuts
	set theps2 = a4020.eps
    endif
     if ( $counter == 20 ) then
	set thevar = nneu
	set theps = a4100.eps
	set thevar2 = nneu_allcuts
	set theps2 = a4120.eps
    endif
     if ( $counter == 21 ) then
	set thevar = qtot
	set theps = a4300.eps
	set thevar2 = qtot_allcuts
	set theps2 = a4320.eps
    endif
    if ( $counter == 22 ) then
	set thevar = nkp
	set theps = a4400.eps
	set thevar2 = nkp_allcuts
	set theps2 = a4420.eps
    endif
    if ( $counter == 23 ) then
	set thevar = nks
	set theps = a4500.eps
	set thevar2 = nks_allcuts
	set theps2 = a4520.eps
    endif
     if ( $counter == 24 ) then
	set thevar = nle
	set theps = a4600.eps
	set thevar2 = nle_allcuts
	set theps2 = a4620.eps
    endif
    
    echo '<tr NOSAVE>' >> $outfile
    echo '<td BGCOLOR="#FFCCCC" NOSAVE>' $thevar '</td>' >> $outfile
    echo '<td><a href='$theps'>leptcuts</a></td>' >> $outfile
    echo '<td><a href='$theps2'>allcuts</a></td>' >> $outfile
    echo '</tr>' >> $outfile

end 
    echo '</table>' >> $outfile
