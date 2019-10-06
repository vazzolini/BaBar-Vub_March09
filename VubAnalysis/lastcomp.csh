#!/bin/csh -f

\rm $thedir/$name-tempcomp.csh
touch $thedir/$name-tempcomp.csh

if ( ! $?sys ) setenv sys 0

if ( $BFARCH == 'Linux2') then
    
    setenv ROOTSYS /afs/slac.stanford.edu/g/babar/package/root/3.01-06/Linux2
    setenv LD_LIBRARY_PATH /afs/slac.stanford.edu/g/babar/package/root/3.01-06/Linux2/lib

    echo "/afs/slac.stanford.edu/g/babar/package/root/3.01-06/Linux2/bin/root -b<<EOF" >> $thedir/$name-tempcomp.csh
    echo " .x rootlogon_lin.C">> $thedir/$name-tempcomp.csh
    echo " .x compS.C($file1,$file2,$dir,$var,$min,$max,$bin,$shift,$smear,$isbch,$cat,$sys)"   >> $thedir/$name-tempcomp.csh
    echo ".q"    >> $thedir/$name-tempcomp.csh

    
else

    setenv ROOTSYS /afs/slac.stanford.edu/g/babar/package/root/3.03-06/SunOS58
    setenv LD_LIBRARY_PATH /afs/slac.stanford.edu/g/babar/package/root/3.03-06/SunOS58/lib


#    echo "/afs/slac.stanford.edu/g/babar/package/root/3.03-06/SunOS58/bin/root -b<<EOF" >> $thedir/$name-tempcomp.csh
    echo "testroot -b<<EOF" >> $thedir/$name-tempcomp.csh
    echo " .x rootlogon_sun.C">> $thedir/$name-tempcomp.csh
    echo " .x comp.C($file1,$file2,$dir,$var,$min,$max,$bin,$shift,$smear,$isbch,$cat,$sys)"   >> $thedir/$name-tempcomp.csh
    echo ".q"    >> $thedir/$name-tempcomp.csh
    echo "EOF"    >> $thedir/$name-tempcomp.csh
    

endif

chmod u+x $thedir/$name-tempcomp.csh
bsub -q long $thedir/$name-tempcomp.csh
#$thedir/$name-tempcomp.csh


#$thedir/$name-tempcomp.csh
#\rm $name-tempcomp.csh
