#! /bin/tcsh -f


#setenv DIRCHAIN datachains2000
#setenv YEAR 2000
#setenv TYPE data
#setenv DIRSTAT '"datachains2000/stat/alldata2000-"'
#setenv DIRROOT '"/u/ec/daniele/scra/alldata2000-"'
#setenv ISDATA 1

set dirchain = $DIRCHAIN
set year = $YEAR
set type = $TYPE
set dirstat = $DIRSTAT
set dirroot = $DIRROOT
set isdata = $ISDATA

countevents.csh $dirchain $type $year
    
set lines = `ls $dirchain | wc  | awk '{print $1}'`
@ lines --
echo $dirchain

\rm temproot

touch temproot

echo #!/bin/sh >> temproot 

echo "ROOTSYS=/afs/slac.stanford.edu/g/babar/package/root/3.01-06/Linux2" >> temproot
echo "export ROOTSYS" >> temproot

echo "LD_LIBRARY_PATH=/afs/slac.stanford.edu/g/babar/package/root/3.01-06/Linux2/lib" >> temproot
echo "export LD_LIBRARY_PATH" >> temproot

echo "exec /afs/slac.stanford.edu/g/babar/package/root/3.01-06/Linux2/bin/root rootlogon_lin.C 'history.C($lines,$dirroot,$dirstat,$isdata)'  -b $*" >> temproot 

chmod u+x temproot

temproot

\rm temproot

