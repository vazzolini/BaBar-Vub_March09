#!/bin/csh -f

setenv filecut exclsettings/settingfit_mcfit_eta.dat
setenv isbch 1
setenv outputdir $my_scratch_dir/eta

# SIGNAL
setenv chain1 exclchains/reso
setenv dir $outputdir/reso/
setenv thedir $outputdir/reso/

rm -rf $thedir
mkdir $thedir

jobanal_eta.csh

# GENERIC
setenv chain1 exclchains/generic
setenv dir $outputdir/generic/
setenv thedir $outputdir/generic/

rm -rf $thedir
mkdir $thedir

jobanal_eta.csh

# COCKTAIL
setenv chain1 exclchains/cocktail
setenv dir $outputdir/cocktail/ 
setenv thedir $outputdir/cocktail/ 
 
rm -rf $thedir 
mkdir $thedir 
 
jobanal_eta.csh 

