#!/bin/csh -f

setenv filecut exclsettings/settingfit_mcfit_etap.dat
setenv isbch 1
setenv outputdir $my_scratch_dir/etap

# SIGNAL
setenv chain1 exclchains/reso
setenv dir $outputdir/reso/
setenv thedir $outputdir/reso/

rm -rf $thedir
mkdir $thedir

jobanal_etap.csh

# GENERIC
setenv chain1 exclchains/generic
setenv dir $outputdir/generic/
setenv thedir $outputdir/generic/

rm -rf $thedir
mkdir $thedir

jobanal_etap.csh

# COCKTAIL 
setenv chain1 exclchains/cocktail 
setenv dir $outputdir/cocktail/  
setenv thedir $outputdir/cocktail/  
  
rm -rf $thedir  
mkdir $thedir  
  
jobanal_etap.csh  

