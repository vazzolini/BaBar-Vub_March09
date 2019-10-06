#!/bin/csh -f

setenv file1 /u/ec/daniele/scra/excl-prod/cocktail.root
setenv filecut settingfit_mcfit_excl.dat

setenv dir scratch/test/
setenv thedir scratch/test/
setenv isbch 1

\rm -rf $thedir
mkdir $thedir

jobanal.csh


setenv file1 /u/ec/daniele/scra/excl-prod/signal-mix.root

setenv dir scratch/testsig/
setenv thedir scratch/testsig/
setenv isbch 1

\rm -rf $thedir
mkdir $thedir

jobanal.csh

