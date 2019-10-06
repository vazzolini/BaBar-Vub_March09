#!/bin/csh -f

setenv cat 7
#setenv file1 '"~/scra/alldata.root"'
#setenv file2 '"~/scra/allcock.root"'

#setenv dir '"scratch/test"'
#setenv thedir scratch/test
#setenv isbch 2

\rm -rf $thedir
mkdir $thedir


jobplots.csh

