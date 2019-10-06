#! /bin/tcsh

# common flags
set flags = "-d -b -cm CM2 -debug -novarfit -small"
set tag = ""

# default fits with old Belle b->sg
set rewflags = "-rew 11 -wF sysWd/belle3d_sum04_0.inc -FFile wfermifile.dat"
set SF = Belle04

# 1d fits: mx
set fittype = ""
time ./subVirFit.pl $flags $rewflags $fittype        -flag 1DMR14${SF}test$tag
set fittype = "-Sun"
time ./subVirFit.pl $flags $rewflags $fittype        -flag 1DMR14Su${SF}test$tag

# 1d fits: pplus
set fittype = "-distfit 2"
time ./subVirFit.pl $flags $rewflags $fittype        -flag 1DPR14${SF}test$tag
set fittype = "-Sun -distfit 2"
time ./subVirFit.pl $flags $rewflags $fittype        -flag 1DPR14Su${SF}test$tag

# 2d fits
set fittype = "-comb -Sun -mxbin 1.7 -q2bin 8"
time ./subVirFit.pl $flags $rewflags $fittype        -flag 2DR14Su${SF}test$tag
set fittype = "-comb -mxbin 1.7 -q2bin 8"
time ./subVirFit.pl $flags $rewflags $fittype        -flag 2DR14${SF}test$tag

exit
