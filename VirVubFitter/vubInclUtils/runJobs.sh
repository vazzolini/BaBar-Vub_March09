#! /bin/tcsh

set command = "./subVirFit.pl"
# common flags
set flags = "-d -b -cm CM2 -debug -novarfit -small -dirS=/nfs/farm/babar/AWGsemilep01/menges/summer06/store.2/chains-1111 -tmodel -notunbinmes=3 -que kanga"
set tag   = ""
set suffix = ""

# default fits with HFAG b->clv+b->sg
set rewflags = "-rew 12 -wF sysWd/HFAGCombWin06_3d_0.txt"
set SF = HFAGComb06

# 1d fits: mx
set fittype = "-Sun"
set tag = "Su"
$command $flags $rewflags $fittype        -flag 1DMR14${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -depl  -flag 1DMR14${tag}${SF}dep${suffix}
$command $flags $rewflags $fittype -Run12 -flag 1DMR12${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -Run3  -flag 1DMR3${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -Run4  -flag 1DMR4${tag}${SF}default${suffix}
set fittype = ""
set tag = ""
$command $flags $rewflags $fittype        -flag 1DMR14${SF}default${suffix}

# 1d fits: pplus
set fittype = "-Sun -distfit 2"
set tag = "Su"
$command $flags $rewflags $fittype        -flag 1DPR14${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -depl  -flag 1DPR14${tag}${SF}dep${suffix}
$command $flags $rewflags $fittype -Run12 -flag 1DPR12${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -Run3  -flag 1DPR3${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -Run4  -flag 1DPR4${tag}${SF}default${suffix}
set fittype = "-distfit 2"
set tag = ""
$command $flags $rewflags $fittype        -flag 1DPR14${SF}default${suffix}

# 2d fits
set fittype = "-comb -Sun -mxbin 1.7 -q2bin 8"
set tag = "Su"
$command $flags $rewflags $fittype        -flag 2DR14${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -depl  -flag 2DR14${tag}${SF}dep${suffix}
$command $flags $rewflags $fittype -Run12 -flag 2DR12${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -Run3  -flag 2DR3${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -Run4  -flag 2DR4${tag}${SF}default${suffix}
set fittype = "-comb -mxbin 1.7 -q2bin 8"
set tag = ""
$command $flags $rewflags $fittype        -flag 2DR14${SF}default${suffix}

# 1d fits: mx 3p
set suffix = "3p"
set tag = "Su"
set fittype = "-Sun -opt 1"
$command $flags $rewflags $fittype        -flag 1DMR14${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -depl  -flag 1DMR14${tag}${SF}dep${suffix}
$command $flags $rewflags $fittype -Run12 -flag 1DMR12${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -Run3  -flag 1DMR3${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -Run4  -flag 1DMR4${tag}${SF}default${suffix}

# 1d fits: pplus 3p
set suffix = "3p"
set tag = "Su"
set fittype = "-Sun -opt 1 -distfit 2"
$command $flags $rewflags $fittype        -flag 1DPR14${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -depl  -flag 1DPR14${tag}${SF}dep${suffix}
$command $flags $rewflags $fittype -Run12 -flag 1DPR12${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -Run3  -flag 1DPR3${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -Run4  -flag 1DPR4${tag}${SF}default${suffix}

# 2d fits: 3p
set fittype = "-comb -Sun -mxbin 1.7 -q2bin 8 -opt 1"
set tag = "Su"
set suffix = "3p"
$command $flags $rewflags $fittype        -flag 2DR14${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -depl  -flag 2DR14${tag}${SF}dep${suffix}
$command $flags $rewflags $fittype -Run12 -flag 2DR12${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -Run3  -flag 2DR3${tag}${SF}default${suffix}
$command $flags $rewflags $fittype -Run4  -flag 2DR4${tag}${SF}default${suffix}

set suffix = ""

# default fits with babar b->sg
set rewflags = "-rew 12 -wF sysWd/babar05_3d_0.txt"
set SF = Babar05

set fittype = "-Sun"
set tag = "Su"
$command $flags $rewflags $fittype        -flag 1DMR14${tag}${SF}default${suffix}
set fittype = "-Sun -distfit 2"
$command $flags $rewflags $fittype        -flag 1DPR14${tag}${SF}default${suffix}
set fittype = "-comb -Sun -mxbin 1.7 -q2bin 8"
$command $flags $rewflags $fittype        -flag 2DR14${tag}${SF}default${suffix}

# default fits with babar b->sg 3p
set suffix = "3p"
set fittype = "-Sun -opt 1"
set tag = "Su"
$command $flags $rewflags $fittype        -flag 1DMR14${tag}${SF}default${suffix}
set fittype = "-Sun -distfit 2 -opt 1"
$command $flags $rewflags $fittype        -flag 1DPR14${tag}${SF}default${suffix}
set fittype = "-comb -Sun -mxbin 1.7 -q2bin 8 -opt 1"
$command $flags $rewflags $fittype        -flag 2DR14${tag}${SF}default${suffix}

# default fits with old Belle b->sg
set rewflags = "-rew 11 -wF sysWd/belle3d_sum04_0.inc -FFile wfermifile.dat"
set SF = "Belle04"

set fittype = "-Sun"
set tag = "Su"
$command $flags $rewflags $fittype        -flag 1DMR14${tag}${SF}default${suffix}
set fittype = "-Sun -distfit 2"
$command $flags $rewflags $fittype        -flag 1DPR14${tag}${SF}default${suffix}
set fittype = "-comb -Sun -mxbin 1.7 -q2bin 8"
$command $flags $rewflags $fittype        -flag 2DR14${tag}${SF}default${suffix}


# default fits with HFAG b->clv+b->sg
set flags = "-d -b -cm CM2 -debug -novarfit -dirS=/nfs/farm/babar/AWGsemilep01/menges/summer06/store.2/chains-1111 -tmodel -notunbinmes=3 -que xlong"
set suffix = "FullNTP"
set rewflags = "-rew 12 -wF sysWd/HFAGCombWin06_3d_0.txt"
set SF = HFAGComb06

# 1d fits: mx
set fittype = "-Sun"
$command $flags $rewflags $fittype        -flag 1DMR14${tag}${SF}${suffix}

exit
