#!/bin/tcsh -m

source envfit.csh

touch  $DIRSETTING/mysettings.dat_$TEST
\rm  $DIRSETTING/mysettings.dat_$TEST
touch  $DIRSETTING/mysettings.dat_$TEST
 
touch  $DIRFILES/myfiles.dat_$TEST
\rm  $DIRFILES/myfiles.dat_$TEST
touch  $DIRFILES/myfiles.dat_$TEST  


echo totalStat         $TOTSTAT  >>   $DIRSETTING/mysettings.dat_$TEST
echo totalStatModel    $TOTSTATMOD >>   $DIRSETTING/mysettings.dat_$TEST
echo BRRatioGenValue   $BRRATIOGEN >>   $DIRSETTING/mysettings.dat_$TEST   
echo BRRatioValueTail  $BRRATIOTAIL >>   $DIRSETTING/mysettings.dat_$TEST   
echo pstarfact         $PSTARCOR >>   $DIRSETTING/mysettings.dat_$TEST   
echo mxCut             $MXCUT >>   $DIRSETTING/mysettings.dat_$TEST   
echo useCB             $USECB >>   $DIRSETTING/mysettings.dat_$TEST   
echo fixMeanValue      $FIXMEAN >>   $DIRSETTING/mysettings.dat_$TEST   
echo fixSigma          $FIXSIGMA >>   $DIRSETTING/mysettings.dat_$TEST   
echo fixArgus1         $FIXARG1 >>   $DIRSETTING/mysettings.dat_$TEST   
echo fixArgus2         $FIXARG >>   $DIRSETTING/mysettings.dat_$TEST   
echo fixCB1            $FIXCB1 >>   $DIRSETTING/mysettings.dat_$TEST   
echo fixCB2            $FIXCB2 >>   $DIRSETTING/mysettings.dat_$TEST   
echo leptonPCut        $LEPTCUT >>   $DIRSETTING/mysettings.dat_$TEST   
echo mnuSqLow          $MNULOW >>   $DIRSETTING/mysettings.dat_$TEST   
echo mnuSqHigh         $MNUHIGH >>   $DIRSETTING/mysettings.dat_$TEST   
echo chLow             $CHLOW >>   $DIRSETTING/mysettings.dat_$TEST   
echo chHigh            $CHHIGH >>   $DIRSETTING/mysettings.dat_$TEST   
echo depl              $DEPL >>   $DIRSETTING/mysettings.dat_$TEST   
echo Btype             $BTYPE >>   $DIRSETTING/mysettings.dat_$TEST   
echo lepttype          $LEPTTYPE >>   $DIRSETTING/mysettings.dat_$TEST   
echo fittotshape       $FITTOTSHAPE >>   $DIRSETTING/mysettings.dat_$TEST   
echo mixcorr           $MIXCORR >>   $DIRSETTING/mysettings.dat_$TEST   
echo fitMC             $FITMC >>   $DIRSETTING/mysettings.dat_$TEST   
echo multifit          $MULTIFIT >>   $DIRSETTING/mysettings.dat_$TEST   
echo blinding          $BLIND >>   $DIRSETTING/mysettings.dat_$TEST 
echo blindsize         $BLINDSIZE >>   $DIRSETTING/mysettings.dat_$TEST 
echo randomseed        $RANDOMSEED >>   $DIRSETTING/mysettings.dat_$TEST 
echo issmearAll        $ISSMEARALL >>   $DIRSETTING/mysettings.dat_$TEST 
echo smearAllMeanValue $SMEARALLMEANVALUE >>   $DIRSETTING/mysettings.dat_$TEST 
echo smearAllSigma     $SMEARALLSIGMA >>   $DIRSETTING/mysettings.dat_$TEST 
echo issmearBkg        $ISSMEARBKG >>   $DIRSETTING/mysettings.dat_$TEST 
echo smearBkgMeanValue $SMEARBKGMEANVALUE >>   $DIRSETTING/mysettings.dat_$TEST 
echo smearBkgSigma     $SMEARBKGSIGMA >>   $DIRSETTING/mysettings.dat_$TEST 
echo dotrkreweight     $DOTRKWEIGHT >>   $DIRSETTING/mysettings.dat_$TEST 
echo doneureweight     $DONEUWEIGHT >>   $DIRSETTING/mysettings.dat_$TEST 
echo doBdecreweight    $DOBDECWEIGHT >>   $DIRSETTING/mysettings.dat_$TEST 
echo doDdecreweight    $DODDECWEIGHT >>   $DIRSETTING/mysettings.dat_$TEST 
echo dobrecoreweight   $DOBRECOWEIGHT >>   $DIRSETTING/mysettings.dat_$TEST 
echo dopstarreweight   $DOPSTARWEIGHT >>   $DIRSETTING/mysettings.dat_$TEST 
echo domm2reweight     $DOMM2WEIGHT >>   $DIRSETTING/mysettings.dat_$TEST 
echo minintpur         $MININTPUR  >>   $DIRSETTING/mysettings.dat_$TEST
echo maxintpur         $MAXINTPUR  >>   $DIRSETTING/mysettings.dat_$TEST
echo run               $RUN  >>   $DIRSETTING/mysettings.dat_$TEST
echo nnpi0             $CUTNNPI0  >>   $DIRSETTING/mysettings.dat_$TEST   
echo doleptplot        $LEPTPLOT  >>   $DIRSETTING/mysettings.dat_$TEST  
echo deltamb           $DELTAMB  >>   $DIRSETTING/mysettings.dat_$TEST
echo deltaa            $DELTAA  >>   $DIRSETTING/mysettings.dat_$TEST
echo fitOption         $FITOPTION >>   $DIRSETTING/mysettings.dat_$TEST
#echo userbinning       $USERBINNING  >>   $DIRSETTING/mysettings.dat_$TEST  
#echo bin1              $BIN1  >>   $DIRSETTING/mysettings.dat_$TEST  
#echo bin2              $BIN2  >>   $DIRSETTING/mysettings.dat_$TEST  
#echo bin3              $BIN3  >>   $DIRSETTING/mysettings.dat_$TEST  
#echo bin4              $BIN4  >>   $DIRSETTING/mysettings.dat_$TEST  
#echo bin5              $BIN5  >>   $DIRSETTING/mysettings.dat_$TEST  
#echo bin6              $BIN6  >>   $DIRSETTING/mysettings.dat_$TEST  
#echo bin7              $BIN7  >>   $DIRSETTING/mysettings.dat_$TEST  
#echo bin8              $BIN8  >>   $DIRSETTING/mysettings.dat_$TEST  
#echo bin9              $BIN9  >>   $DIRSETTING/mysettings.dat_$TEST  
#echo bin10             $BIN10  >>   $DIRSETTING/mysettings.dat_$TEST  
echo   >>   $DIRSETTING/mysettings.dat_$TEST 

echo fileVubTotal $FILEVUBTOTAL >>  $DIRFILES/myfiles.dat_$TEST 
echo fileVubDstar $FILEVUBDSTAR >>  $DIRFILES/myfiles.dat_$TEST 
echo fileVubDc $FILEVUBDC >>  $DIRFILES/myfiles.dat_$TEST 
echo fileVubDstar0 $FILEVUBDSTAR0 >>  $DIRFILES/myfiles.dat_$TEST 
echo fileVubD0 $FILEVUBD0 >>  $DIRFILES/myfiles.dat_$TEST 
echo fileVcb $FILEVCB >>  $DIRFILES/myfiles.dat_$TEST 
echo fileData $FILEDATA >>  $DIRFILES/myfiles.dat_$TEST 
echo  >>  $DIRFILES/myfiles.dat_$TEST 

rd $DIR
mkdir $DIR
bsub -q long -o $DIR/vubfit.out ../bin/SunOS58/VubFit -F $DIRFILES/myfiles.dat_$TEST -C $DIRSETTING/mysettings.dat_$TEST -D $DIR -P $TEST  -Mes mesparsetting.dat 
#-Fermi

   
