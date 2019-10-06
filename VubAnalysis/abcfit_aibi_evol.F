C
C setup y0 parameters and their errors (dy0/correlation)
C
      SUBROUTINE AIBI_EVOL(PARIDIN,ITEVOLKMN,NJET,NUP,P_COR,Y0,DY0
     $     ,V,IERR)
*****************************************************
* Routine to read and return covariance matrix and starting values
* of the fit parameters depending on choice of parameterisation
*****************************************************
      IMPLICIT NONE

#include "abcfit_setup.inc"
#include "abcfit_aibi_evol.inc"
#include "abcfit_bmatrix.inc"
#include "abcfit.inc"

      INTEGER PARIDIN(NPARTICLES),ITEVOLKMN(3),NJET,NUP,IERR
      REAL P_COR(4,NPARTICLES)
      DOUBLE PRECISION Y0(PLEN),DY0(PLEN),V(PLEN,PLEN),PCOR(4,NPARTICLES
     &     )
C
C BY REPLACING THE INCLUDE FILE "abcfit_aibi_evol.inc"
C WITH A ".inc" FILE MADE BY MAKE_PARAM YOU CAN OVERRIDE THE
C RUN-TIME READIN OF A DAT-FILE. _BUT_ YOU CAN THEN _ONLY_
C HAVE ONE PARAMETERISATION PER SUBROUTINE.....
C
      INTEGER NPAR_FILE,NJET_FILE,YDEF,COVDEF,PARID(NPARTICLES)
      INTEGER NFUNC(5,NUMFUNC,MAXNUMOFBLOCKS),PARIDBLOCK(NPARTICLES
     &     ,MAXNUMOFBLOCKS),IPARAM(NPARTICLES)
      DOUBLE PRECISION PBOUND(2,NUMFUNC,MAXNUMOFBLOCKS),CBOUND(2,NUMFUNC
     &  ,MAXNUMOFBLOCKS)
      DOUBLE PRECISION PFUNC(NFPAR,NUMFUNC,MAXNUMOFBLOCKS)
      INTEGER NBIN1(MAXNUMOFBLOCKS),NBIN2(MAXNUMOFBLOCKS)
     &  ,PNUMREF(MAXNUMOFBLOCKS)
      DOUBLE PRECISION DBIN1(2,MAXNUMOFBLOCKS),DBIN2(2,MAXNUMOFBLOCKS)
      DOUBLE PRECISION COR(NUMPAR,NUMPAR,NMATR,MAXNUMOFBLOCKS)

      LOGICAL USED
      INTEGER ITEVOLKMN_OLD(3,MAXNUMOFBLOCKS),UNT,I,ITEMP(NPARTICLES)
     &  ,NFBLOCK(MAXNUMOFBLOCKS),NUMCOV(MAXNUMOFBLOCKS),J,H,K,L
     &     ,IBLOCK,PARID_OLD(MAXNUMOFBLOCKS),NPARAMS,LTAG
     &     ,PARIDFILE(NPARTICLES,MAXNUMOFBLOCKS),ISTAT,FLEN
     &     ,NJET_OLD(MAXNUMOFBLOCKS)

      INTEGER INIT1,LENOCC
      PARAMETER (INIT1=3*MAXNUMOFBLOCKS)

      DOUBLE PRECISION WORK(PLEN),P,TH,PARAM,FUNC_SYMMETRIC
     $     ,FUNC_ASYMMETRIC,FUNC_XDFUNC,VAL_DATA,VAL_MC

      CHARACTER*10 NUMBERS
      CHARACTER*20 ABCTAG_OLD(MAXNUMOFBLOCKS)
      CHARACTER*50 FILENAME
      CHARACTER*72 LINE

      EXTERNAL FUNC_SYMMETRIC,FUNC_ASYMMETRIC,FUNC_XDFUNC

      DATA LTAG /0/
      DATA ITEVOLKMN_OLD /INIT1*-2/, PARID_OLD /MAXNUMOFBLOCKS*-2/,
     &     NJET_OLD /MAXNUMOFBLOCKS*-2/
      DATA NUMBERS /'0123456789'/

      SAVE

      V(1,1)=-1.0D0
      IERR=1
      H=0
      DO I=1,NJET
        PARID(I)=MOD(PARIDIN(I),100)
        H=H+PARID(I)*100**(I-1)
        IPARAM(I)=NPAR(MOD(PARIDIN(I),10))
        IF (MOD(PARIDIN(I)/10,10).GT.1) IPARAM(I)=IPARAM(I)+1
      ENDDO
      IF (NBLOCK.GT.0) THEN
        USED=.TRUE.
        DO IBLOCK=1,NBLOCK
          IF (ITEVOLKMN(1).NE.ITEVOLKMN_OLD(1,IBLOCK).OR.ITEVOLKMN(2).NE
     &      .ITEVOLKMN_OLD(2,IBLOCK).OR.ITEVOLKMN(3).NE.ITEVOLKMN_OLD(3
     &         ,IBLOCK).OR.H.NE.PARID_OLD(IBLOCK).OR.NJET.NE
     &         .NJET_OLD(IBLOCK).OR.ABCTAG.NE.ABCTAG_OLD(IBLOCK)) THEN
            USED=.FALSE.
          ELSE
            USED=.TRUE.
            GOTO 100
          ENDIF
        ENDDO
        IBLOCK=IBLOCK-1
 100    CONTINUE
      ELSE
        USED=.FALSE.
        IBLOCK=1
      ENDIF
C
C ASKING FOR NEW PARAMETERISATION???
C

      IF (.NOT.USED) THEN
C
C YES -- MAKE FILENAME AND PERFORM READIN
C
        ISTAT = 1
        NBLOCK=MIN(MAXNUMOFBLOCKS,NBLOCK+1)
        IBLOCK=NBLOCK
        PARID_OLD(IBLOCK)=H
        ITEVOLKMN_OLD(1,IBLOCK)=ITEVOLKMN(1)
        ITEVOLKMN_OLD(2,IBLOCK)=ITEVOLKMN(2)
        ITEVOLKMN_OLD(3,IBLOCK)=ITEVOLKMN(3)
        ABCTAG_OLD(IBLOCK)=ABCTAG
        NJET_OLD(IBLOCK)=NJET
        FLEN=18+NJET
        
        DO WHILE (ISTAT.GT.0.AND.ISTAT.LE.NJET)
          DO I=1,NJET
            PARIDFILE(I,IBLOCK)=PARID(I)
          ENDDO
          LINE=' '
          H=1
          DO I=1,NJET
            IF (PARIDFILE(I,IBLOCK).GT.19) THEN
              LINE=LINE(1:H)//NUMBERS(MOD(PARIDFILE(I,IBLOCK)/10,10)
     &             +1:MOD(PARIDFILE(I,IBLOCK)/10,10)+1)
              H=H+1
            ENDIF
            LINE=LINE(1:H)//NUMBERS(MOD(PARIDFILE(I,IBLOCK),10)
     &           +1:MOD(PARIDFILE(I,IBLOCK),10)+1)
            H=H+1
          ENDDO
          LTAG=MAX(0,INDEX(ABCTAG,' ')-1)
          IF (LTAG.GT.0) THEN
            FILENAME='aibi_evol_'//ABCTAG(1:LTAG)//'_'//LINE(2:H)/
     &           /'_'//NUMBERS(ITEVOLKMN(3)+1:ITEVOLKMN(3)+1)//'0'/
     &           /NUMBERS(ITEVOLKMN(1)+1:ITEVOLKMN(1)+1)//'.dat'
            LTAG=LTAG+1
          ELSE
            FILENAME='aibi_evol_'//LINE(2:H)//'_'/
     &           /NUMBERS(ITEVOLKMN(3)+1:ITEVOLKMN(3)+1)//'0'/
     &           /NUMBERS(ITEVOLKMN(1)+1:ITEVOLKMN(1)+1)//'.dat'
          ENDIF
          FLEN=17+H+LTAG
C
C GET UNIT TO USE
C
          UNT=90
          USED=.TRUE.
          DO WHILE (USED)
            INQUIRE(UNIT=UNT,OPENED=USED)
            UNT=UNT+1
          ENDDO
          UNT=UNT-1
C
C CHECK FILE EXISTENCE
C
          INQUIRE(FILE=FILENAME(1:18+H-1+LTAG),EXIST=USED)
          IF (.NOT.USED) THEN
            WRITE(6,*) 'AIBI_EVOL ERROR: Evolution file ',FILENAME(1:18
     &           +H-1+LTAG),' NOT found!!!!'
            NBLOCK=NBLOCK-1
            GOTO 999
          ENDIF
          ISTAT=0
        ENDDO
C
        OPEN(UNIT=UNT,FILE=FILENAME(1:18+H-1+LTAG),FORM='FORMATTED',
     &    STATUS='OLD')
C
C FIRST READ COMMENT LINES
C       
       

C...JUST A HACK
        DO I=1,20
          READ(UNT,'(a)',ERR=900) LINE
c          WRITE(6,*) LINE
        ENDDO
C
C READ SOME INFO ON THE FILE...
C
        READ(UNT,'(a)',ERR=900) LINE
        READ(UNT,*,ERR=900) NJET_FILE
        READ(UNT,'(a)',ERR=900) LINE
        READ(UNT,*,ERR=900) ITEMP
        NPAR_FILE=0
        DO I=1,NJET_FILE
          IF (ITEMP(I).GE.0) NPAR_FILE=NPAR_FILE+1
        ENDDO

        IF (NJET_FILE.NE.NJET) THEN
          WRITE(6,*) 'AIBI_EVOL ERROR: NJET inconsistency: ',
     $          NJET_FILE,NJET
          NBLOCK=NBLOCK-1
          GOTO 999
        ENDIF
C
C GET READIN VERSION
C
        READ(UNT,'(a)',ERR=900) LINE
        READ(UNT,*,ERR=900) YDEF,COVDEF
C
C NOW LOOP AND READ FUNCTION DEFINITIONS AND PARAMETERS FOR Y0 AND DY0
C
        IF (YDEF.LE.1) THEN
          NFBLOCK(IBLOCK)=0
          DO I=1,NPAR_FILE
C
C READ ENERGY AND COS_THETA RANGES FOR THIS PARAMETERISATION BLOCK
C
            READ(UNT,'(a)',ERR=900) LINE
            READ(UNT,*,ERR=900) ITEMP(I),(WORK(H),H=1,4)
            L=NPAR(MOD(ITEMP(I),10))
            IF (MOD(ITEMP(I)/10,10).GT.1) L=L+1
C
C NOW READ EACH FUNCTION-BLOCK
C
            DO J=1,2*L
              NFBLOCK(IBLOCK)=NFBLOCK(IBLOCK)+1
              IF (NFBLOCK(IBLOCK).GT.NUMFUNC) THEN
                WRITE(6,*)
     $            'AIBI_EVOL ERROR: Too many function blocks'
     &            //'given: ',NFBLOCK(IBLOCK)
                NBLOCK=NBLOCK-1
                GOTO 999
              ENDIF
              DO H=1,2
                PBOUND(H,NFBLOCK(IBLOCK),IBLOCK)=WORK(H)
                CBOUND(H,NFBLOCK(IBLOCK),IBLOCK)=WORK(2+H)
              ENDDO
C
C NFUNC(1,J)=FUNCTYP,NFUNC(2,J)=PARNUM,NFUNC(3,J)=STARTJET,
C NFUNC(4,J)=ENDJET,NFUNC(5,J)=#PARAMS
C
              READ(UNT,'(a)',ERR=900) LINE
              READ(UNT,*,ERR=900) (NFUNC(H,NFBLOCK(IBLOCK),IBLOCK),H=1,5
     &          )
              IF (NFUNC(1,NFBLOCK(IBLOCK),IBLOCK).GE.4) THEN
                WRITE(6,*)
     $            'AIBI_EVOL ERROR: Unknown parameterisation '
     &            //'function referenced: ',NFUNC(1,NFBLOCK(IBLOCK)
     &            ,IBLOCK)
                NBLOCK=NBLOCK-1
                GOTO 999
              ENDIF
              IF (NFUNC(1,NFBLOCK(IBLOCK),IBLOCK).LE.0) GOTO 200
              IF (NFUNC(5,NFBLOCK(IBLOCK),IBLOCK).GT.NFPAR) THEN
                WRITE(6,*)
     $            'AIBI_EVOL ERROR: Too many functions parameters'
     &            //'given: ',NFUNC(5,NFBLOCK(IBLOCK),IBLOCK)
                NBLOCK=NBLOCK-1
                GOTO 999
              ENDIF
              READ(UNT,'(a)',ERR=900) LINE
              DO H=1,NFUNC(5,NFBLOCK(IBLOCK),IBLOCK),3
                READ(UNT,*,ERR=900) (PFUNC(K,NFBLOCK(IBLOCK),IBLOCK),K=H
     &            ,MIN(H+2,NFUNC(5,NFBLOCK(IBLOCK),IBLOCK)))
              ENDDO
            ENDDO
          ENDDO
 200      CONTINUE
        ELSE
          WRITE(6,*)
     $          'AIBI_EVOL ERROR: Unknown read format for Y0 and DY0'
     &          //' requested: ',YDEF
          NBLOCK=NBLOCK-1
          GOTO 999
        ENDIF
C
C NOW READ DEFINITIONS AND PARAMETERS FOR CORRELATION MATRICES
C
        NUMCOV(IBLOCK)=0
        IF (ITEVOLKMN(2).EQ.1) THEN
          IF (COVDEF.LE.1) THEN
            READ(UNT,'(a)',ERR=900) LINE
            READ(UNT,*,ERR=900) PNUMREF(IBLOCK),NBIN1(IBLOCK)
     &        ,NBIN2(IBLOCK)
            IF (NBIN1(IBLOCK)*NBIN2(IBLOCK).GT.NMATR) THEN
              WRITE(6,*)
     $          'AIBI_EVOL ERROR: Too many correlation'
     &          //'-matrices requested: ',NBIN1(IBLOCK)*NBIN2(IBLOCK)
              NBLOCK=NBLOCK-1
              GOTO 999
            ENDIF
            READ(UNT,'(a)',ERR=900) LINE
            READ(UNT,*,ERR=900) (DBIN1(J,IBLOCK),J=1,2),(DBIN2(J,IBLOCK)
     &        ,J=1,2)
            NPARAMS=0
            DO I=1,NJET
              NPARAMS=NPARAMS+IPARAM(I)
            ENDDO
            DO I=1,NBIN1(IBLOCK)*NBIN2(IBLOCK)
              READ(UNT,'(a)',ERR=900) LINE
              READ(UNT,*,ERR=900) ((COR(J,H,I,IBLOCK),J=1,NPARAMS)
     &          ,H=1,NPARAMS)
            ENDDO
            NUMCOV(IBLOCK)=NBIN1(IBLOCK)*NBIN2(IBLOCK)
          ELSE
            WRITE(6,*)
     $          'AIBI_EVOL ERROR: Unknown read format for correlation'
     &          //'-matrix requested: ',COVDEF
            NBLOCK=NBLOCK-1
            GOTO 999
          ENDIF
        ENDIF
        CLOSE(UNT)
      ENDIF
C
C NOW CONTRUCT Y0 AND DY0
C
      NPARAMS=0
      DO I=1,NJET-NUP
        DO J=1,NFBLOCK(IBLOCK)
          IF (I.GE.NFUNC(3,J,IBLOCK).AND.I.LE.NFUNC(4,J,IBLOCK)) THEN
            P = SQRT(DBLE(P_COR(1,I)*P_COR(1,I)+P_COR(2,I)*P_COR(2,I)
     &        +P_COR(3,I)*P_COR(3,I)))
            TH = 1.0D0
            IF (P.GT.0.0D0) TH = DBLE(P_COR(3,I))/P

COB NEW VARIABLES E*COSTH FOR BABAR
            
            IF (SMOOTHVARFLAG(1:LENOCC(SMOOTHVARFLAG)).EQ.'PT') THEN
               TH=DBLE(P_COR(3,I))
            ENDIF
            IF (SMOOTHVARFLAG(1:LENOCC(SMOOTHVARFLAG)).EQ.'P') THEN
               TH=P
            ENDIF
            IF (SMOOTHVARFLAG(1:LENOCC(SMOOTHVARFLAG)).EQ.'E') THEN
               TH=DBLE(P_COR(4,I))
            ENDIF
            IF (SMOOTHVARFLAG(1:LENOCC(SMOOTHVARFLAG)).EQ.
     &           'PvsABS(TH)') THEN
               TH=ABS(TH)
            ENDIF
             IF (SMOOTHVARFLAG(1:LENOCC(SMOOTHVARFLAG)).EQ.'MMISS') THEN
               TH=DBLE(MMISS)
            ENDIF
            
            IF (SMOOTHVARFLAG.EQ.' ') TH=SIGN(MIN(1.0D0,ABS(TH)),TH)

            TH=MIN(MAX(TH,CBOUND(1,J,IBLOCK)),CBOUND(2,J,IBLOCK))
            IF (TH.EQ.CBOUND(1,J,IBLOCK)) IERR=3
            IF (TH.EQ.CBOUND(2,J,IBLOCK)) IERR=4
            P = MIN(MAX(DBLE(P_COR(4,I)),PBOUND(1,J,IBLOCK)),PBOUND(2,J
     &        ,IBLOCK))

            IF (SMOOTHVARFLAG.EQ.'PvsABS(TH)') THEN
               P = SQRT(DBLE(P_COR(1,I)*P_COR(1,I)+P_COR(2,I)*P_COR(2,I)
     &        +P_COR(3,I)*P_COR(3,I)))
            ENDIF

            IF (P.EQ.PBOUND(1,J,IBLOCK)) IERR=3
            IF (P.EQ.PBOUND(2,J,IBLOCK)) IERR=4
            PARAM=0.0D0
            IF (NFUNC(1,J,IBLOCK).EQ.2) PARAM=FUNC_SYMMETRIC(TH,P
     &        ,PFUNC(1,J,IBLOCK))
            IF (NFUNC(1,J,IBLOCK).EQ.1) PARAM=FUNC_ASYMMETRIC(TH,P
     &        ,PFUNC(1,J,IBLOCK))

            IF (NFUNC(1,J,IBLOCK).EQ.3) PARAM=FUNC_XDFUNC(TH,
     &           PFUNC(1,J,IBLOCK),1)
            IF (NFUNC(1,J,IBLOCK).EQ.3) VAL_MC=FUNC_XDFUNC(TH,
     &           PFUNC(1,J,IBLOCK),2)
            IF (NFUNC(1,J,IBLOCK).EQ.3) VAL_DATA=FUNC_XDFUNC(TH,
     &           PFUNC(1,J,IBLOCK),3)
        
            IF (NFUNC(2,J,IBLOCK).LE.IPARAM(I)) THEN
C
C Y0
C            
              Y0(NPARAMS+NFUNC(2,J,IBLOCK))=PARAM
            ELSE
C     
C     DY0
C
COB DATA ONLY
c               PARAM=VAL_DATA/VAL_MC*PARAM
              DY0(NPARAMS+NFUNC(2,J,IBLOCK)-IPARAM(I))=PARAM

            ENDIF
          ENDIF
        ENDDO
        IF (((PARIDIN(I).GE.100.AND.PARIDIN(I).LT.200).OR.(PARIDIN(I).GE
     &       .300.AND.PARIDIN(I).LT.400)).AND.MOD(PARIDIN(I)/10,10).GT.0
     &       ) THEN
          Y0(NPARAMS+2)=Y0(NPARAMS+4)
          DY0(NPARAMS+2)=DY0(NPARAMS+4)
        ENDIF
        NPARAMS=NPARAMS+NPARAM(I)
      ENDDO

      IF (NUP.GT.0) THEN
        NPARAMS=0
        H=0
        DO I=1,NJET
          ITEMP(I)=PARIDIN(I)
          IF (I.LE.NJET-NUP) NPARAMS=NPARAMS+NPARAM(I)
          H=H+NPARAM(I)
          DO J=1,4
            PCOR(J,I)=DBLE(P_COR(J,I))
          ENDDO
        ENDDO
        CALL CVP2Y(NJET,ITEMP,PCOR,WORK,I)
        DO I=NPARAMS+1,H
          Y0(I)=WORK(I)
          DY0(I)=1.0D0
        ENDDO
      ENDIF

      NPARAMS=0
      DO I=1,NJET
        NPARAMS=NPARAMS+NPARAM(I)
      ENDDO
      IF (NUMCOV(IBLOCK).LE.0) THEN
        DO I=1,NPARAMS  
          V(I,I) = DY0(I)*DY0(I)
        ENDDO
        IF (ITEVOLKMN(2).EQ.1) THEN
          WRITE(6,*) 'AIBI_EVOL WARNING:'
     &       //' Only diagonal covariance matrix available'
          GOTO 999
        ENDIF
      ENDIF
C
C NOW DETERMINE CORRELATION MATRIX
C
      IF (NUMCOV(IBLOCK).GT.0.AND.ITEVOLKMN(2).EQ.1) THEN
        P = SQRT(DBLE(P_COR(1,PNUMREF(IBLOCK))*P_COR(1,PNUMREF(IBLOCK))
     &    +P_COR(2,PNUMREF(IBLOCK))*P_COR(2,PNUMREF(IBLOCK))+P_COR(3
     &    ,PNUMREF(IBLOCK))*P_COR(3,PNUMREF(IBLOCK))))
        TH = 1.0D0
        IF (P.GT.0.0D0) TH = DBLE(P_COR(3,PNUMREF(IBLOCK)))/P
        TH=SIGN(MIN(1.0D0,ABS(TH)),TH)
        CALL GETCORBIN(NBIN1(IBLOCK),NBIN2(IBLOCK),DBIN1(1,IBLOCK)
     &    ,DBIN2(1,IBLOCK),P,TH,H)
        IF (H.LE.0) THEN
          WRITE(6,*) 'AIBI_EVOL ERROR: No correlation bin found'
          GOTO 999
        ENDIF
C     
C FINALLY BUILD FULL COVARIANCE MATRIX
C
        DO I=1,NPARAMS
          DO J=1,NPARAMS
            V(I,J)=COR(I,J,H,IBLOCK)*DY0(I)*DY0(J)
          ENDDO
        ENDDO
        WRITE(6,*)
     &       'AIBI_EVOL WARNING -- Correlation matrix setup is '
     &       //'not completely implemented'
        WRITE(6,*)
     &       '                     Longitudinal params only is missing'

      ENDIF
      RETURN
 900  CONTINUE
      WRITE(6,*) 'AIBI_EVOL ERROR: Error reading from file '
     $     ,FILENAME(1:FLEN)
 999  RETURN
      END
C
C DUMMY ROUTINE
C
      SUBROUTINE AIBI_EVOL_USER(PARID,ITEVOLKMN,NJET,NUP,P_COR
     $     ,Y0,DY0,V,IERR)
*****************************************************
* USER Routine to return covariance matrix and starting values
* of the fit parameters depending on choice of parameterisation
*****************************************************
      IMPLICIT NONE

#include "abcfit_setup.inc"

      INTEGER PARID(NPARTICLES),ITEVOLKMN(3),NJET,NUP,IERR
      REAL P_COR(4,NPARTICLES)
      DOUBLE PRECISION Y0(PLEN),DY0(PLEN),V(PLEN,PLEN)
C
C INTERNAL VARIABLES
C
      INTEGER I,J,ITYPPN,ITEVOL,ITETRUE
      DOUBLE PRECISION WORK(PLEN,2)
C
C DEFAULT PRINTOUT (UNCOMMENT IF NOTHING ELSE IS DONE)
C
      WRITE(6,*) 'AIBI_EVOL_USER ERROR: Dummy version of AIBI_EVOL_USER'
     &     //' called'
      V(1,1)=-1.0D0
      IERR=2
C
 999  RETURN
      END
C**************************************************************************
C
C GENERAL SUPPORT ROUTINES FOR AIBI EVOLUTION FUNCTIONS..
C
C**************************************************************************
      SUBROUTINE GETCORBIN(NBIN1,NBIN2,DBIN1,DBIN2,P,TH,IBIN)
      IMPLICIT NONE

      INTEGER NBIN1,NBIN2,IBIN
      DOUBLE PRECISION P,TH
      DOUBLE PRECISION DBIN1(2),DBIN2(2)

      REAL RTEMP

      IBIN=0

      IF (NBIN1.GT.0) THEN
        RTEMP = (DBIN1(2)-DBIN1(1))/REAL(NBIN1)
        IF (P-DBIN1(1).GE.0.0.AND.P-DBIN1(2)
     &    .LE.0.0) IBIN = INT((P-DBIN1(1))/RTEMP)+1
      ENDIF

      IF (NBIN2.GT.0) THEN
        RTEMP = (DBIN2(2)-DBIN2(1))/REAL(NBIN2)
        IF (TH-DBIN2(1).GE.0.0.AND.TH-DBIN2(2)
     &    .LE.0.0) THEN
          IF (IBIN.GT.0) THEN
            IBIN =IBIN+INT((TH-DBIN2(1))/RTEMP)*NBIN1
          ELSE
            IF (NBIN1.EQ.0) IBIN = INT((TH-DBIN2(1))
     &        /RTEMP)+1
          ENDIF
        ELSE
          IBIN=0
        ENDIF
      ENDIF

 999  RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FUN_E_DEPENDENCE (E,PAR)
C-----------------------------------------------------------------------
      IMPLICIT NONE

      DOUBLE PRECISION PAR(6),E
      FUN_E_DEPENDENCE = PAR(1)+PAR(2)*E+PAR(3)*E**2+PAR(4)*E**3+
     +     PAR(5)*E**4+PAR(6)*E**5
      END      
      
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FUNC_XDFUNC (X,PAR,ITAG)
C-----------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8 X
      REAL*8 PAR(18)
      REAL*8 PAR0(6)
      EXTERNAL FUN_E_DEPENDENCE
      REAL*8 FUN_E_DEPENDENCE
      INTEGER I,ITAG,ISTART
      
      IF (ITAG.EQ.1) ISTART=0
      IF (ITAG.EQ.2) ISTART=6
      IF (ITAG.EQ.3) ISTART=12
      DO I = 1,6
         PAR0(I) = PAR(I+ISTART)
      ENDDO
      
      FUNC_XDFUNC = FUN_E_DEPENDENCE(X,PAR0)

      END


C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FUNC_SYMMETRIC (TH,E,PAR)
C-----------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8 TH,E
      REAL*8 PAR(18)
      REAL*8 PAR0(6),PAR6(6),PAR8(6),TH_0,TH_6,TH_8
      EXTERNAL FUN_E_DEPENDENCE
      REAL*8 FUN_E_DEPENDENCE
      INTEGER I
      
      DO I = 1,6
         PAR0(I) = PAR(I)
         PAR6(I) = PAR(I+6)
         PAR8(I) = PAR(I+12)
      ENDDO
         
      TH_0 = FUN_E_DEPENDENCE(E,PAR0)
      TH_6 = FUN_E_DEPENDENCE(E,PAR6)
      TH_8 = FUN_E_DEPENDENCE(E,PAR8)
      FUNC_SYMMETRIC = TH_0+TH_6*TH**6+TH_8*TH**8
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION FUNC_ASYMMETRIC (TH,E,PAR)
C-----------------------------------------------------------------------
      IMPLICIT NONE

      REAL*8 TH,E
      REAL*8 PAR(18)
      REAL*8 PAR0(6),PAR3(6),PAR5(6),TH_0,TH_3,TH_5
      EXTERNAL FUN_E_DEPENDENCE
      REAL*8 FUN_E_DEPENDENCE
      INTEGER I
      
      DO I = 1,6
         PAR0(I) = PAR(I)
         PAR3(I) = PAR(I+6)
         PAR5(I) = PAR(I+12)
      ENDDO
         
      TH_0 = FUN_E_DEPENDENCE(E,PAR0)
      TH_3 = FUN_E_DEPENDENCE(E,PAR3)
      TH_5 = FUN_E_DEPENDENCE(E,PAR5)
      FUNC_ASYMMETRIC = TH_0+TH_3*TH**3+TH_5*TH**5

      END   











