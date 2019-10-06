#define ABC_BMATRIX_DBG 0
#define ABC_CVY2P_DBG 0
#define ABC_CVP2Y_DBG 0
C============================================================
C B-MATRIX CALCULATION ROUTINES
C============================================================
      SUBROUTINE BMATRIX(IMODE,NJET,IPAR,Y0)
      IMPLICIT NONE

#include "abcfit_setup.inc"
#include "abcfit_bmatrix.inc"

      INTEGER IMODE,NJET,IPAR(NPARTICLES)
      DOUBLE PRECISION Y0(PLEN)

      INTEGER I,I1,J,MODE,NPARAMS
      DOUBLE PRECISION ROW(4*NPARTICLES),WORK(4*NPARTICLES+3,4
     &     *NPARTICLES+3),WORK1(NBFMX,4*NPARTICLES+3),P(4,NPARTICLES)
      DOUBLE PRECISION BETA,DBETADX

      EXTERNAL BETA,DBETADX
C
      MODE=IMODE
      IF (MODE.GT.10) MODE=MODE-10
C
C SET UP FOR CONSTRAINTS AND DERIVATIVES
C
      NPARAMS=0
      DO I=1,NJET
        NPARAMS=NPARAMS+NPARAM(I)
      ENDDO
      ALPHA_IJ = ALPHAIJ
      DALPHA_IJ = 1.0D0
      ALPHA_KL = ALPHAKL
      DALPHA_KL = 1.0D0
      IF (IMODE.GT.10) THEN
        ALPHA_IJ = BETA(1,ALPHAIJ)
        DALPHA_IJ = DBETADX(1,ALPHAIJ)
        ALPHA_KL = 0.0D0
        DALPHA_KL = 0.0D0
        IF (MODE.EQ.5.OR.MODE.EQ.7) THEN
          ALPHA_KL = BETA(2,ALPHAKL)
          DALPHA_KL = DBETADX(2,ALPHAKL)
        ENDIF
      ENDIF
C
      CALL CVY2P(NJET,IPAR,Y0,P)
C
C CALCULATE CONSTRAINTS
C
      CALL FY(MODE,NJET,P)
C 
C     BUILD B-MATRIX USING THE "CHAIN-RULE" OF DERIVATIVE
C
      DO I=1,NJET
        PIR(I)=P(1,I)**2+P(2,I)**2+P(3,I)**2
        EIR(I)=P(4,I)
        MIR(I)=SQRT(MAX(0.0D0,P(4,I)**2-PIR(I)))
        PIR(I)=SQRT(PIR(I))
        BIER(I)=DLOG(MAX(1.0D-32,MAX(1.0D0,EIR(I)/PIR(I))-1.0D0))
      ENDDO
C
      CALL DP_DPARAM(NJET,IPAR,P,WORK)
      WORK(4*NJET+1,NPARAMS+1)=1.0D0 ! ALPHAIJ TERM
      WORK(4*NJET+2,NPARAMS+2)=1.0D0 ! ALPHAKL TERM
C
      J=0                       ! NUMBER OF ALPHA-PARAMETERS
      IF (MODE.EQ.3) J=1
      IF (MODE.EQ.5.OR.MODE.EQ.7) J=2
#if ABC_BMATRIX_DBG > 0
      PRINT *,'ABC_BMATRIX_DEBUG: Dp_Dparam-matrix'
      DO I=1,4*NJET+J
        PRINT '(20f10.4)',(WORK(I,I1),I1=1,NPARAMS)
      ENDDO
#endif
      CALL DF_DP(MODE,NJET,P,WORK1)
#if ABC_BMATRIX_DBG > 0
      PRINT *,'ABC_BMATRIX_DEBUG: Df_Dp-matrix'
      DO I=1,NCONS
        PRINT '(20f12.5)',(WORK1(I,I1),I1=1,4*NJET+J)
      ENDDO
#endif
      CALL DMMLT(NCONS,4*NJET+J,NPARAMS+J,WORK1(1,1),WORK1(1,2),WORK1(2
     &  ,1),WORK(1,1),WORK(1,2),WORK(2,1),B(1,1),B(1,2),B(2,1),ROW)
#if ABC_BMATRIX_DBG > 0
      PRINT *,'ABC_BMATRIX_DEBUG: Df_Dparam-matrix'
      DO I=1,NCONS
        PRINT '(20f12.5)',(B(I,I1),I1=1,NPARAMS)
      ENDDO
#endif
 999  RETURN
      END

      SUBROUTINE DEFUVEC(IPARAM,NJET,P_REC,IERR)
      IMPLICIT NONE

#include "abcfit_setup.inc"
#include "abcfit_bmatrix.inc"

      INTEGER IPARAM(NPARTICLES),NJET,IERR
      REAL P_REC(4,NPARTICLES)

      INTEGER I,J
      REAL RPMD(6,NPARTICLES)
      DOUBLE PRECISION SCALB,DSQRT,PMD(6,NPARTICLES)

      IERR=0

      DO I = 1,NJET
        IF (IPARAM(I).GE.0) THEN
          IF ( P_REC(4,I).EQ.0.0 ) THEN
            IERR = -5
            GOTO 999
          ENDIF
C
C INPUT IS IN REAL -- WE CALCULATE MASS IN REAL _AND THEN_ CONVERT
C
          RPMD(1,I) = P_REC(1,I)
          RPMD(2,I) = P_REC(2,I)
          RPMD(3,I) = P_REC(3,I)
          RPMD(4,I) = P_REC(4,I)
          RPMD(6,I) = RPMD(1,I)**2+RPMD(2,I)**2+RPMD(3,I)**2
          RPMD(5,I) = SQRT(MAX(0.0,RPMD(4,I)**2-RPMD(6,I)))
          RPMD(6,I) = SQRT(RPMD(6,I))
          DO J=1,6
            PMD(J,I) = DBLE(RPMD(J,I))
          ENDDO
          PIM(I)  = DBLE(RPMD(6,I))
          PTIM(I) = DBLE(SQRT(RPMD(1,I)**2+RPMD(2,I)**2))
          EIM(I)  = DBLE(RPMD(4,I))
          MIM(I)  = DBLE(RPMD(5,I))
          BIM(I)  = DLOG(MAX(1.0D-32,MAX(1.0D0,DBLE(RPMD(4,I)/RPMD(6,I))
     &         )-1.0D0))
        ENDIF
      ENDDO
c
c then Construct reference frame (if needed) for BMATRIX routine
c
      
      DO I = 1,NJET
        IF (MOD(IPARAM(I),10).GE.0.AND.MOD(IPARAM(I),10).LE.1.OR
     &       .IPARAM(I).EQ.4) THEN ! JUST A HACK
c r vector
          URX(I) = PMD(1,I)/PMD(6,I)
          URY(I) = PMD(2,I)/PMD(6,I)
          URZ(I) = PMD(3,I)/PMD(6,I)
c phi vector = z X ur
          IF (MOD(IPARAM(I)/100,10).EQ.1.OR.MOD(IPARAM(I)/100,10).EQ.3)
     &         THEN
            UTHX(I) = 0.D0
            UTHY(I) = 0.D0
            UTHZ(I) = 0.D0
            UPHX(I) = 0.D0
            UPHY(I) = 0.D0
            UPHZ(I) = 0.D0
          ELSE
            IF ( PMD(6,I).EQ.PMD(3,I) ) THEN
              UPHX(I) = 1.D0
              UPHY(I) = 0.D0
              UPHZ(I) = 0.D0
            ELSE
              UPHX(I) = -URY(I)
              UPHY(I) = URX(I)
              UPHZ(I) = 0.D0
c Normalize it
              SCALB = DSQRT(UPHX(I)**2+UPHY(I)**2+UPHZ(I)**2)
              IF ( SCALB.EQ.0D0 ) SCALB = 1.D-8
              UPHX(I) = UPHX(I)/SCALB
              UPHY(I) = UPHY(I)/SCALB
              UPHZ(I) = UPHZ(I)/SCALB
            ENDIF
c theta vector = ur X uphi
            UTHX(I) = -(URY(I)*UPHZ(I)-URZ(I)*UPHY(I))
            UTHY(I) = -(URZ(I)*UPHX(I)-URX(I)*UPHZ(I))
            UTHZ(I) = -(URX(I)*UPHY(I)-URY(I)*UPHX(I))
            SCALB = DSQRT(UTHX(I)**2+UTHY(I)**2+UTHZ(I)**2)
            IF ( SCALB.EQ.0D0 ) SCALB = 1.D-8
            UTHX(I) = UTHX(I)/SCALB
            UTHY(I) = UTHY(I)/SCALB
            UTHZ(I) = UTHZ(I)/SCALB
          ENDIF
        ENDIF
      ENDDO
 999  RETURN
      END
C#################################################
C PARAMETER-TRANSLATION INTO MOMENTA (AND VICE VERSE) AND DERIVATIVES
C#################################################
      SUBROUTINE CVY2P(NJET,IPAR,Y,P)
      IMPLICIT NONE
C
#include "abcfit_setup.inc"
#include "abcfit_bmatrix.inc"

      INTEGER NJET,IPAR(NPARTICLES)
      DOUBLE PRECISION P(4,NPARTICLES),Y(PLEN),TEMP(2)

      INTEGER I,J,K
C
      J=1
      DO I = 1, NJET
        K=MOD(IPAR(I),10)
        IF (I.GT.1) J=J+NPARAM(I-1)
        IF (K.EQ.0) THEN     ! ALEPH-PARAMETRISATION
          P(1,I)=Y(J)*URX(I)*PIM(I)+Y(J+1)*UTHX(I)+Y(J+2)*UPHX(I
     &      )                                           
          P(2,I)=Y(J)*URY(I)*PIM(I)+Y(J+1)*UTHY(I)+Y(J+2)*UPHY(I
     &      )                                           
          P(3,I)=Y(J)*URZ(I)*PIM(I)+Y(J+1)*UTHZ(I)+Y(J+2)*UPHZ(I
     &      )
        ELSEIF (K.EQ.4) THEN ! JUST A HACK P_FIT=a*P_RECO
          P(1,I)=Y(J)*URX(I)*PIM(I)
          P(2,I)=Y(J)*URY(I)*PIM(I)                     
          P(3,I)=Y(J)*URZ(I)*PIM(I)
        ELSEIF (K.EQ.1) THEN ! NJK(DELPHI)-PARAMETRISATION
          P(1,I)=EXP(Y(J))*URX(I)*PIM(I)+Y(J+1)*UTHX(I)+Y(J+2)
     &      *UPHX(I)                                        
          P(2,I)=EXP(Y(J))*URY(I)*PIM(I)+Y(J+1)*UTHY(I)+Y(J+2)
     &      *UPHY(I)                                        
          P(3,I)=EXP(Y(J))*URZ(I)*PIM(I)+Y(J+1)*UTHZ(I)+Y(J+2)
     &      *UPHZ(I)
        ELSEIF (K.EQ.2) THEN ! MOMENTUM AND DIRECTION
          P(1,I) = Y(J)*SIN(Y(J+1))*COS(Y(J+2))
          P(2,I) = Y(J)*SIN(Y(J+1))*SIN(Y(J+2))
          P(3,I) = Y(J)*COS(Y(J+1))
        ELSEIF (K.EQ.3) THEN ! 3-MOMENTUM
          P(1,I) = Y(J)
          P(2,I) = Y(J+1)
          P(3,I) = Y(J+2)
       ELSEIF (K.EQ.3) THEN     ! PT,COSTHETA,PHI
          P(1,I) = Y(J)*COS(Y(J+2))
          P(2,I) = Y(J)*SIN(Y(J+2))
          P(3,I) = Y(J)*TAN(Y(J+1))
       ENDIF
        TEMP(1) = P(1,I)**2+P(2,I)**2+P(3,I)**2
        K=MOD(IPAR(I)/10,10)
        IF (K.EQ.0) THEN ! SCALED MASS
          P(4,I) = EIM(I)*SQRT(TEMP(1))/PIM(I)
          IF (SQRT(MAX(P(4,I)**2-TEMP(1),0.0D0)).LE.0.0D0) P(4,I)
     &         =SQRT(TEMP(1))
        ELSEIF (K.EQ.1) THEN !FIXED MASS
          P(4,I) = SQRT(TEMP(1)+MIM(I)**2)
        ELSEIF (K.EQ.2) THEN ! FIRST 4TH PARAMETER
          P(4,I) = SQRT(TEMP(1)+Y(J+NPARAM(I)-1)**2*MIM(I)**2)
        ELSEIF (K.EQ.3) THEN ! SECOND 4TH PARAMETER
          P(4,I) = SQRT(TEMP(1)+EXP(2.0D0*Y(J+NPARAM(I)-1))*MIM(I)**2)
        ELSEIF (K.EQ.4) THEN    ! THIRD 4TH PARAMETER
          P(4,I) = (EXP(Y(J+NPARAM(I)-1)+BIM(I))+1.0D0)*SQRT(TEMP(1))
        ENDIF
#if ABC_CVY2P_DBG > 0
        PRINT'(1X,A,I2,A,4F12.5)','ABC_CVY2P_DBG: P(',I,')=',P(1,I),P(2
     &       ,I),P(3,I),P(4,I)
        PRINT'(1X,A,I1,A,4F12.5)','ABC_CVY2P_DBG: Y(J->J+'
     &       ,NPARAM(I)-1,')=',(Y(K),K=J,J
     &       +NPARAM(I)-1)
#endif
      ENDDO
       
  999 RETURN
      END
C
      SUBROUTINE CVP2Y(NJET,IPAR,P,Y,IERR)
      IMPLICIT NONE
C
#include "abcfit_setup.inc"
#include "abcfit_bmatrix.inc"

      INTEGER NJET,IPAR(NPARTICLES),IERR
      DOUBLE PRECISION P(4,NPARTICLES),Y(PLEN),DTEMP

      INTEGER I,J,K
C
      IERR=0

      J=1
      DO I = 1, NJET
        K=MOD(IPAR(I),10)
        IF (I.GT.1) J=J+NPARAM(I-1)
        IF (K.EQ.0) THEN     ! ALEPH-PARAMETRISATION
          Y(J) = (P(1,I)*URX(I)+P(2,I)*URY(I)+P(3,I)*URZ(
     &      I))/PIM(I) 
          IF (Y(J).LT.0.0D0) IERR=1
          Y(J+1) = P(1,I)*UTHX(I)+P(2,I)*UTHY(I)+P(3,I)*UTHZ(I)
          Y(J+2) = P(1,I)*UPHX(I)+P(2,I)*UPHY(I)+P(3,I)*UPHZ(I)
        ELSEIF (K.EQ.4) THEN ! JUST A HACK P_FIT=a*P_RECO
          Y(J) = (P(1,I)*URX(I)+P(2,I)*URY(I)+P(3,I)*URZ(
     &      I))/PIM(I) 
          IF (Y(J).LT.0.0D0) IERR=1
        ELSEIF (K.EQ.1) THEN ! NJK(DELPHI)-PARAMETRISATION
          Y(J) = (P(1,I)*URX(I)+P(2,I)*URY(I)+P(3,I)*URZ(
     &      I))/PIM(I)
          IF (Y(J).LE.0.0D0) THEN
            IERR=1
            Y(J) = -1.0D-50
          ELSE
            Y(J) = LOG(ABS(Y(J)))
          ENDIF
          Y(J+1) = P(1,I)*UTHX(I)+P(2,I)*UTHY(I)+P(3,I)*UTHZ(I)
          Y(J+2) = P(1,I)*UPHX(I)+P(2,I)*UPHY(I)+P(3,I)*UPHZ(I)
        ELSEIF (K.EQ.2) THEN ! MOMENTUM AND DIRECTION
          Y(J) = SQRT(P(1,I)**2+P(2,I)**2+P(3,I)**2)
          IF (Y(J).LE.0.0D0) THEN
            IERR=1
            Y(J) = 1.0D-50
            Y(J+1) = 1.0D0
            Y(J+2) = 0.0D0
          ELSE
            Y(J+1) = ACOS(P(3,I)/Y(J))
            Y(J+2) = ATAN2(P(2,I),P(1,I))
          ENDIF
        ELSEIF (K.EQ.3) THEN ! 3-MOMENTUM
          Y(J) = P(1,I)
          Y(J+1) = P(2,I)
          Y(J+2) = P(3,I)
        ENDIF
        K=MOD(IPAR(I)/10,10)
        IF (K.EQ.2) THEN ! FIRST 4TH PARAMETER
          IF (MIM(I).LE.0.0D0) THEN
            Y(J+NPARAM(I)-1) = 1.0D0
            IERR=1
          ELSE
            Y(J+NPARAM(I)-1) = SQRT(MAX(0.0D0,P(4,I)**2-P(1,I)**2-P(2,I)
     &           **2-P(3,I)**2))/MIM(I)
          ENDIF
        ELSEIF (K.EQ.3) THEN ! SECOND 4TH PARAMETER
          IF (MIM(I).LE.0.0D0) THEN
            Y(J+NPARAM(I)-1) = 1.0D0
            IERR=1
          ELSE
            IF (P(4,I)**2-P(1,I)**2-P(2,I)**2-P(3,I)**2.LE.0.0D0) THEN
              Y(J+NPARAM(I)-1) = 0.0D0
              IERR=1
            ELSE
              Y(J+NPARAM(I)-1) = LOG(SQRT(MAX(1.0D-32,P(4,I)**2-P(1,I)
     &             **2-P(2,I)**2-P(3,I)**2))/MIM(I))
            ENDIF
          ENDIF
        ELSEIF (K.EQ.4) THEN ! THIRD 4TH PARAMETER
          DTEMP=SQRT(P(1,I)**2+P(2,I)**2+P(3,I)**2)
          IF (DTEMP.GT.0.0D0) THEN
            Y(J+NPARAM(I)-1) = LOG(MAX(1.0D-32,P(4,I)/SQRT(P(1,I)**2+P(2
     &           ,I)**2+P(3,I)**2)-1.0D0))-LOG(MAX(1.0D-32,EIM(I)/PIM(I)
     &           -1.0D0))
          ELSE
            Y(J+NPARAM(I)-1) = 1.0D32
            IERR=1
          ENDIF
          IF (P(4,I).LE.DTEMP.OR.EIM(I).LE.PIM(I)) IERR=1
        ENDIF
#if ABC_CVP2Y_DBG > 0
        PRINT'(1X,A,I2,A,4F12.5)','ABC_CVP2Y_DBG: P(',I,')=',P(1,I),P(2
     &       ,I),P(3,I),P(4,I)
        PRINT'(1X,A,I1,A,4F12.5)','ABC_CVP2Y_DBG: Y(J->J+'
     &       ,NPARAM(I)-1,')=',(Y(K),K=J,J
     &       +NPARAM(I)-1)
        PRINT'(1X,A,I1)','ABC_CVP2Y_DBG: IERR=',IERR
#endif
      ENDDO
  999 RETURN
      END
C
      SUBROUTINE DP_DPARAM(NJET,IPAR,P,DPDY)
      IMPLICIT NONE

#include "abcfit_setup.inc"
#include "abcfit_bmatrix.inc"

      INTEGER IPAR(NPARTICLES),NJET
      DOUBLE PRECISION P(4,NPARTICLES),DPDY(4*NPARTICLES+3,4*NPARTICLES
     &     +3)

      INTEGER I,J,H,K
      DOUBLE PRECISION SCALE

      CALL DMSET(4*NJET+3,4*NJET+3,0.0D0,DPDY(1,1),DPDY(1,2),DPDY(2,1))
      J=1
      DO I = 1, NJET
        K=MOD(IPAR(I),10)
        IF (I.GT.1) J=J+NPARAM(I-1)
        H=(I-1)*4
        IF (K.EQ.0) THEN     ! ALEPH-PARAMETRISATION
          DPDY(1+H,J) = PIM(I)*URX(I)
          DPDY(2+H,J) = PIM(I)*URY(I)
          DPDY(3+H,J) = PIM(I)*URZ(I)
          IF (IPAR(I)/100.NE.1.AND.IPAR(I)/100.NE.3) THEN
            DPDY(1+H,J+1) = UTHX(I)
            DPDY(2+H,J+1) = UTHY(I)
            DPDY(3+H,J+1) = UTHZ(I)
            DPDY(1+H,J+2) = UPHX(I)
            DPDY(2+H,J+2) = UPHY(I)
            DPDY(3+H,J+2) = UPHZ(I)
          ENDIF
        ELSEIF (K.EQ.4) THEN  ! JUST A HACK P_FIT=a*P_RECO
          DPDY(1+H,J) = PIM(I)*URX(I)
          DPDY(2+H,J) = PIM(I)*URY(I)
          DPDY(3+H,J) = PIM(I)*URZ(I)
        ELSEIF (K.EQ.1) THEN ! NJK(DELPHI)-PARAMETRISATION
          DPDY(1+H,J) = EXP(AI(I))*PIM(I)*URX(I)
          DPDY(2+H,J) = EXP(AI(I))*PIM(I)*URY(I)
          DPDY(3+H,J) = EXP(AI(I))*PIM(I)*URZ(I)
          IF (IPAR(I)/100.NE.1.AND.IPAR(I)/100.NE.3) THEN
            DPDY(1+H,J+1) = UTHX(I)
            DPDY(2+H,J+1) = UTHY(I)
            DPDY(3+H,J+1) = UTHZ(I)
            DPDY(1+H,J+2) = UPHX(I)
            DPDY(2+H,J+2) = UPHY(I)
            DPDY(3+H,J+2) = UPHZ(I)
          ENDIF
        ELSEIF (K.EQ.2) THEN ! MOMENTUM AND DIRECTION
          DPDY(1+H,J) = SIN(BI(I))*COS(CI(I))
          DPDY(2+H,J) = SIN(BI(I))*SIN(CI(I))
          DPDY(3+H,J) = COS(BI(I))
          IF (IPAR(I)/100.NE.1.AND.IPAR(I)/100.NE.3) THEN
            DPDY(1+H,J+1) = AI(I)*COS(BI(I))*COS(CI(I))
            DPDY(2+H,J+1) = AI(I)*COS(BI(I))*SIN(CI(I))
            DPDY(3+H,J+1) = -AI(I)*SIN(BI(I))
            DPDY(1+H,J+2) = -AI(I)*SIN(BI(I))*SIN(CI(I))
            DPDY(2+H,J+2) =  AI(I)*SIN(BI(I))*COS(CI(I))
            DPDY(3+H,J+2) =  0.0D0
          ENDIF
        ELSEIF (K.EQ.3) THEN ! 3-MOMENTUM
          DPDY(1+H,J) = 1.0D0
          IF (IPAR(I)/100.NE.1.AND.IPAR(I)/100.NE.3) THEN
            DPDY(2+H,J+1) = 1.0D0
            DPDY(3+H,J+2) = 1.0D0
          ENDIF
       ELSEIF (K.EQ.4) THEN     ! PT,COSTHETA,PHI   
          DPDY(1+H,J) = COS(CI(I)) 
          DPDY(2+H,J) = SIN(CI(I))
          DPDY(3+H,J) = 0.0D0
          IF (IPAR(I)/100.NE.1.AND.IPAR(I)/100.NE.3) THEN
            DPDY(1+H,J+1) = 0.0D0
            DPDY(2+H,J+1) = 0.0D0
            DPDY(3+H,J+1) = PIM(I)
            DPDY(1+H,J+2) = -AI(I)*SIN(CI(I))
            DPDY(2+H,J+2) =  AI(I)*COS(CI(I))
            DPDY(3+H,J+2) =  0.0D0
          ENDIF
        ENDIF
        K=MOD(IPAR(I)/10,10)
        IF (K.EQ.1) THEN ! FIXED MASS
          DPDY(4+H,J) = (DPDY(1+H,J)*P(1,I)+DPDY(2+H,J)*P(2,I)
     &         +DPDY(3+H,J)*P(3,I))/EIR(I)
          DPDY(4+H,J+1) =(DPDY(1+H,J+1)*P(1,I)+DPDY(2+H,J+1)*P(2,I
     &         )+DPDY(3+H,J+1)*P(3,I))/EIR(I)
          DPDY(4+H,J+2) =(DPDY(1+H,J+2)*P(1,I)+DPDY(2+H,J+2)*P(2,I
     &         )+DPDY(3+H,J+2)*P(3,I))/EIR(I)
        ELSE
          SCALE=EIM(I)/PIM(I)
          IF (K.EQ.4) SCALE=1.0D0+EXP(DI(I)+BIM(I))
          DPDY(4+H,J) = SCALE*(DPDY(1+H,J)*P(1,I)+DPDY(2+H,J)*P(2,I)
     &         +DPDY(3+H,J)*P(3,I))/PIR(I)
          DPDY(4+H,J+1) =SCALE*(DPDY(1+H,J+1)*P(1,I)+DPDY(2+H,J+1)*P(2,I
     &         )+DPDY(3+H,J+1)*P(3,I))/PIR(I)
          DPDY(4+H,J+2) =SCALE*(DPDY(1+H,J+2)*P(1,I)+DPDY(2+H,J+2)*P(2,I
     &         )+DPDY(3+H,J+2)*P(3,I))/PIR(I)
        ENDIF
        IF (K.EQ.2) THEN ! FIRST 4TH PARAMETER
          DPDY(4+H,J+NPARAM(I)-1) = DI(I)*MIM(I)**2/EIR(I)
        ELSEIF (K.EQ.3) THEN ! SECOND 4TH PARAMETER
          DPDY(4+H,J+NPARAM(I)-1) = EXP(2.0D0*DI(I))*MIM(I)**2/EIR(I)
        ELSEIF (K.EQ.4) THEN ! THIRD 4TH PARAMETER
          DPDY(4+H,J+NPARAM(I)-1) = EXP(DI(I)+BIM(I))*PIR(I)
        ENDIF
      ENDDO
 999  RETURN
      END
C#################################################
C ROUTINES TO CALCULATE CONSTRAINTS AND THEIR DERIVATIVES
C#################################################
      SUBROUTINE FY(IMODE,NJET,PFIT)
      IMPLICIT NONE

#include "abcfit_setup.inc"
#include "abcfit_bmatrix.inc"

      INTEGER IMODE,NJET
      DOUBLE PRECISION PFIT(4,NPARTICLES)

      INTEGER ICONS,I,J,K
      DOUBLE PRECISION DTEMP(4)
     
      CALL DVSET(NBFMX,0.0D0,FX(1),FX(2))
      
      DO I = 1, NJET
        FX(1) = FX(1) + PFIT(1,I) ! px-conservation
        FX(2) = FX(2) + PFIT(2,I) ! py-conservation
        FX(3) = FX(3) + PFIT(3,I) ! pz-conservation
        FX(4) = FX(4) + PFIT(4,I) ! e-conservation
      ENDDO
c      PRINT*,'FX IN',(FX(I),I=1,4)
c      PRINT*,'PSYS ',(PSYS(I),I=1,4)
      ICONS=1
      IF (ABS(PSYS(1)).LE.ABS(PSYS(4))) THEN
        FX(ICONS) = FX(1) - PSYS(1)
        ICONS=ICONS+1
      ENDIF
      IF (ABS(PSYS(2)).LE.ABS(PSYS(4))) THEN
        FX(ICONS) = FX(2) - PSYS(2)
        ICONS=ICONS+1
      ENDIF
      IF (ABS(PSYS(3)).LE.ABS(PSYS(4))) THEN
        FX(ICONS) = FX(3) - PSYS(3)
        ICONS=ICONS+1
      ENDIF
      IF (PSYS(4).GT.0.0D0) THEN
        FX(ICONS) = FX(4) - PSYS(4)
        ICONS=ICONS+1
      ENDIF
c      PRINT*,'FX OUT',(FX(I),I=1,4)
      
C
C MASS CONSTRAINTS
C
      IF (IMODE.GE.2.AND.IMODE.LE.7) THEN
        J=1
        DO I=1,4
          DTEMP(I)=0.0D0
        ENDDO
        DO WHILE (MINDEX(J,1).GT.0)
          DO I=1,4
            DTEMP(I)=DTEMP(I)+PFIT(I,MINDEX(J,1))
          ENDDO
          J=J+1
        ENDDO
        FX(5) = SQRT(DTEMP(4)**2-DTEMP(1)**2-DTEMP(2)**2-DTEMP(3)**2
     &    )
     
        IF (IMODE.GE.4) THEN
          J=1
          DO I=1,4
            DTEMP(I)=0.0D0
          ENDDO
          DO WHILE (MINDEX(J,2).GT.0)
            DO I=1,4
              DTEMP(I)=DTEMP(I)+PFIT(I,MINDEX(J,2))
            ENDDO
            J=J+1
          ENDDO
          FX(6) = SQRT(DTEMP(4)**2-DTEMP(1)**2-DTEMP(2)**2-DTEMP(3
     &      )**2)
        ENDIF
        IF (IMODE.EQ.2) FX(ICONS)=FX(5)-PSYS(5)
        IF (IMODE.EQ.3) FX(ICONS)=ALPHA_IJ*FX(5)-PSYS(5)
        IF (IMODE.EQ.4) FX(ICONS)=FX(5)-FX(6)
COB NEW EQUAL MASS CONSTAINT
        IF (IMODE.EQ.9) FX(ICONS)=FX(5)-FX(6)-PSYS(5)
        IF (IMODE.EQ.5) FX(ICONS)=ALPHA_IJ*FX(5)-ALPHA_KL*FX(6)
        IF (IMODE.EQ.6) THEN
          FX(ICONS)=FX(5)-PSYS(5)
          FX(ICONS+1)=FX(6)-PSYS(6)
        ENDIF
        IF (IMODE.EQ.7) THEN
          FX(ICONS)=ALPHA_IJ*FX(5)-PSYS(5)
          FX(ICONS+1)=ALPHA_KL*FX(6)-PSYS(6)
        ENDIF
      ELSEIF (IMODE.EQ.8) THEN
C
C EQUAL ENERGY CONSTRAINT
C
        DO K=1,2
          DTEMP(K)=0.0D0
          J=1
          DO WHILE (EINDEX(J,K).GT.0)
            DTEMP(K)=DTEMP(K)+PFIT(4,EINDEX(J,K))
            J=J+1
          ENDDO
        ENDDO
        FX(ICONS) = DTEMP(1)-DTEMP(2)

      ENDIF
      
 999  RETURN
      END
C
      SUBROUTINE DF_DP(IMODE,NJET,P,DFDP)
      IMPLICIT NONE

#include "abcfit_setup.inc"
#include "abcfit_bmatrix.inc"

      INTEGER IMODE,NJET
      DOUBLE PRECISION P(4,NPARTICLES),DFDP(NBFMX,4*NPARTICLES+3)

      INTEGER ICONS,I,J

      CALL DMSET(NBFMX,4*NJET+2,0.0D0,DFDP(1,1),DFDP(1,2),DFDP(2,1))
      
      DO I = 1, NJET
        ICONS=1
        IF (ABS(PSYS(1)).LE.ABS(PSYS(4))) THEN
          DFDP(ICONS,4*(I-1)+1) = 1.0d0 ! px-conservation
          ICONS=ICONS+1
        ENDIF
        IF (ABS(PSYS(2)).LE.ABS(PSYS(4))) THEN
          DFDP(ICONS,4*(I-1)+2) = 1.0d0 ! py-conservation
          ICONS=ICONS+1
        ENDIF
        IF (ABS(PSYS(3)).LE.ABS(PSYS(4))) THEN
          DFDP(ICONS,4*(I-1)+3) = 1.0d0 ! pz-conservation
          ICONS=ICONS+1
        ENDIF
        IF (PSYS(4).GT.0.0D0) THEN
          DFDP(ICONS,4*(I-1)+4) = 1.0D0 ! e-conservation
          ICONS=ICONS+1
        ENDIF
      ENDDO
C
C MASS CONSTRAINTS
C
      IF (IMODE.GE.2.AND.IMODE.LE.7) THEN
        CALL DM_DP(NJET,ICONS,MINDEX(1,1),P,DFDP)
        IF (IMODE.GE.4) THEN
          CALL DM_DP(NJET,ICONS+1,MINDEX(1,2),P,DFDP)
          IF (IMODE.EQ.4) THEN
            DO I=1,4*NJET
              DFDP(ICONS,I)=DFDP(ICONS,I)-DFDP(ICONS+1,I)
            ENDDO
          ENDIF
        ENDIF      
        IF (IMODE.EQ.3.OR.IMODE.EQ.5.OR.IMODE.EQ.7) THEN
          DO I=1,4*NJET
            DFDP(ICONS,I)=ALPHA_IJ*DFDP(ICONS,I)
            IF (IMODE.EQ.5) DFDP(ICONS,I)=DFDP(ICONS,I)-ALPHA_KL
     &        *DFDP(ICONS+1,I)
            IF (IMODE.EQ.7) DFDP(ICONS+1,I)=ALPHA_KL*DFDP(ICONS+1,I)
          ENDDO
          DFDP(ICONS,4*NJET+1) = DALPHA_IJ*DFDP(ICONS,4*NJET+1)
          IF (IMODE.EQ.5) DFDP(ICONS,4*NJET+2) = - DALPHA_KL*DFDP(ICONS
     &      +1,4*NJET+1)
          IF (IMODE.EQ.7) THEN
            DFDP(ICONS+1,4*NJET+2) =  DALPHA_KL*DFDP(ICONS+1,4*NJET+1)
            DFDP(ICONS+1,4*NJET+1) = 0.0D0
          ENDIF
        ENDIF
      ELSE IF (IMODE.EQ.8) THEN
C
C EQUAL ENERGY CONSTRAINT
C
        DO I=1,2
          J=1
          DO WHILE (EINDEX(J,I).GT.0)
            DFDP(ICONS,4*EINDEX(J,I)) = DBLE(3-2*I)
            J=J+1
          ENDDO
        ENDDO
      ENDIF
 999  RETURN
      END

      SUBROUTINE DM_DP(NJET,IROW,INDEX,P,DMDP)
      IMPLICIT NONE

#include "abcfit_setup.inc"

      INTEGER NJET,IROW,INDEX(NPARTICLES)
      DOUBLE PRECISION P(4,NPARTICLES),DMDP(NBFMX,4*NPARTICLES+3)

      INTEGER I,J
      DOUBLE PRECISION DTEMP(4),DTEMP1
           
      DO I=1,4
        DTEMP(I)=0.0D0
      ENDDO

      I=1
      DO WHILE (INDEX(I).GT.0.AND.INDEX(I).LT.8)
        DO J=1,4
          DTEMP(J)=DTEMP(J)+P(J,INDEX(I))
        ENDDO
        I=I+1
      ENDDO

      DTEMP1=SQRT(MAX(0.0D0,DTEMP(4)**2-DTEMP(1)**2-DTEMP(2)**2
     &  -DTEMP(3)**2))
      IF (DTEMP1.GT.0.0D0) THEN
        I=1
        DO WHILE (INDEX(I).GT.0.AND.INDEX(I).LT.NPARTICLES)
          DO J=1,3
            DMDP(IROW,(INDEX(I)-1)*4+J)=-DTEMP(J)/DTEMP1
          ENDDO
          DMDP(IROW,(INDEX(I)-1)*4+4)=DTEMP(4)/DTEMP1
          I=I+1
        ENDDO
      ENDIF

      DMDP(IROW,NJET*4+1)=DTEMP1
      
      RETURN
      END
*================================================================
*               Transformation function Gauss-Breit-Wigner
*================================================================
*
* This is original MATHKINE code.
*
*================================================================
*       G=Beta(alpha) function to use Breit Wigner fit
*               See ALEPH note ALEPH 95-129
*       ALPHA  = Gaussian variable
*================================================================
      FUNCTION BETA(J,X)
      IMPLICIT NONE

#include "abcfit_setup.inc"
#include "abcfit_bmatrix.inc"

      INTEGER J
      DOUBLE PRECISION  X,E,FUN_SSQRT,BETA,PI,DERF,DABS

      PARAMETER (PI=3.141592653589793238)
      E = EXP(1.D0)
      DELTA = GAMMA(J)*FUN_SSQRT(4*BETA0**2 - GAMMA(J)**2)

      IF ( DELTA.EQ.0.D0 ) DELTA = 1.D-12
      LAMBDA =0.25*PI + 0.5*ATAN((2.*BETA0**2 - GAMMA(J)**2)/DELTA)

      XZETA =ATAN(GAMMA(J)**2/DELTA) + LAMBDA*DERF(0.7071067811865475*(-
     &     ALPHA0 + X)/SIGMA)

      IF (DABS(XZETA-PI/2D0).LT.1D-12) THEN
        IF (XZETA.LT.PI/2D0) XZETA = PI/2D0-1D-12
        IF (XZETA.GE.PI/2D0) XZETA = PI/2D0+1D-12
      ENDIF
      BETA = FUN_SSQRT(BETA0**2 - 0.5*GAMMA(J)**2 + 0.5*DELTA*TAN(XZETA)
     $     )
      IF ( BETA.EQ.0D0 ) BETA = 1.D-12
      RETURN

      END
*================================================================
*       dBeta(alpha)/d(alpha) derivative of the BETAfunction
*               See ALEPH note ALEPH 95-129
*       ALPHA  = Gaussian variable
*================================================================
      FUNCTION DBETADX(J,X)
      IMPLICIT NONE

#include "abcfit_setup.inc"
#include "abcfit_bmatrix.inc"

      INTEGER NERROR,J
      DOUBLE PRECISION  X,E,BETA,DBETADX,GX,PI,DERF,DABS,DSQRT
      DOUBLE PRECISION FUN_SECANTE,DX
      PARAMETER (PI=3.141592653589793238)

      DATA NERROR/0/

      E = EXP(1.D0)
      DELTA = GAMMA(J)*DSQRT(4*BETA0**2 - GAMMA(J)**2)

      IF ( DELTA.EQ.0.D0 ) DELTA = 1.D-8
      LAMBDA =0.25*PI + 0.5*ATAN((2.*BETA0**2 - GAMMA(J)**2)/DELTA)

      XZETA =ATAN(GAMMA(J)**2/DELTA) + LAMBDA*DERF(0.7071067811865475*(-
     &ALPHA0 + X)/SIGMA)

      IF (DABS(XZETA-PI/2D0).LT.1D-8) THEN
        IF ( NERROR.LE.10 ) THEN
          PRINT *,'ABCFIT ERROR: INVALID TAN ARG ',XZETA
          IF ( NERROR.EQ.10 )
     &      PRINT *,'STOP SMALL TANGEANTE ERROR OUTPUT'
          NERROR = NERROR+1
        ENDIF
        IF (XZETA.LT.PI/2D0) XZETA = PI/2D0-1D-8
        IF (XZETA.GE.PI/2D0) XZETA = PI/2D0+1D-8
      ENDIF
      DX = (-alpha0 + x)**2/(2.*sigma**2)
      IF ( DX.GE.20.D0 ) THEN
        DBETADX = 0.D0
      ELSE
        GX = BETA(J,X)
        DBETADX = 0.3535533905932737*DELTA*LAMBDA*FUN_SECANTE(XZETA)**2/
     &       (E**DX*GX*DSQRT(PI)*SIGMA)
      ENDIF
      RETURN
      END
*================================================================
*               Compute the secante: 1/cos(x)
*================================================================
      FUNCTION FUN_SECANTE(X)
      DOUBLE PRECISION FUN_SECANTE,X,DCOS,DABS
      FUN_SECANTE = DCOS(X)
      IF ( DABS(FUN_SECANTE).GT.1.D-8 ) THEN
        FUN_SECANTE = 1.D0/FUN_SECANTE
      ELSE
        IF ( FUN_SECANTE.GT.0D0 ) THEN
          FUN_SECANTE = 1.D8
        ELSE
          FUN_SECANTE = -1.D8
        ENDIF
      ENDIF
      RETURN
      END
*================================================================
*               Safe double precision SQRT
*================================================================
      DOUBLE PRECISION FUNCTION FUN_SSQRT(X)
      IMPLICIT NONE
      DOUBLE PRECISION X,OLDVAL
      INTEGER NSMALL,NBIG
      DATA OLDVAL,NSMALL,NBIG/1.D0,1,1/
      SAVE OLDVAL,NSMALL,NBIG
      IF ( X.LT.0.D0 ) THEN
        IF ( X.NE.OLDVAL ) THEN
          OLDVAL = X
          IF ( X.GT.-1.D-4 ) THEN
            IF ( NSMALL.LE.10 ) THEN
              NSMALL = NSMALL + 1
              PRINT *,'FUN_SSQRT ERROR: NEGATIVE SQRT ',X
            ENDIF
            IF ( NSMALL.EQ.10 ) THEN
              PRINT *,'STOP FUN_SSQRT SMALL ERROR OUTPUT'
            ENDIF
          ELSE
            IF ( NBIG.LE.10 ) THEN
              NBIG = NBIG + 1
              PRINT *,'FUN_SSQRT ERROR: NEGATIVE SQRT ',X
            ENDIF
            IF ( NBIG.EQ.10 ) THEN
              PRINT *,'STOP FUN_SSQRT BIG ERROR OUTPUT'
            ENDIF
          ENDIF
        ENDIF
        FUN_SSQRT = -DSQRT(-X)
      ELSE
        FUN_SSQRT = DSQRT(X)
      ENDIF
      RETURN
      END





