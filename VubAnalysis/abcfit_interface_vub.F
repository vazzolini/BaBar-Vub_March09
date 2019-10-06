      SUBROUTINE ABCFIT_INTERFACE_VUB(ISMEAR,ILTYP,CVAL,P_REC,P_FIT,
     &     CHI2T,PROBCHI,IERR,ISV)

C-----------------------------------------------------------------------
      IMPLICIT NONE 
C------------------ABCFIT INCLUDE FILE(S)------------------ 
#include "abcfit_setup.inc"
#include "abcfit.inc"
#include "abcfit_bmatrix.inc"
#include "abcfit_event-by-event.inc"
C------------INPUT PARAMETERS FOR INTERFACE--------------
      REAL    CVAL(4),P_REC(4,NPARTICLES)
      INTEGER ILTYP,ISMEAR
C------------OUTPUT PARAMETERS FROM INTERFACE--------------
      REAL    CHI2T, PROBCHI ,P_FIT(4,NPARTICLES)
      INTEGER IERR, ISV
C--------------INPUT PARAMETERS FOR ABCFIT--------------
C
      INTEGER ITF,ITYPP(NPARTICLES),ITEVOL,NUP
      REAL    M0(3),G0(3)
C--------------OUTPUT PARAMETERS FOR ABCFIT-------------
      INTEGER NDF
C--------------OUTPUT PARAMETERS FOR CY02PM-------------
      INTEGER NDIM
      REAL ERRMTR(PLEN,PLEN),PULL(PLEN),COVP(4*NPARTICLES,4*NPARTICLES)
     &     ,COVM(2,2)
C     SET VALUES FOR MASS HYPOTHESES
      REAL MPI,MMU,MK0,MKC,MASS(NPARTICLES),PI
      EXTERNAL PROB 
      REAL PROB
      DATA MPI/0.13967D0/
      DATA MMU/0.10566D0/
      DATA MK0/0.49767D0/
      DATA MKC/0.493677D0/
      DATA PI/3.1415927D0/
      INTEGER I,J,K,L
#include "ntparam.inc"
      INTEGER UNIT
      CHARACTER*80 NAME
      INTEGER ICALL
      DATA ICALL /0/
C     
C     SET NUMBER OF UNMEASURED PARTICLES LOCATED AT THE END OF P_REC ARRAY
C     
      NUP=1 ! MISSING NEUTRINO
      
C     
C     SET TOTAL NUMBER OF PARTICLES IN FIT
C     
      
      NJET=4 ! (BRECO,L,X,NEUTRINO)
C     
C     WHICH CONSTRAINS TO USE??
C     
      ITF=4                     ! EQUAL MASS CONSTAINT
C For 3C fit, use the following and input correct G0(I) for ITF=7 (see below)
C      ITF=6                     ! (M1=M0(1) ; M2=M0(2))
C      ITF=7                     ! (M1~M0(1) ; M2~M0(2)) mit Gauss verschmiert
C SET MASSINDEX := WHICH PARICLES FORM A RESONANCE
      
C     PARTICLE 1 GIVES B1=BRECO
      MINDEX(1,1)=1             ! B RECO CAND ONLY
      MINDEX(2,1)=0
C     PARTICLE 2-4 GIVES B2=BSEMIL
      MINDEX(1,2)=2             ! LEPTON
      MINDEX(2,2)=3             ! X-SYSTEM
      MINDEX(3,2)=4             ! NEUTRINO
      MINDEX(4,2)=0
      
C     
C     WHICH PARAMETRISATION DATA TO USE? DO WE WANT NON-DIAGONAL
C COVARIANCE? TRUE OR RECO BINNING? WHICH *TAG* SHOULD WE USE?
C     HERE WE DECIDE TO USE A FILE WITH TAG=1, TRUE BINNING AND DIAGONAL
C COVARIANCE MATRIX
C     
C     +---- TRUE=1
C     |+--- DIAGONAL=0
C     ||+-- TAG=1
C     VVV
C     ITEVOL=
      
      IF (ILTYP.EQ.2) THEN      ! ELECTRON
         ITEVOL=002
         MASS(2)=0.000511       ! SET MASS OF ELECTRON ZERO
      ELSEIF (ILTYP.EQ.1) THEN  ! MUON
         ITEVOL=001
         MASS(2)=MMU
      ENDIF   
C     
C     SET VALUES OF MASS CONSTRAINTS (IF USED -- WE SET THEM ANYWAY)
C
      DO I=1,3
         M0(I)=5.279
         G0(I)=0.030            ! JUST A GUESS
      ENDDO
C     
C     WHICH KIND OF PARAMETRISATION? WE TRY ALEPH (DELPHI IS RECOMMENDED)
C     AND HAVE FULL CORRECTION (BOTH TRANSVERSE AND LONGITUDINAL) AND
C     THE FITTED JET MASSES SCALE WITH FITTED MOMENTA
C     
C      +---- SCALED JET-MASSES=0 (OLD MATHKINE WAY)
C      |+--- FULL CORRECTION=0
C      ||+-- ALEPH=0
C      VVV
CITYPP=

C ------ BRECO CANDIATE      
      ITYPP(1)=20               ! 4TH PARAMTER ALA ALEPH (30 ALA DELPHI)
C ------ LEPTON
      ITYPP(2)=10               ! LEPTON FIXED MASS 
C ------ X-System
      ITYPP(3)=0                ! SIMPLE RESCALING : EFIT=PFIT * ERECO/PRECO
C ------ NEUTRINO
      ITYPP(4)=0                ! JUST TO MAKE SURE Mv STAY ZERO IN THE FIT
      MASS(4)=0.0
      
C MAKE SURE THAT MASS HYPOTHESIS FOR LEPTON AND NEUTRINO
C IS SET CORRECT AND GIVE "GOOD STARTING VALUES" FOR THE
C NEUTRINO      
      CALL VZERO(P_FIT,4*NUP)
      DO J=1,NJET-NUP
         DO I=1,4
            P_FIT(I,1)=P_FIT(I,1)+P_REC(I,J)
         ENDDO   
      ENDDO
C MISSING MOMENTUM   
      DO I=1,3
         P_REC(I,NJET)=CVAL(I)-P_FIT(I,1) 
      ENDDO 
C MISSING ENERGY
      P_FIT(4,2)=CVAL(4)-P_FIT(4,1)
C MISSING MASS**2
      MMISS= P_FIT(4,2)**2 - 
     &            (P_REC(1,I)**2+P_REC(2,I)**2+P_REC(3,I)**2) 
      IF ( ISV.NE.0) THEN 
         PRINT*,'MMISS=',MMISS
      ENDIF
C SET MASS HYPOTHEIS FOR LEPTON AND NEUTRINO
C FIRST LEPTON (POS 2)
         P_REC(4,2)=SQRT( (P_REC(1,2)**2+P_REC(2,2)**2+P_REC(3,2)**2)
     &        +MASS(2)**2)
C SECOND NEUTRINO (POS 4)
         P_REC(4,4)=SQRT( (P_REC(1,4)**2+P_REC(2,4)**2+P_REC(3,4)**2)
     &        +MASS(4)**2)
 
      IF ( ISV.NE.0) THEN 
         DO I=3,4
            PRINT*,'mass',I,MASS(I)
         ENDDO
      ENDIF
     
C---------PARAMETER INITIALISATION FOR ABCFIT-----------
C DON'T CALL AIBI_EVOL_USER TO GET PARAMETERISATION
C USE AIBI_EVOL FROM ABCFIT TO GIVE PARAMETRISATION DATA
      USER_AIBI_EVOL=.FALSE.

      SMOOTHVARFLAG='MMISS'
      IF (ISMEAR.eq.0) THEN
         ABCTAG='mmiss '
      ENDIF
      IF (ISMEAR.eq.-1) THEN
         ABCTAG='flat '
C        Reset to expected value
         ISMEAR=0
      ENDIF


      IF (ISMEAR.GT.100) THEN
         IF (ISMEAR .EQ. 101) THEN
            SMOOTHVARFLAG='MMISS'
            ABCTAG = 'mmiss '
         ENDIF
         IF (ISMEAR .EQ. 102) THEN
            SMOOTHVARFLAG=' '
            ABCTAG = 'costh '
         ENDIF
         IF (ISMEAR .EQ. 103) THEN
            SMOOTHVARFLAG='E'
            ABCTAG = 'energ '
         ENDIF
         ISMEAR = 0
      ENDIF


C CHECK FOR SMEARING OPTION!
      IF (ISMEAR.GT.0) THEN
         NJET=1
         NUP=0
         ITF=1
         ITEVOL=001
         SMOOTHVARFLAG='E'
         IF (ISMEAR.EQ.1) THEN ! DC B0 CAND SMEARING
            ABCTAG='dc '
         ELSEIF (ISMEAR.EQ.2) THEN ! DSTAR B0 CAND SMEARING
            ABCTAG='dstar '
         ELSEIF (ISMEAR.EQ.3) THEN ! D0 B+- CAND SMEARING
            ABCTAG='d0 '
         ELSEIF (ISMEAR.EQ.4) THEN ! DSTAR0 B+- CAND SMEARING
            ABCTAG='dstar0 '
         ENDIF
         ISMEAR=1
      ENDIF  
C       PRINT*,'ITF,CVAL',ITF,CVAL
C       PRINT*,'P_REC',P_REC
       IF ( ISV.NE.0) THEN 
         PRINT *, 'BEFORE THE KINEMATIC FIT:'
         DO I=1,4
            PRINT*,'PFIT',I,(P_FIT(J,I),J=1,4)
         ENDDO
       ENDIF
       CALL ABCFIT(ISMEAR,NJET,NUP,P_REC,ITF,CVAL,M0,G0,ITYPP,
     &      ITEVOL,P_FIT,CHI2T,NDF,IERR)
       IF ( ISV.NE.0) THEN 
         PRINT *, 'AFTER THE KINEMATIC FIT:'
         DO I=1,4
            PRINT*,'PFIT',I,(P_FIT(J,I),J=1,4)
         ENDDO
       ENDIF
C SET PARAMTERS BACK FOR NORMA FITTING
       IF (ISMEAR.EQ.1) THEN
          ABCTAG=' '
          SMOOTHVARFLAG=' '
       ENDIF   
      
C     GET COV MATRIX FOR P AND M AFTER FIT IF FIT IS OK
      PROBCHI = -1
      IF (IERR.GE.0) THEN
         IF (CHI2T.GT.0) THEN
           PROBCHI = PROB(CHI2T, NDF)
         ENDIF
         CALL CY02PM(MINDEX,I,ERRMTR,PULL,COVP,COVM) ! ERROR
      ENDIF   
      
      
      RETURN
      END

      
C--------------------------------------------------------------------
      SUBROUTINE OPEN_FILE(UNIT,NAME,FLAG)
C--------------------------------------------------------------------

      IMPLICIT NONE     

      LOGICAL LEXIST
      INTEGER IENERGY,FLAG,UNIT
      CHARACTER*80 NAME
     


C     WRITE THE OUTPUT TO FILE
      
C     CHECK IF FILE EXISTS

      INQUIRE (FILE=NAME,EXIST=LEXIST)
      
      IF (LEXIST) THEN
        
         PRINT*,'FILE ',NAME,' ALREADY EXISTS'
         IF (FLAG.EQ.1.OR.FLAG.EQ.2) THEN
           PRINT*,'THE PROGRAM WILL BE TERMINATED!!!!!!!!!! '
           PRINT*,'IF YOU WANT TO OVERWRITE THE FILE'
           PRINT*,'=> JUST DELETE THE FILE ON YOU OWN'
           STOP
         ENDIF
         
         OPEN (UNIT,FILE=NAME,ACCESS='SEQUENTIAL'
     & ,FORM='FORMATTED',STATUS='OLD')
       ELSE
         
         IF (FLAG.NE.2) THEN
           
C     OPEN A NEW FILE AND KEEP IT
           
           PRINT*,'FILE ',NAME,' DOES NOT EXIST'
           PRINT*,'OPEN NEW FILE AND KEEP IT'
           
           OPEN (UNIT,FILE=NAME,ACCESS='SEQUENTIAL'
     &          ,FORM='FORMATTED',STATUS='NEW')
  
         ELSEIF (FLAG.EQ.2) THEN
 
           PRINT*,'OPTION DISP= IS NOT SUPPORTED IN LINUX!'
           PRINT*,'SORRY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           STOP
          
C     OPEN A NEW FILE AND DELETE IT AFTER USE
           
c           PRINT*,'FILE ',NAME,' DOES NOT EXIST'
c           PRINT*,'OPEN NEW FILE AND DELETE IT AFTER USE'

c           OPEN (UNIT,FILE=NAME,ACCESS='SEQUENTIAL'
c     &          ,FORM='FORMATTED',STATUS='NEW',DISP='DELETE')   
           
                      
         ENDIF

       ENDIF  
    

      RETURN
      END

C-----------------------------------------
      SUBROUTINE WHICH_UNIT(UNT)    
C-----------------------------------------      
C
C GET UNIT TO USE
C
      IMPLICIT NONE
      
C
C     OUTPUT
C
      INTEGER UNT
C
C     INTERNAL
C      
      LOGICAL USED
      
      UNT=90
      USED=.TRUE.
      DO WHILE (USED)
        INQUIRE(UNIT=UNT,OPENED=USED)
        UNT=UNT+1
      ENDDO
      UNT=UNT-1
      
      RETURN
      END
      






