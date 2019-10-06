C=============================================================
C     Routine to book, fill and output ntuples with
C     information for making parameterisations
C     ALPHA-independent version
C
C IMPORTANT: Remember to call NTPARAM WITH IFLAG=1 AFTER _LAST_
C            event to terminate ntuple properly.
C=============================================================
      SUBROUTINE NTPARAM(IFLAG,NJET_IN,PTRUE_IN,PRECO_IN,PID_IN
     &  ,EVFLAG_IN,CUT_IN)
C
C Input:
C     IFLAG      = .NE.1 NO MEANING
C                   1  TERMINATE NTUPLE
C
C     NJET_IN    = NUMBER OF PARTICLES
C     PTRUE_IN   = _TRUE_ 4-MOMENTA (PX,PY,PZ,E) OF UP TO 7 PARTICLES
C     PRECO_IN   = _RECO_ 4-MOMENTA (PX,PY,PZ,E) OF UP TO 7 PARTICLES
C     PID_IN     = PARTICLE-FLAG OF UP TO 7 PARTICLES (E.G. 1123 for jjlv)
C     EVFLAG_IN  = INTEGER FLAG TO TAG EVENT-TYPE
C     CUT        = ADDITIONAL REAL PARAMETER FOR EVENT-SELECTION (E.G. NN-OUTPUT)
C
      IMPLICIT NONE

      INTEGER IFLAG

      INTEGER NJET_IN,EVFLAG_IN,PID_IN(7)
      REAL PTRUE_IN(4,7),PRECO_IN(4,7),CUT_IN

#include "ntparam.inc"

      LOGICAL FIRST,USED
      INTEGER I,J,UNT
      CHARACTER*100 OLDDIR
      CHARACTER*500 CTEMP

      DATA FIRST /.TRUE./
      SAVE FIRST,UNT

C      PRINT*,'IFLAG,NJET_IN,PTRUE_IN,PRECO_IN,PID_IN
C     &  ,EVFLAG_IN,CUT_IN',IFLAG,NJET_IN,PTRUE_IN,PRECO_IN,PID_IN
C     &  ,EVFLAG_IN,CUT_IN
      
      CALL HCDIR(OLDDIR,'R')

      IF (FIRST) THEN
        FIRST=.FALSE.
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

C        PRINT*,'UNT,IREC,I',UNT,IREC,I
        CALL HROPEN(UNT,'PAROUT','ntparam.hist','N',IREC,I)
C        PRINT*,'BOOKIMG',I
        
        IF (I.NE.0) THEN
          PRINT*,'**********************************************'
          PRINT*,'**********************************************'
          PRINT*,'ERROR BY OPENING FILE ntparam.his'
          PRINT*,'**********************************************'
          PRINT*,'**********************************************'
          STOP
        ENDIF

        CALL HBNT(ID_IN,'PAROUT',' ')
        CTEMP=CHFORM_CFPARM
        I=INDEX(CTEMP,' ')-1
        CALL HBNAME(ID_IN,'CFPARM',NJET,CTEMP(1:I))
      ELSE
        CALL HCDIR('//PAROUT',' ')
      ENDIF

      IF (IFLAG.EQ.1) THEN
        I=0
        CALL HROUT(ID_IN,I,' ')
        CALL HREND('PAROUT')
        CLOSE (UNT)
        GOTO 999
      ENDIF

      CALL VZERO(NJET,NCFPARM)

      NJET=NJET_IN
      EVFLAG=EVFLAG_IN
      DO I=1,NJET
        PID(I)=PID_IN(I)
        DO J=1,4
          PTRUE(J,I)=PTRUE_IN(J,I)
          PRECO(J,I)=PRECO_IN(J,I)
        ENDDO
      ENDDO
      CUT=CUT_IN
      
      CALL HFNT(ID_IN)
 999  CONTINUE
      CALL HCDIR(OLDDIR,' ')

      
      RETURN
      END






