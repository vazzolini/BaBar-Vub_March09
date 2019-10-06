      PROGRAM MAIN
      INTEGER I,IERR,ILTYP,ISMEAR,ISV
      REAL P_REC(4,4),P_FIT(4,4),CHI,CVAL(4),PROBCHI
C PARTICLE ONE
      P_REC(1,1)= -0.360311     ! PX
      P_REC(2,1)=9.39542E-02        ! PY
      P_REC(3,1)=  2.97556      ! PZ
      P_REC(4,1)= 6.07137       ! E
C PARTICLE TWO
      P_REC(1,2)= -0.891153     ! PX
      P_REC(2,2)= -0.921733    ! PY
      P_REC(3,2)= 0.859916      ! PZ
      P_REC(4,2)= 1.54376       ! E
C PARTICLE THREE
      P_REC(1,3)=   0.561458    ! PX
      P_REC(2,3)=  0.714287    ! PY
      P_REC(3,3)=   -0.159553     ! PZ
      P_REC(4,3)=   2.20856     ! E
C PARTICLE FOUR
      P_REC(1,4)= 0.574537     ! PX
      P_REC(2,4)=  0.117605     ! PY
      P_REC(3,4)=   2.19973     ! PZ
      P_REC(4,4)=    2.27656     ! E
      
      CVAL(1)=  -0.110255              ! SUM(PX)=CVAL(1)
      CVAL(2)= 0.              ! SUM(PY=CVAL(2)
      CVAL(3)= 5.87722                ! SUM(PZ=CVAL(3)
      CVAL(4)=  12.1018            ! SUM(E)=ETOT=CVAL(4)
      ISMEAR=0
      ILTYP=1
      CALL  ABCFIT_INTERFACE_VUB(ISMEAR,ILTYP,CVAL,P_REC,P_FIT,CHI2T,
     &   PROBCHI,IERR,ISV)
      
      PRINT*,'IERR',IERR,CHI2T
     
      END
      




