CC HBOOK STUFF
C      INTEGER UNIT,ISTAT
C      INTEGER HMEMOR,HSIZE
C      PARAMETER (HSIZE=5000000)
C      COMMON/PAWC/HMEMOR(HSIZE)
C      INTEGER IQUEST
C      COMMON /QUEST/ IQUEST(100)
C      SAVE
C      CALL HLIMIT(HSIZE)


      INTEGER MAXNUMOFPAR
      PARAMETER (MAXNUMOFPAR=7)

      INTEGER ID_IN
      PARAMETER (ID_IN=1234)
      
      INTEGER IREC
      PARAMETER (IREC=1024)

      INTEGER NCFPARM
      PARAMETER (NCFPARM=1+MAXNUMOFPAR+1+1+2*4*MAXNUMOFPAR) ! NUMBER OF I*4 WORDS
      CHARACTER*500 CHFORM_CFPARM
      PARAMETER (CHFORM_CFPARM=
     &  'NJET[0,7]:I,'//
     &  'PID(NJET)[0,100]:I,'//
     &  'EVFLAG[-200,200]:I,'//
     &  'CUT:R,'//
     &  'PTRUE(4,NJET):R,'//
     &  'PRECO(4,NJET):R '
     &  )

      INTEGER NJET,PID,EVFLAG
      REAL CUT,PTRUE,PRECO

      COMMON /CFPARM/ NJET,PID(MAXNUMOFPAR),EVFLAG,CUT,PTRUE(4
     &  ,MAXNUMOFPAR),PRECO(4,MAXNUMOFPAR)


      