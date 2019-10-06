      INTEGER NPAR(-1:19)
      DATA NPAR /0,
     &     3,3,3,3,1,3,0,0,0,0, ! BASIC
     &     1,1,1,1,1,1,0,0,0,0/ ! 4TH PARAMETER

      SAVE NPAR

      INTEGER MINDEX,EINDEX,IPAM,NUMJET,NCONS,NPARAM
      COMMON /BMATRXI/ MINDEX(NPARTICLES,3),EINDEX(NPARTICLES,3)
     &     ,IPAM(NPARTICLES),NUMJET,NCONS,NPARAM(NPARTICLES)

      DOUBLE PRECISION Pim,Ptim,Eim,MIM,BIM,PIR,EIR,MIR,BIER
      DOUBLE PRECISION urx,ury,urz
      DOUBLE PRECISION uthx,uthy,uthz
      DOUBLE PRECISION uphx,uphy,uphz
      DOUBLE PRECISION ai,bi,ci,di
      DOUBLE PRECISION MZ,GZ
      DOUBLE PRECISION BETA0,ALPHA0,GAMMA,DELTA,LAMBDA,SIGMA,XZETA
      DOUBLE PRECISION alphaij,alphakl,ALPHA_IJ,ALPHA_KL,DALPHA_IJ
     $     ,DALPHA_KL
      DOUBLE PRECISION  PSYS,B,FX

      COMMON /BMATRXD/ 
     &     Pim(NPARTICLES),Ptim(NPARTICLES),Eim(NPARTICLES)
     &	   ,MIM(NPARTICLES)
     &     ,BIM(NPARTICLES),PIR(NPARTICLES),EIR(NPARTICLES)
     &     ,MIR(NPARTICLES),BIER(NPARTICLES),urx(NPARTICLES)
     &     ,ury(NPARTICLES),urz(NPARTICLES),uthx(NPARTICLES)
     &     ,uthy(NPARTICLES),uthz(NPARTICLES),uphx(NPARTICLES)
     &     ,uphy(NPARTICLES),uphz(NPARTICLES),ai(NPARTICLES)
     &     ,bi(NPARTICLES),ci(NPARTICLES),di(NPARTICLES),MZ(3),GZ(3)
     &     ,BETA0,ALPHA0,GAMMA(2),DELTA,LAMBDA,SIGMA,XZETA,alphaij
     &     ,alphakl,ALPHA_IJ,ALPHA_KL,DALPHA_IJ,DALPHA_KL,PSYS(NBFMX)
     &     ,B(NBFMX,PLEN),FX(NBFMX)
