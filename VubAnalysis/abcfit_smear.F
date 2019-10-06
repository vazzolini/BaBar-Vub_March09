C===============================================================C
C                                                               C
C     AAAA     BBBB      CCCC     FFFFFF   IIIIII   TTTTTTTT    C
C    A    A    B   B    C    C    F          II        TT       C
C    A    A    B   B   C          F          II        TT       C
C    AAAAAA    BBBB    C          FFFFF      II        TT       C
C    A    A    B   B   C          F          II        TT       C
C    A    A    B   B    C    C    F          II        TT       C
C    A    A    BBBB      CCCC     F        IIIIII      TT       C
C                                                               C
C===============================================================C
C
C                            BY
C
C                 OLIVER BUCHMULLER AND #5
C
C                FORMERLY KNOWN AS MATHKINE
C
C
C
C HISTORY:
C
C Version 1.0 : August 1998, Oliver Buchmuller and #5
C               -Major upgrade/rewriting of mathkine 1.7.
C               -Exact implementation of unmeasured objects (i.e.
C                semi-leptonic fit -- neutrino is now an UNMEASURED object)
C               -Full covariance matrix of fitted objects implemented
C                (output is in CY02PM which also gives covariance matrix of
C                 fitted momenta and fitted masses).
C               -True parameterisation binning in the input parameters
C               -Full input covariance matrix is implemented everywhere.
C               -B-matrix calculation from mathkine has been replaced
C                completely to allow more general calculation in terms of
C                constrains, parameterisations and number of particles
C               -Numerous parameterisations implemented and it is now easy to
C                to introduce additional
C
C ABCFIT 1.x is backward compatible in _MAIN_ functionality
C AND _FULL_ reliability with the old mathkine 1.x.
C
C MAIN CODE IS NOW ALPHA INDEPENDENT -- ALPHA INTERFACE ROUTINES ARE PROVIDED.
C
C===============================================================C
      SUBROUTINE ABCFIT(ISMEAR_IN,NJET,NUP,P_REC,ITF,CVAL,M0,G0,PARID,
     & ITEVOL,P_FIT,CHI2T,NDF,IERR)
*==============================================================
*  In
*       NJET        = Number of objects in event =<7
*       NUP         = Number of Unmeasured Particles (e.g. neutrinos)
*                     in event (must be last in p_rec)
*                     (e.g. nup=0 for qqqq and nup=1 for qqlv)
*       P_REC(4,7)  = four-momentum of up to 7 objects
*       ITF         = mn
*         ---> n    = 1 Four momentum and energy conservation,
*                        Som{j=1,njet}p_fit(j) = cval
*              n    = 2 Energy momentum conservation
*                       + m(1,2) = m0(1)
*                       G0 is not used
*              n    = 3 Energy momentum conservation
*                       + alpha(1,2)*m(1,2) = m0(1)
*                       sigma(alpha) = G0(1)/m0(1)
*              n    = 4 Energy momentum conservation
*                       + m(1,2)-m(3,4) = 0
*              n    = 5 Energy momentum conservation
*                       + alpha(1,2)*m(1,2)-alpha(3,4)*m(3,4) = 0
*                       sigma(alpha(1,2)) = G0(1)/m0(1)
*                       sigma(alpha(3,4)) = G0(2)/m0(2)
*              n    = 6 Energy momentum conservation
*                       + m(1,2)=m0(1)  m(3,4)=m0(2)
*              n    = 7 Energy momentum conservation
*                       + alpha(1,2)*m(1,2)=m0(1)  alpha(3,4)*m(3,4)=m0(2)
*              n    = 8 Energy momentum conservation
*                       + E(1,2)=E(3,4) "equal energy"
*         -->  m    = 0 In case of mass constraint, gaussian fit is used
*                   = 1 In case of mass constraint, Breit Wigner fit is
*                       used, i.e the alpha,s parameters have a Breit-Wigner
*                       ditribution
*       CVAL(4)     = value of 4-momentum constrains Px,Py,Pz,E
*                     (e.g. 0,0,0,ecm). Can be used to switch of some of the
*                     constraints:
*                     If abs(cval(1)).gt.abs(cval(4)) Px constraint is not
*                     applied. Same holds for cval(2) (Py) and cval(3) (Pz).
*                     If cval(4).lt.0.0 E constraint is not applied, BUT
*                     remember to set |value| such that you don't switch of
*                     Px, Py and Pz constraints by accident....
*       M0(3)       = value of up to 3 mass-constrains (e.g. mw,mw)
*       G0(3)       = width of up to 3 mass-constrains (e.g. w-width)
*       PARID(7)    = kmn
*                n  = 0  Fit parameterisation a la ALEPH
*                   = 1  Fit parameterisation a la DELPHI
*                   = 2  Fit parameterisation using P,theta,phi 
*                   = 3  Fit parameterisation using Px,Py,Pz
*                m  = 0  Fitted jet energies scale with fitted momenta (old
*                        mathkine)
*                     1  Fitted jet energies is determined from fixed input 
*                        jet mass and fitted momenta
*                     2  Fit jet mass using a 4'th parameter a la ALEPH:
*                               M_fit = d*M_reco
*                     3  Fit jet mass using a 4'th parameter a la DELPHI:
*                               M_fit = exp(d)*M_reco
*                     4  Fit jet mass using a 4'th parameter:
*                               E_fit = (exp(d+log(E_reco/P_reco-1))+1)*P_fit
*                k  = 0  Full correction on jet momentum (// and T)
*                   = 1  Only correction on jet momentum (//)
*                   = 2  Full correction on jet momentum, but the longitudinal
*                        parameter is assumed to be unmeasured
*                   = 3  switch on k=1 _and_ k=2
*       ITEVOL      = kmn
*                n  > 0  Evolution with flavour/jet characteristics:
*                        The file 
*
*                           aibi_evol_<itypp*njet>_<itevol-k0n>.dat
*
*                        is used for parametrisation definition and parameter
*                        evolution e.g. aibi_evol_0000_001.dat for itypp=0
*                        (ALEPH) and njet=4 and itevol=001 (diagonal reco-
*                        binning)
*
*                m  = 0  Diagonal input covariance matrix
*                   = 1  Full non-diagonal input covariance matrix
*                k  = 0  binning in reco values for correction-factors
*                   = 1  binning in true values for correction-factors
*
*  Out
*       P_FIT(4,7)  = fitted four-momentum of up to NPARTICLES objects
*       CHI2T       = total chi2
*       NDF         = number of degree of freedom
*       IERR        = error flag
*
*  ERROR OUTPUT:
*       IERR        <  0 results are irrelevant
*                   >= 0 results fulfils constrains and can be used
*
*                   = -1 number of input jets/objects do not match the
*                        definition of one or more constrain(s)
*                   = -2 number of unmeasured parameters (=3*NUP) exeeds
*                        number of constrains
*                   = -3 undefined parameterisation requested
*                   = -4 aibi_evol(_user) routine failed to provide
*                        starting values and/or correlation matrix
*                   = -5 One measured jet has zero energy
*                   = -6 One fitted particle has ZERO longitudinal component
*                   = -7 One fitted particle has longitudinal component
*                        larger than PMAX (default=10000 GeV)
*                   = -8 Error in inverting BVBT (=Wb) (constraint "solution")
*                   = -9 Error in inverting Wa (unmeasured particle "solution")
*                   = -10 chi2 overflow
*                   = -11 Maximum number of iteration reached without
*                         fulfilment of constrains (with FXMAX in GeV
*                         (default 0.0001)). This is certainly
*                         a bad chi2 event.
*                   = -12 One fitted particle has mass greater than 
*                         PMAX (default=10000 GeV)
*                   = -20 True binning fit failed -- reco binning is used
*                   = -21 Too many iterations for true binning fit -- reco
*                         binning is used
*                   = -22 True binning fit failed with lower bound violation
*                         -- reco binning is used
*                   = -23 True binning fit failed with upper bound violation
*                         -- reco binning is used
*                   = -99 Kind of fit does not exist
*                   = 1 Error matrix invertion failed when v is non diagonal
*                       v is forced to be diagonal, and fit is performed
*                   = 2 Maximum iteration reached, result could be used
*                       but this is certainly a bad chi2 event
*                       The maximum number of iteration is NITERMX
*                       (default is 50).
*                   = 3 At least one measured AND one unmeasured object
*                       has been reversed by the fit
*                   = 4 At least one measured object has been reversed
*                       by the fit
*                   = 5 At least one measured object has been reversed
*                       by the fit
*                   = 6 Input object with beta==1 supposed to be fit
*                        with a "mass-related" 4'th parameter is fit
*                        with fixed beta==1 assumption
*                   = 8 True binning fit okay -- BUT with lower bound violation
*                   = 9 True binning fit okay -- BUT with upper bound violation
*
      IMPLICIT NONE

#include "abcfit_setup.inc"
#include "abcfit.inc"

c Input parameters
      INTEGER NJET,NUP,ITF,PARID(NPARTICLES),ITEVOL
      REAL    P_REC(4,NPARTICLES),CVAL(4),M0(3),G0(3)
c Output parameters
      INTEGER NDF,IERR
      REAL    P_FIT(4,NPARTICLES),CHI2T
c internal variables
      LOGICAL BOUND(2)
      INTEGER I,J,ITER,ITEVOLRECO,IERTRU
      REAL CHI2TOLD,DCHI2,P_COR(4,NPARTICLES)
      INTEGER ISMEAR_IN
      COMMON/SMEAR/ISMEAR
      INTEGER ISMEAR

      ISMEAR=ISMEAR_IN
C
C TOO MANY JETS?
C
      IF (NJET.GT.NPARTICLES) THEN
        IERR=-1
        RETURN
      ENDIF


C======================================================================
C      ENTRY ABCFIT_FULL(NJET,NUP,P_REC,ITF,CVAL,M0,G0,PARID,ITEVOL,P_FIT
C     &     ,CHI2T,NDF,IERR)
C
C ONLY NEW STUFF IS: ITYPP->PARID(NPARTICLES)
C
C WHERE PARID(I)=ITYPP for particle/jet I
C
C======================================================================
C
C True binning????
C
      IF (ITEVOL.GE.100) THEN
         ITEVOLRECO=MOD(ITEVOL,100) ! go for y0-corrections binned in true info
      ELSE
         ITEVOLRECO=ITEVOL      ! go for y0-corrections binned in reco info
      ENDIF
      BOUND(1)=.FALSE.
      BOUND(2)=.FALSE.
      IERTRU=0
      ITER=0
C
C save P_REC in P_COR(point corrections) for y0-corrections binned in true info
C
      DO I=1,NJET
        DO J=1,4
          P_COR(J,I)=P_REC(J,I)
        ENDDO
      ENDDO
C
C MAIN ITERATIVE MINIMISATION IF TRUE VALUES ARE WANTED
C            
 1234 CONTINUE

      CALL ABCFITKERNEL(1,NJET,NUP,P_REC,PARID,ITF,CVAL,M0,G0,P_COR
     &     ,ITEVOLRECO,P_FIT,CHI2T,NDF,IERR)
      
      IF (IERR.EQ.-11) CALL
     &     ABCFITRECOVER(NJET,NUP,P_REC,PARID,ITF,CVAL,M0,G0,P_COR
     &     ,ITEVOLRECO,P_FIT,CHI2T,NDF,IERR)

      IF (IERR.LT.0.OR.CHI2T.GE.CHI2MAX) GOTO 123 ! fit fails
      IF (ITER.GT.NITERMXTRUE) GOTO 123 ! fit with true info fails use reco info      
      IF (ITEVOL.LT.100) THEN
        IF (IERR.GE.8.AND.IERR.LE.9) IERR=0
        GOTO 123                ! use only reco info
      ENDIF

      CHI2TOLD=0.0
      DCHI2=100.0
      DO WHILE (ABS(DCHI2).GT.1D-2.AND.ITER.LE.NITERMXTRUE) ! fit using true info

        CHI2TOLD=CHI2T ! pick up the chi2 from the 'reco' fit
        DO I=1,NJET
          DO J=1,4
            P_COR(J,I)=P_FIT(J,I)
          ENDDO
        ENDDO

        CALL ABCFITKERNEL(1,NJET,NUP,P_REC,PARID,ITF,CVAL,M0,G0
     &    ,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)

        IF (IERR.EQ.-11) CALL ABCFITRECOVER(NJET,NUP
     &       ,P_REC,PARID,ITF,CVAL,M0,G0,P_COR,ITEVOLRECO,P_FIT,CHI2T
     &       ,NDF,IERR)

        IF (IERR.EQ.8) BOUND(1)=.TRUE. ! CHECK LOWER BOUND VIOLATION
        IF (IERR.EQ.9) BOUND(2)=.TRUE. ! CHECK UPPER BOUND VIOLATION
        IF (IERR.LT.0.OR.CHI2T.GE.CHI2MAX) THEN ! 'true' fit fails use 'reco' fit
          IERTRU=IERR
          ITER=201
        ENDIF

        DCHI2=ABS(CHI2TOLD-CHI2T)  ! monitor the diff in chi2
        ITER=ITER+1
        
      ENDDO

      IF (BOUND(1)) IERR = 8    ! LOWER BOUND VIOLATION
      IF (BOUND(2)) IERR = 9    ! UPPER BOUND VIOLATION

      IF (ITER.GT.NITERMXTRUE) GOTO 1234 ! no convergence of 'true' fit use reco
 123  CONTINUE 
      IF (ITER.GT.NITERMXTRUE.AND.IERTRU.NE.0) IERR = -20 ! True fit failed
      IF (ITER.GT.NITERMXTRUE.AND.IERTRU.EQ.0) THEN
        IERR = -21               ! NUMBER OF ITERATIONS EXEEDED
        IF (BOUND(1)) IERR = -22 ! LOWER BOUND VIOLATION
        IF (BOUND(2)) IERR = -23 ! UPPER BOUND VIOLATION
      ENDIF

      RETURN
      END

c======================================================================
      SUBROUTINE ABCFITKERNEL(IFLAG,NJET,NUP,P_REC,PARID,ITF,CVAL,M0,G0
     &  ,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)
c======================================================================

      IMPLICIT NONE

#include "abcfit_setup.inc"
#include "abcfit.inc"
#include "abcfit_bmatrix.inc"

c Input parameters
      INTEGER IFLAG,NJET,NUP,ITF,ITEVOL,PARID(NPARTICLES)
      REAL    P_REC(4,NPARTICLES),CVAL(4),M0(3),G0(3),P_COR(4,NPARTICLES
     &     )
c Output parameters
      INTEGER NDF,IERR
      REAL    P_FIT(4,NPARTICLES),CHI2T
c Internal variables
      INTEGER ITFMN(2),ITEVOLKMN(3)
      DOUBLE PRECISION Y(PLEN),Y0(PLEN),DY0(PLEN),V(PLEN,PLEN),VINV(PLEN
     &     ,PLEN)
      DOUBLE PRECISION VBT(PLEN,NBFMX),BVBT(NBFMX,NBFMX),WA(PLEN,PLEN)
     &     ,WAWB(NBFMX,NBFMX),YMY0(PLEN),Y_OLD(PLEN)
      DOUBLE PRECISION WORK(PLEN),Z(PLEN),WORK2(PLEN,PLEN)
      DOUBLE PRECISION CHI2TX,CHI2MN,DCHI2,F2
      INTEGER I,J,K,NID,NBF,NITER,IERT,NITERR,IBOUND,NIDREAL
     &  ,NPARAMS,NFREE,NLFREE,LFREEP(NPARTICLES)

      COMMON/SMEAR/ISMEAR
      INTEGER ISMEAR
      REAL RVEC(6)
C
C parse user flags and report setup
C
      CALL ABPUF(NJET,ITF,PARID,ITEVOL,ITFMN,ITEVOLKMN)
c valid parameterisations??
      NUMJET=0
      IERR = 0
      DO I=1,NJET
        IF (NPAR(MOD(PARID(I),10)).LE.0.OR.NPAR(MOD(PARID(I)/10,10)+10)
     &       .LE.0) IERR=-3
C
C setup particle parameterisation
C
        IPAM(I)=PARID(I)
C
C ALWAYS FIT NUPS WITH ALEPH PARAMETERISATION
C
c        IF (I.GT.NJET-NUP) IPAM(I)=PARID(I)-MOD(PARID(I),10)
      ENDDO
      IF ( IERR.NE.0  ) GOTO 999
c
c YES! now prepare for the fit...
C
C IS THE CONSTRAINS VALID FOR THE GIVEN NUMBER OF JETS
      IF (NJET.LT.2) IERR=-1
c determine number of parameters and constraints from user input
c
      IERR = -99
C
      IF (MINDEX(1,1).LE.0.OR.MINDEX(1,1).GE.NPARTICLES) THEN
        MINDEX(1,1)=1
        MINDEX(2,1)=2           ! SET LIST OF PARTICLES TO FORM MASS1
        MINDEX(3,1)=0           ! END OF LIST
        MINDEX(1,2)=3
        MINDEX(2,2)=4           ! SET LIST OF PARTICLES TO FORM MASS2
        MINDEX(3,2)=0           ! END OF LIST
        MINDEX(1,3)=0           ! SET LIST OF PARTICLES TO FORM MASS3
      ENDIF
      IF (EINDEX(1,1).LE.0.OR.EINDEX(1,1).GE.NPARTICLES) THEN
        EINDEX(1,1)=1
        EINDEX(2,1)=2           ! SET LIST OF PARTICLES TO FORM ENERGY-OBJECT1
        EINDEX(3,1)=0           ! END OF LIST
        EINDEX(1,2)=3
        EINDEX(2,2)=4           ! SET LIST OF PARTICLES TO FORM ENERGY-OBJECT2
        EINDEX(3,2)=0           ! END OF LIST
        EINDEX(1,3)=0           ! SET LIST OF PARTICLES TO FORM ENERGY-OBJECT3
      ENDIF
C
      NID=0
      IF ( ABS(ITFMN(1)).EQ.1 ) THEN
        NBF = 4
        IERR = 0
      ENDIF
      IF ( ABS(ITFMN(1)).EQ.2.OR.ABS(ITFMN(1)).EQ.4 ) THEN
        NBF = 5
        IERR = 0
      ENDIF
      IF ( ABS(ITFMN(1)).EQ.3 ) THEN
        NID = 1
        NBF = 5
        IERR = 0
      ENDIF
      IF ( ABS(ITFMN(1)).EQ.5 ) THEN
        NID = 2
        NBF = 5
        IERR = 0
      ENDIF
      IF ( ABS(ITFMN(1)).EQ.6 ) THEN
        NBF = 6
        IERR = 0
      ENDIF
      IF ( ABS(ITFMN(1)).EQ.7 ) THEN
        NID = 2
        NBF = 6
        IERR = 0
      ENDIF
      IF ( ABS(ITFMN(1)).EQ.8 ) THEN
        NBF = 5
        IERR = 0
      ENDIF
      IF ( IERR.NE.0  ) GOTO 999
C
C DO WE WANT TO SWITCH OFF SOME 4-MOMENTUM CONSTRAINTS??
C
      DO I=1,3
        IF (ABS(CVAL(I)).GT.ABS(CVAL(4))) NBF=NBF-1
      ENDDO
      IF (CVAL(4).LE.0.0) THEN
        NBF=NBF-1
      ENDIF
C
c
c setup commonblocks for b-matrix calculation
c
c Check and put momentum input vector in common 
c AND construct reference frame (if needed) for BMATRIX routine
c
      CALL DEFUVEC(IPAM,NJET,P_REC,IERR)
      IF (IERR.NE.0) GOTO 999
C
C CONSISTENCY CHECK OF VECTORS
C
      DO I=1,NJET
        IF (MIM(I).LE.0.0D0.AND.MOD(IPAM(I)/10,10).GT.1) IPAM(I)
     &       =IPAM(I)-10*MOD(IPAM(I)/10,10)
      ENDDO
C
      NDF   = NBF
      NCONS = NBF
      NFREE = 0
      NLFREE = 0
      NPARAMS = 0
      DO I=1,NJET
        NPARAM(I)=NPAR(MOD(IPAM(I),10))
        IF (MOD(IPAM(I)/10,10).GT.1) NPARAM(I)=NPARAM(I)+1
        IF (MOD(IPAM(I)/100,10).EQ.1.OR.MOD(IPAM(I)/100,10).EQ.3)
     &       NPARAM(I)=NPARAM(I)-2
        IF (I.LE.NJET-NUP.AND.(MOD(IPAM(I)/100,10).EQ.2.OR.MOD(IPAM(I)
     &       /100,10).EQ.3)) THEN
          IF (MOD(IPAM(I)/10,10).LT.1) THEN
            LFREEP(NLFREE+1)=NPARAMS+1
          ELSE
            LFREEP(NLFREE+1)=NPARAMS+1
            LFREEP(NLFREE+2)=NPARAMS+4
          ENDIF
          NLFREE=NLFREE+NPARAM(I)
          IF (MOD(IPAM(I)/100,10).EQ.2) NLFREE=NLFREE-2
        ENDIF
        IF (I.GT.NJET-NUP) NFREE=NFREE+NPARAM(I)
        NPARAMS=NPARAMS+NPARAM(I)
      ENDDO
      NFREE=NFREE+NLFREE
      NID=NID+NPARAMS
c too many unmeasured particle-variables
c      PRINT*,'IERR,NFREE,NDF',IERR,NFREE,NDF
      IF (IERR.EQ.0.AND.NFREE.GE.NDF) IERR=-2
      IF (IERR.NE.0) GOTO 999
C
C initialise constrain values
C
      DO I=1,3
        IF ( M0(I).NE.0 ) THEN
          MZ(I)    = DBLE(M0(I))
        ELSE
          MZ(I)    = 1.D-8
        ENDIF
        GZ(I)      = DBLE(G0(I))
      ENDDO
      DO I=1,4
        PSYS(I)=DBLE(CVAL(I))
      ENDDO
      PSYS(5)=MZ(1)
      PSYS(6)=MZ(2)
      PSYS(7)=MZ(3)
C
C setup y0 parameters and their errors (dy0/correlation)
C
      IBOUND=0
      I=NJET
c      IF (IFLAG.LT.0) I=ABS(IFLAG)
      IF (USER_AIBI_EVOL) THEN
        CALL AIBI_EVOL_USER(PARID,ITEVOLKMN,I
     &        ,NUP,P_COR,Y0,DY0,V,IERR)
      ELSE
        IF (ISMEAR.EQ.1)  SMOOTHVARFLAG='E'
        CALL AIBI_EVOL(PARID,ITEVOLKMN,I,NUP
     &        ,P_COR,Y0,DY0,V,IERR)
        
        IF (ISMEAR.EQ.1) THEN 
           CALL RNORML(RVEC,6)
           Y(1)=1D0+RVEC(1)*DY0(1)
           Y(2)=RVEC(1)*DY0(2)
           Y(3)=RVEC(1)*DY0(3)
           Y(4)=1D0+RVEC(4)*DY0(4)
C           PRINT*,'Y',(REAL(Y(I)),I=1,4)
C           DO I=1,NJET
C              PRINT*,'PTRUE',(P_REC(J,I),J=1,4)
C           ENDDO   
           CALL CVY2P(NJET,IPAM,Y,DPFIT)
           DO I=1,NJET
              DO J=1,4
                 P_FIT(J,I)=REAL(DPFIT(J,I))
              ENDDO
           ENDDO
C           DO I=1,NJET
C              PRINT*,'PSME',(P_FIT(J,I),J=1,4)
C           ENDDO 

           RETURN
        ENDIF   

        IF (IERR.GE.3.AND.IERR.LE.4) IBOUND=IERR-2
      ENDIF

 
c valid y0/dy0/v output???        
      IERR = 0
      IF (V(1,1).LE.0.0D0) IERR=-4
      IF ( IERR.NE.0  ) GOTO 999
c
c Remember that NUP are always fit with ALEPH
c
      IF (NFREE.GT.0) THEN
        J=0
        DO I=1,NJET
          IF (I.GT.NJET-NUP) Y0(J+1)=1.0D0
          IF (MOD(IPAM(I)/100,10).EQ.2.OR.MOD(IPAM(I)/100,10).EQ.3) Y0(J
     &         +1)=1.0D0
          J=J+NPARAM(I)
        ENDDO
      ENDIF
C
C setup alpha-parameters of constrains
C
c Breit wigner initialisation
c in principle alpha0 and sigma can be whatever
c you want, but small sigma gives less mathematical
c overflow
      ALPHA0 = 1.
      BETA0 = 1.
      SIGMA  = 0.0001
c
      IF ( ABS(ITFMN(1)).EQ.3 ) THEN
         IF ( ABS(ITFMN(2)).EQ.0 ) THEN
c Gaussian constraint
            Y0(NPARAMS+1)  = 1.D0
            DY0(NPARAMS+1) = GZ(1)/MZ(1)
         ELSE
c Breit Wigner constraint
c In principle the final result should not depend on
c the width of the gaussian parameter. So we can
c keep the same definition as for Gaussian case.
c Remarks: The smaller dy0 is, the less negative SQRT
            Y0(NPARAMS+1)  = ALPHA0
            DY0(NPARAMS+1) = SIGMA
            GAMMA(1)       = GZ(1)/MZ(1)
         ENDIF
         V(NPARAMS+1,NPARAMS+1)=DY0(NPARAMS+1)**2
      ENDIF
      IF ( ABS(ITFMN(1)).EQ.5.OR.ABS(ITFMN(1)).EQ.7 ) THEN
c Gaussian constraint
        IF ( ABS(ITFMN(2)).EQ.0 ) THEN
           Y0(NPARAMS+1)  = 1.D0
           Y0(NPARAMS+2)  = 1.D0
           DY0(NPARAMS+1) = GZ(1)/MZ(1)
           DY0(NPARAMS+2) = GZ(2)/MZ(2)
        ELSE
c Breit Wigner constraint
           Y0(NPARAMS+1)  = ALPHA0
           Y0(NPARAMS+2)  = ALPHA0
           DY0(NPARAMS+1) = SIGMA
           DY0(NPARAMS+2) = SIGMA
           GAMMA(1)       = GZ(1)/MZ(1)
           GAMMA(2)       = GZ(2)/MZ(2)
        ENDIF
        V(NPARAMS+1,NPARAMS+1)=DY0(NPARAMS+1)**2
        V(NPARAMS+2,NPARAMS+2)=DY0(NPARAMS+2)**2
      ENDIF
c
c we need to invert V if it is not diagonal
c
      IF ( ITEVOLKMN(2).NE.0 ) THEN
        CALL UCOPY(V,VINV,PLEN*PLEN*2)
        CALL DINV(NID,VINV,PLEN,WORK,IERT)
        IF ( IERT.NE.0 ) THEN
c error invertion/v is forced to be diagonal
          WRITE(6,*)
     &       'Error in ABCFITKERNEL: Failed to invert covariance matrix'
          ITEVOLKMN(2)=0
          CALL VZERO(V,PLEN*PLEN*2)
          CALL VZERO(VINV,PLEN*PLEN*2)
          DO I = 1,NID
            V(I,I) = DY0(I)**2
            IF (V(I,I).NE.0D0 ) VINV(I,I) = 1.D0/V(I,I)
          ENDDO
          IF (IERR.EQ.0.OR.IERR.GE.1) IERR=1
        ENDIF
      ENDIF
C
      IF (IFLAG.LE.0) THEN
        K=0
        DO I=1,NJET
          DO J=1,4
            DPFIT(J,I)=DBLE(P_REC(J,I))
          ENDDO
          IF (I.LT.ABS(IFLAG)) K=K+NPARAM(I)
        ENDDO
        CALL CVP2Y(NJET,IPAM,DPFIT,Z,I)
        DO I=MAX(1,ABS(IFLAG)),NJET
          DO J=1,NPARAM(I)
            Y0(K+J)=Z(K+J)
            IF (IFLAG.LT.0) THEN
              DY0(K+J)=0.5D0*10.0D0**(DBLE(MIN(2,J))-1.0D0)
              V(K+J,K+J) = DY0(K+J)**2
            ELSE
              DY0(K+J)=10.0D0*DY0(K+J)
              V(K+J,K+J) = DY0(K+J)**2
            ENDIF
          ENDDO
          K=K+NPARAM(I)
        ENDDO
        PRINT*,'Y0: ',Y0
        PRINT*,'DY0: ',DY0
        
      ENDIF
      if (iflag.gt.1) then
        PRINT*,'Y0: ',Y0
        PRINT*,'dY0: ',dY0
        print*,'ipam: ',ipam
      endif
c
c Init fitting procedure
c
      CHI2T = CHI2MAX
      CHI2MN = 1.0E10
      CHI2TX  = CHI2MN
      NIDREAL=NID
      DO I = 1,NID
        Y(I) = Y0(I)
        DO J = 1,NBF
          B(J,I) = 0D0
        ENDDO
      ENDDO
c
c Overwrite default starting point for recovery mode
c      
      IF (IFLAG.EQ.2) THEN
        DO I=1,NJET
          DO J=1,4
            DPFIT(J,I)=DBLE(P_FIT(J,I))
          ENDDO
        ENDDO
        CALL CVP2Y(NJET,IPAM,DPFIT,Y,I)
      ENDIF
c
c Start fitting procedure
c
c      PRINT*,'FX START',FX
      DO 990 NITER=1,NITERMX
       
         
 
        J=1
        DO I = 1,NJET
c
c prepare for bmatrix routines
c
          IF (I.GT.1) J=J+NPARAM(I-1)
          AI(I)   = Y(J)
C
C CHECK LOWER/UPPER BOUND ON AI
C
          IF (I.LE.NJET-NUP.AND.MOD(IPAM(I)/100,10).NE.2.AND.MOD(IPAM(I)
     &         /100,10).NE.3) THEN
            
            IF (MOD(IPAM(I),10).EQ.0) THEN
              IF (PIM(I)*ABS(AI(I)).LE.FXMAX) THEN
c too LOW value for ALEPH parameterisation
                IERR = -6
                GOTO 999
              ENDIF
              IF (PIM(I)*AI(I).GE.PMAX) THEN
c too LARGE value for ALEPH parameterisation
                IERR = -7
                GOTO 999
              ENDIF
            ELSEIF (MOD(IPAM(I),10).EQ.1) THEN
              IF (LOG(PIM(I))+AI(I).LE.LOG(FXMAX)) THEN
c too LOW value for DELPHI parameterisation
                IERR = -6
                GOTO 999
              ENDIF
              IF (LOG(PIM(I))+AI(I).GE.LOG(PMAX)) THEN
c too LARGE value for DELPHI parameterisation
                IERR = -7
                GOTO 999
              ENDIF
            ENDIF
            IF (MOD(IPAM(I)/10,10).EQ.3) THEN
              IF (LOG(MIM(I))+DI(I).GE.LOG(PMAX)) THEN
c too LARGE value for EXP-type MASS-parameterisations
                IERR = -12
                GOTO 999
              ENDIF
            ELSEIF (MOD(IPAM(I)/10,10).EQ.4) THEN
              IF (LOG(PIM(I))+DI(I).GE.LOG(PMAX)) THEN
c too LARGE value for EXP-type MASS-parameterisations
                IERR = -12
                GOTO 999
              ENDIF
            ENDIF
          ENDIF
          IF (MOD(IPAM(I)/100,10).EQ.1.OR.MOD(IPAM(I)/100,10).EQ.3) THEN
            BI(I)   = 0.0D0
            CI(I)   = 0.0D0
            IF (MOD(IPAM(I)/10,10).GT.1) THEN
              DI(I)   = Y(J+1)
            ELSE
              DI(I)   = 0.0D0
            ENDIF
          ELSE
            BI(I)   = Y(J+1)
            CI(I)   = Y(J+2)
            DI(I)   = Y(J+3)
          ENDIF
        ENDDO
        IF ( ABS(ITFMN(1)).EQ.3 ) THEN
          ALPHAIJ = Y(NPARAMS+1)
        ENDIF
        IF ( ABS(ITFMN(1)).EQ.5.OR.ABS(ITFMN(1)).EQ.7 ) THEN
          ALPHAIJ  = Y(NPARAMS+1)
          ALPHAKL  = Y(NPARAMS+2)
        ENDIF
c     
c Compute B matrix and constraint function fx(y)
c
        IF (ABS(ITFMN(2)).NE.0) THEN
          CALL BMATRIX(10+ITFMN(1),NJET,IPAM,Y)
        ELSE
          CALL BMATRIX(ITFMN(1),NJET,IPAM,Y)
        ENDIF

c        if (iflag.gt.1) then
c          print*,'niter,fx ',niter,fx,nparams
c          PRINT*,'Y0: ',(Y0(I),I=1,nparams)
c          PRINT*,'Y: ',(Y(I),I=1,nparams)
c        endif

c        PRINT*,'NITER,FX',NITER,FX
c        IF (niter.eq.2) stop
c
c compute global constrain fulfillment estimator
c
        F2=0.0D0
        DO I=1,NBF
          F2 = F2 + FX(I)**2
        ENDDO
c
c remap parameters etc. and put all unmeasured parameters at the end
c
        IF (NFREE.GE.1) THEN
          IF (NLFREE.GE.1) THEN
            DO I=NLFREE,1
              CALL REMAP(ITEVOLKMN(2)*10,NID,NBF,LFREEP(I),LFREEP(I),Y
     &             ,Y0,DY0,B,V)
            ENDDO
          ENDIF
          IF (NID.GT.NPARAMS) CALL REMAP(ITEVOLKMN(2)*10
     &         ,NID,NBF,NPARAMS-NFREE+1,NPARAMS-NLFREE,Y,Y0,DY0,B,V)
          NID=NID-NFREE
        ENDIF       
C
C now start the tedious matrix calculations
C
c Full calculation work with non diagonal matrix
c
c compute vbt = VB{T} FOR SPECIAL CASE
c
        CALL DMMLT(NID,NID,NBF,V(1,1),V(1,2),V(2,1),B(1,1),
     &      B(2,1),B(1,2),VBT(1,1),VBT(1,2),VBT(2,1),WORK)
c
c compute bvbt = B*VB{T}
c
        CALL DMMLT(NBF,NID,NBF,B(1,1),B(1,2),B(2,1),
     &      VBT(1,1),VBT(1,2),VBT(2,1),BVBT(1,1),BVBT(1,2),BVBT(2,1),
     &      WORK)
c
c invert bvb{T}
c
c      PRINT*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',NITER  
c      PRINT *,'ABC_DEBUG: bVBT'
c      DO I=1,NBF
c        PRINT '(20f12.5)',(BVBT(I,j),j=1,NBF)
c      ENDDO
c      PRINT*,'v: ',(v(i,i),i=1,nparams)
C      PRINT *,'ABC_DEBUG: B'
c      DO I=1,NBF
c        PRINT '(20f12.5)',(B(I,j),j=1,nparams)
c      ENDDO
c      PRINT*,'NBFMX',NBFMX
        CALL DINV(NBF,BVBT,NBFMX,WORK,IERT)

c        PRINT *,'ABC_DEBUG: bVBT invert',IERT
c      DO I=1,NBF
c        PRINT '(20f12.5)',(BVBT(I,j),j=1,NBF)
c      ENDDO
        IF ( IERT.EQ.-1) THEN
c  matrix could not be inverted. In some cases, one row
c  of b is zero (essentially at the first step, when
c  b{i} and c{i} are equal to zero)
          IF ( NITER.EQ.1 ) THEN
            IF (NFREE.GE.1) THEN
              NID=NIDREAL
              IF (NID.GT.NPARAMS) CALL REMAP(ITEVOLKMN(2)
     &             *10,NID,NBF,NPARAMS-NFREE+1,NID-NFREE+NLFREE,Y,Y0,DY0
     &             ,B,V)
              DO I=1,NLFREE
                CALL REMAP(ITEVOLKMN(2)
     &             *10,NID,NBF,LFREEP(I),NID-1,Y,Y0,DY0,B,V)
              ENDDO
            ENDIF
            J=1
            DO I = 1,NJET
              IF (I.GT.1) J=J+NPARAM(I-1)
              DO K=1,NPARAM(I)-1
                Y(J+K)=DY0(J+K)
              ENDDO
            ENDDO
            GOTO 990
          ENDIF
          IERR = -8
          GOTO 999
        ENDIF
c        PRINT*,'FX END',FX
c
c compute VB{T}(BVB{T})-1
c
        CALL DMMLT(NID,NBF,NBF,VBT(1,1),VBT(1,2),VBT(2,1),BVBT(1,1),
     &      BVBT(1,2),BVBT(2,1),VBT(1,1),VBT(1,2),VBT(2,1),WORK) 
c      PRINT *,'ABC_DEBUG: VBT'
c      DO I=1,NBF
c        PRINT '(20f12.5)',(VBT(j,i),j=1,NID)
c      ENDDO

        IF (NFREE.GE.1) THEN
c
c and some other things for the general case
c
c compute WbA
c      
          CALL DMMLT(NBF,NBF,NFREE,BVBT(1,1),BVBT(1,2)
     &          ,BVBT(2,1),B(1,NID+1),B(1,NID+2),B(2,NID+1),WA(1,1),WA(1
     &          ,2),WA(2,1),WORK)
c
c compute A{T} WA for general case
c
          CALL DMMLT(NFREE,NBF,NFREE,B(1
     &         ,NID+1),B(2,NID+1),B(1,NID+2),WA(1,1),WA(1,2),WA(2,1)
     &         ,WA(1,1),WA(1,2),WA(2,1),WORK)
           
          CALL DINV(NFREE,WA,PLEN,WORK,IERT)
          IF ( IERT.EQ.-1) THEN
c matrix could not be inverted. In some case, one raw
c of b is zero (essentially at the first step, when
c b{i} and c{i} are equal to zero)
            IF ( NITER.EQ.1 ) THEN
              NID=NIDREAL
              IF (NID.GT.NPARAMS) CALL REMAP(ITEVOLKMN(2)
     &             *10,NID,NBF,NPARAMS-NFREE+1,NID-NFREE+NLFREE,Y,Y0,DY0
     &             ,B,V)
              DO I=1,NLFREE
                CALL REMAP(ITEVOLKMN(2)
     &             *10,NID,NBF,LFREEP(I),NID-1,Y,Y0,DY0,B,V)
              ENDDO
              J=1
              DO I = 1,NJET
                IF (I.GT.1) J=J+NPARAM(I-1)
                DO K=1,NPARAM(I)-1
                  Y(J+K)=DY0(J+K)
                ENDDO
              ENDDO
              GOTO 990
            ENDIF
            IERR = -9
            GOTO 999
          ENDIF
C
C COMPUTE A{T} WB
C
          CALL DMMLT(NFREE,NBF,NBF,B(1,NID+1),B(2,NID+1)
     &         ,B(1,NID+2),BVBT(1,1),BVBT(1,2),BVBT(2,1),WORK2(1,1)
     &         ,WORK2(1,2),WORK2(2,1),WORK)
C
C COMPUTE WA A{T} WB (SAVED IN WAWB)
C
          CALL DMMLT(NFREE,NFREE,NBF
     &         ,WA(1,1),WA(1,2),WA(2,1),WORK2(1,1),WORK2(1,2),WORK2(2,1)
     &         ,WAWB(1,1),WAWB(1,2),WAWB(2,1),WORK)
C
C CALCULATE NEW VBT FOR GENERAL CASE (= VBT*(1-A*WAWB))
C
          CALL DMMLT(NBF,NFREE,NBF,B(1,NID+1),B(1,NID+2)
     &         ,B(2,NID+1),WAWB(1,1),WAWB(1,2),WAWB(2,1),WORK2(1,1)
     &         ,WORK2(1,2),WORK2(2,1),WORK)
          CALL DMMLT(NID,NBF,NBF,VBT(1,1),VBT(1,2),VBT(2,1),WORK2(1,1),
     &         WORK2(1,2),WORK2(2,1),WORK2(1,1),WORK2(1,2),WORK2(2,1),
     &         WORK)
          CALL DMSUB(NID,NBF,VBT(1,1),VBT(1,2),VBT(2,1),WORK2(1,1),
     &         WORK2(1,2),WORK2(2,1),VBT(1,1),VBT(1,2),VBT(2,1))
        ENDIF
c
c compute fx(y{l}) = B(y{l}-y{0})-f{y{l}}
c
        DO I = 1,NIDREAL
          YMY0(I) = Y(I)-Y0(I)
          Y_OLD(I)= Y(I)
          Y(I)    = Y0(I)
        ENDDO
        
        CALL DMMLS(NBF,NIDREAL,1,B(1,1),B(1,2),B(2,1),YMY0(1),YMY0(2)
     &       ,YMY0(2),FX(1),FX(2),FX(2))
c
c compute y{l+1} = y{0}+vb{T}(bvb{T})-1*fx(y{l})
c
        
        CALL DMMLA(NID,NBF,1,VBT(1,1),VBT(1,2),VBT(2,1),
     &      FX(1),FX(2),FX(2),Y(1),Y(2),Y(2))
       
        IF (NFREE.GE.1) CALL DMMLA(NFREE,NBF,1,WAWB(1,1)
     &       ,WAWB(1,2),WAWB(2,1),FX(1),FX(2),FX(2),Y(NID+1),Y(NID+2)
     &       ,Y(NID+2))
c        
c Compute chi2
c
 200    CHI2TX = 0.D0
        DO I = 1,NIDREAL
          YMY0(I) = Y(I)-Y0(I)
        ENDDO
c        PRINT*,'NITERNITERNITERNITER',NITER
c        PRINT*,'FX',FX
c        PRINT*,'Y   ',(Y(I),I=1,3)
c        PRINT*,'Y0  ',(Y0(I),I=1,3)
c        PRINT*,'YMY0',(YMY0(I),I=1,3)
c        PRINT*,'DY0 ',(DY0(I),I=1,3)
        IF ( ITEVOLKMN(2).EQ.0 ) THEN
          DO I = 1,NID
            IF ( V(I,I).NE.0 ) CHI2TX = CHI2TX + (YMY0(I))**2/V(I,I)
          ENDDO
        ELSE
c Non diagonal matrix: full calculation
c    chi2 = (y-y0)^{T} V{-1} (y-y0)
c z = V{-1}(y-y0)
          CALL DMMLT(NID,NID,1,VINV(1,1),VINV(1,2),VINV(2,1),YMY0(1),
     &      YMY0(2),YMY0(2),Z(1),Z(2),Z(2),WORK)
          DO I = 1,NID
            CHI2TX = CHI2TX + YMY0(I)*Z(I)
          ENDDO
        ENDIF
        
        DCHI2=ABS(CHI2TX-CHI2MN)
        IF (CHI2TX.NE.0.0D0) DCHI2=DCHI2/ABS(CHI2TX)
c        IF ( (DCHI2.LE.1.0D-8).OR
c     &       .NITER.EQ.NITERMX) THEN
c        print*,'chi2tx,CHI2MN,niter: ',chi2tx,CHI2MN
c     &       ,niter,dchi2,f2,FXMAX*FXMAX*DBLE(NBF*NBF)
c
c now check if convergence criteria is fulfilled
c        
        IF ( (DCHI2.LE.CHIMX.AND.F2.LE.FXMAX*FXMAX*DBLE(NBF*NBF)).OR
     &       .NITER.EQ.NITERMX) THEN

          
      
c
c accuracy reached
c fill all output parameters and leave the routine
c
          IF ( CHI2TX.GE.DBLE(CHI2MAX) ) THEN
            CHI2TX = DBLE(CHI2MAX)
            IERR = -10
            GOTO 999
          ENDIF
C
C COMPUTE ERROR MATRIX
C
          CALL DMMLT(NBF,NID,NID,B(1,1),B(1,2),B(2,1),V(1,1),V(1,2),
     &        V(2,1),WORK2(1,1),WORK2(1,2),WORK2(2,1),WORK)
          CALL DMMLT(NID,NBF,NID,VBT(1,1),VBT(1,2),VBT(2,1),WORK2(1,1),
     &        WORK2(1,2),WORK2(2,1),DERRMTR(1,1),DERRMTR(1,2),
     &        DERRMTR(2,1),WORK)
          IF (NFREE.GE.1) THEN
            CALL DMMLT(NFREE,NBF,NID,WAWB(1,1),WAWB(1,2
     &            ),WAWB(2,1),WORK2(1,1),WORK2(1,2),WORK2(2,1),WORK2(1,1
     &            ),WORK2(1,2),WORK2(2,1),WORK)
            DO I=1,NFREE
              DO J=1,NFREE
                DERRMTR(NID+I,NID+J)=WA(I,J)
              ENDDO
            ENDDO
            DO I=1,NID
              DO J=1,NFREE
                DERRMTR(NID+J,I)=-WORK2(J,I)
                DERRMTR(I,NID+J)=-WORK2(J,I)
              ENDDO
            ENDDO
          ENDIF
C
C PULL CALCULATION
C
          DO I=1,NID
            IF (DERRMTR(I,I).GT.0.0D0) THEN
              DPULL(I)=YMY0(I)/SQRT(DERRMTR(I,I))
            ELSE
              DPULL(I)=-1.0D10
            ENDIF
          ENDDO
C
C FINALISE ERROR MATRIX
C
          CALL DMSUB(NID,NID,V(1,1),V(1,2),V(2,1),DERRMTR(1,1),
     &        DERRMTR(1,2),DERRMTR(2,1),DERRMTR(1,1),DERRMTR(1,2),
     &        DERRMTR(2,1))
          IF (NFREE.GE.1) THEN
            DO I=NID+1,NIDREAL
               DPULL(I)=0.0D0
            ENDDO
            NID=NIDREAL
            IF (NID.GT.NPARAMS) THEN
              CALL REMAP(ITEVOLKMN(2)*10,NID,NBF,NPARAMS-NFREE+1,NID
     &             -NFREE+NLFREE,Y,Y0,DY0,B,V)
              CALL REMAP(-10,NID,NBF,NPARAMS-NFREE+1,NID
     &             -NFREE+NLFREE,DPULL,Y0,DY0,B,DERRMTR)
            ENDIF
            DO I=1,NLFREE
              CALL REMAP(ITEVOLKMN(2)*10,NID,NBF,LFREEP(I),NID-1,Y,Y0
     &             ,DY0,B,V)
              CALL REMAP(ITEVOLKMN(2)*10,NID,NBF,LFREEP(I),NID-1,DPULL
     &             ,Y0,DY0,B,DERRMTR)
            ENDDO
          ENDIF
C
C DONE -- WHAT ELSE DO YOU DREAM OFF???
C
c compute fitted jets
c
         
          CALL CVY2P(NJET,IPAM,Y,DPFIT)

          DO I = 1,NJET
c
c save output fitted parameters
c
            IF ( MOD(IPAM(I),10).EQ.0.AND.AI(I).LT.0.0D0.AND.I.LE.NJET
     &           -NUP.AND.(IERR.EQ.0.OR.IERR.GE.4) ) IERR = 4
            IF ( MOD(IPAM(I),10).EQ.0.AND.AI(I).LT.0.0D0.AND.I.GT.NJET
     &           -NUP.AND.IERR.EQ.4 ) IERR = 3
            IF ( MOD(IPAM(I),10).EQ.0.AND.AI(I).LT.0.0D0.AND.I.GT.NJET
     &           -NUP.AND.(IERR.GE.5.OR.IERR.EQ.0) ) IERR = 5
            IF ( IPAM(I).NE.PARID(I).AND.(IERR.GE.6.OR.IERR.EQ.0).AND.I
     &           .LE.NJET-NUP ) IERR = 6
            DO J = 1,4
              WORK(J) = DPFIT(J,I)
            ENDDO
            WORK(5) = SQRT(WORK(1)**2+WORK(2)**2+WORK(3)**2)
            IF ( WORK(4).LT.WORK(5) ) DPFIT(4,I) = WORK(5)
          ENDDO
          IF ( NITER.EQ.NITERMX.AND.(IERR.EQ.0.OR.IERR.GE.2)) IERR = 2
          IF (IERR.EQ.2.AND.F2.GT.FXMAX*FXMAX*DBLE(NBF*NBF)) IERR = -11
c          print*,'niter,fx FINAL',niter,fx,nid,nparams,nfree,nup
c     &         ,lfreep(1)
c        PRINT*,'Y0: ',Y0
c        PRINT*,'Y: ',Y
          GOTO 999
        ENDIF
c continue fitting procedure
C
C HAS THE FIT GOT DIFFICULTIES TO CONVERGE??
C THEN LET'S TRY SOME BLACK MAGIC
C
c        IF (NITER.EQ.30.and.iflag.gt.1) THEN
c          DO I=1,NID
c            Y(I) = (1.0D0-0.05D0)*Y(I)+0.05D0
c     &           *Y_OLD(I)
c          ENDDO
c        ENDIF
C GOT IT???
        IF ( CHI2TX.LT.0D0 ) THEN
          PRINT *,'ABCFIT ERROR: CHI2 IS NEGATIV !!!'
          CHI2MN = -CHI2TX
        ELSE
          CHI2MN = CHI2TX
        ENDIF
c
c now remap unmeasured particles back to their old position
c
        IF (NFREE.GE.1) THEN
          NID=NIDREAL
          IF (NID.GT.NPARAMS) CALL REMAP(ITEVOLKMN(2)
     &         *10,NID,NBF,NPARAMS-NFREE+1,NID-NFREE+NLFREE,Y,Y0,DY0
     &         ,B,V)
          DO I=1,NLFREE
            CALL REMAP(ITEVOLKMN(2)
     &           *10,NID,NBF,LFREEP(I),NID-1,Y,Y0,DY0,B,V)
          ENDDO
        ENDIF
 990  ENDDO
 999  CONTINUE
      NITERR = NITER
      NITMK  = NITER
      NDF    = NBF
      IF (NFREE.GE.1) NDF=NDF-NFREE
      IF ( IERR.GE.0 ) THEN
c------------------only real output---------
        IF (IERR.EQ.0) THEN
          IF (IBOUND.EQ.1) IERR=8
          IF (IBOUND.EQ.2) IERR=9
        ENDIF
        CHI2T = ABS(REAL(CHI2TX))
        DO I=1,NID
          DPARAM(I)=Y(I)
        ENDDO
        DO I=1,NJET
          DO J=1,4
            P_FIT(J,I)=REAL(DPFIT(J,I))
          ENDDO
        ENDDO
        NUMJET = NJET
      ELSE
        CHI2T = CHI2MAX
      ENDIF

      RETURN
      END
C
C SPECIAL RECOVER ROUTINE
C
      SUBROUTINE ABCFITRECOVER(NJET,NUP,P_REC,PARID,ITF,CVAL,M0,G0,P_COR
     &     ,ITEVOL,P_FIT,CHI2T,NDF,IERR)
*********************************************************************
*
* THE STRATEGY FOR RECOVERY IS TO GET AN ALTERNATIVE STARTING POINT
* (NOT USED IN CHI2 DEF OR AS INPUT FOR EVOLUTION -- ONLY AS STARTING
* POINT => DOES NOT BIAS THE RESULT)
* 
* NORMALLY THE STARTING POINT OF ABCFIT IS CHI2=0 WITH F2 > 0 (I.E.
* CONSTRAINS ARE NOT SATISFIED). WITH THE ALTERNATIVE STARTING POINT IT
* IS OPPOSITE: F2 = 0 AND CHI2 > 0. THE ALTERNATIVE STARTING IS
* DETERMINED BY AN EXPLICIT SOLUTION OF A SIMPLIFIED SYSTEM.
*
*********************************************************************
      IMPLICIT NONE

#include "abcfit_setup.inc"
#include "abcfit.inc"
#include "abcfit_bmatrix.inc"

C INPUT PARAMETERS
      INTEGER NJET,NUP,ITF,PARID(NPARTICLES),ITEVOL
      REAL    P_REC(4,NPARTICLES),CVAL(4),M0(3),G0(3),P_COR(4,NPARTICLES
     &     )
C OUTPUT PARAMETERS
      INTEGER NDF,IERR
      REAL    P_FIT(4,NPARTICLES),CHI2T
C
C INTERNAL
C
      REAL CVAL_SAVE(4)
      INTEGER I,J,IG0
      INTEGER PARID_SAVE(NPARTICLES),NUP_SAVE,NJET_SAVE
     &     ,NSPLIT(NPARTICLES),ISPLIT,IERR_SAVE,IOPTION,IOPTMAX,ITF_SAVE
     &     ,NBF,IDEBUG

      PARAMETER (IOPTMAX=11,IDEBUG=1)

      INTEGER STATUS(IOPTMAX+1)

      COMMON/FILL/FILL
      REAL FILL(10)

C
c      PRINT*,'*********************************************'
c      PRINT*,'*      TRYING RECOVERY OF A NASTY EVENT     *'
c      PRINT*,'*      ERROR-CODE ON ENTRY: ',IERR,'    *'
c      PRINT*,'*********************************************'
      RETURN
C
C    SAVE THE INPUT INFO
C
      NUP_SAVE=NUP
      NJET_SAVE=NJET
      ITF_SAVE=ITF
      IERR_SAVE=IERR
      DO I=1,4
        CVAL_SAVE(I)=CVAL(I)
      ENDDO
      J=0      
      DO I=1,NJET
        PARID_SAVE(I)=PARID(I)
        IF (I.GT.NJET-NUP) J=J+NPARAM(I)
      ENDDO
      NBF=NDF+J
      DO I=1,IOPTMAX+1
        STATUS(I)=-1000
      ENDDO
      IF (IDEBUG.GT.0) THEN
        PRINT*,'NJET,NBF,NUP,ITF:',NJET,NBF,NUP,ITF
        DO I=1,NJET
          PRINT'(1X,A,4F12.5)','PRECO PX,PY,PZ,E: ',P_REC(1,I),P_REC(2,I
     &         ),P_REC(3,I),P_REC(4,I)
        ENDDO
      ENDIF

      FILL(1)=REAL(ITF)
      FILL(2)=REAL(IERR)
C
C TRY ALL DIFFERENT RECOVERY OPTIONS ...
C
      IOPTION=1
      DO WHILE (IOPTION.LE.IOPTMAX.AND.IERR.LT.0)
        IF ((ITF.EQ.1.OR.ITF.EQ.4.OR.ITF.EQ.8).AND.IOPTION.LE.1) THEN
          IF (NBF.GT.NJET-NUP) NUP=0
C
C THIS OPTION _MUST_ BE THE FIRST...
C
C DO we need to add and extra dummy particle(s) to fulfill the
c constraints
C
          ISPLIT=0
          IF (NBF.GT.NJET) THEN
            DO I=1,NJET
              IF (NPARAM(I).GE.4) THEN
                ISPLIT=ISPLIT+1
                NSPLIT(ISPLIT)=I
              ENDIF
            ENDDO
C
C We can only add dummy particles if we have at least one particle
C which has a mass variable
C
            IF (ISPLIT.GT.0) THEN
              
              NJET=NBF
            ENDIF
          ENDIF
          IF (NBF.LE.NJET) THEN
C
C Now set ITYPP for a linear system (scaling assymption)
C
            DO I=1,NJET
              PARID(I)=100
            ENDDO
            
            CALL ABCFITKERNEL(-NJET_SAVE,NJET,NUP,P_REC,PARID,ITF,CVAL
     &           ,M0,G0,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)

            STATUS(IOPTION)=IERR
            
            IOPTION = 2

            DO I=1,NJET
              PARID(I)=PARID_SAVE(I)
            ENDDO

            IF (ISPLIT.GT.0) THEN
              
              NJET=NJET_SAVE
            ENDIF
          ELSE
            IOPTION = 2
          ENDIF
          IF (NBF.GT.NJET-NUP_SAVE) NUP=NUP_SAVE
        ELSEIF (IOPTION.LE.1) THEN
          IOPTION = 2
        ENDIF
C
C DEFAULT EXPECTATION
C        
        IF (IERR.LT.0.AND.IOPTION.EQ.2) THEN
          I=ITEVOL
          ITEVOL=MOD(ITEVOL/10,10)*10+MOD(ITEVOL,10)
          CALL ABCFITKERNEL(0,NJET,NUP,P_REC,PARID,ITF,CVAL,M0
     &         ,G0,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)
          ITEVOL=I
          STATUS(IOPTION)=IERR
          IOPTION=3
        ENDIF
C
C NO NUPS IF ANY!
C        
        IF (IERR.LT.0.AND.IOPTION.EQ.3) THEN
          IF (NUP.GT.0) THEN
            DO I=NJET-NUP+1,NJET
              PARID(I)=001
            ENDDO
            NUP=0
            CALL ABCFITKERNEL(-NJET_SAVE+NUP_SAVE,NJET,NUP,P_REC,PARID
     &           ,ITF,CVAL,M0,G0,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)
            NUP=NUP_SAVE
            DO I=NJET-NUP+1,NJET
              PARID(I)=PARID_SAVE(I)
            ENDDO
            STATUS(IOPTION)=IERR
            IOPTION=4
          ELSE
            IOPTION=4
          ENDIF
        ENDIF
C
C FLIPPING BETWEEN EQUAL ENERGY AND EQUAL MASS
C        
        IF (IERR.LT.0.AND.IOPTION.EQ.4) THEN
          IF (ITF.EQ.4.OR.ITF.EQ.8) THEN
            ITF=MAX(4,MOD(ITF_SAVE*2,16))
            CALL ABCFITKERNEL(1,NJET,NUP,P_REC,PARID,ITF
     &           ,CVAL,M0,G0,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)
            STATUS(IOPTION)=IERR
            IOPTION=5
            ITF=ITF_SAVE
          ELSE
            IOPTION=5
          ENDIF
        ENDIF
C
C SEMI-VOODOONIC
C

        IF (IERR.LT.0.AND.IOPTION.EQ.5) THEN
          IF (ITF.GE.5.AND.ITF.LE.7) THEN
            ITF=4
            CALL ABCFITKERNEL(1,NJET,NUP,P_REC,PARID,ITF
     &           ,CVAL,M0,G0,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)
            ITF=ITF_SAVE
            STATUS(IOPTION)=IERR
            IOPTION=6
          ELSE
            IOPTION=6
          ENDIF
        ENDIF
C
C FULLY VOODOONIC
C
        IF (IERR.LT.0.AND.IOPTION.EQ.6) THEN
          IF (ITF.NE.1) THEN
            ITF=1
            CALL ABCFITKERNEL(1,NJET,NUP,P_REC,PARID,ITF
     &           ,CVAL,M0,G0,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)
            ITF=ITF_SAVE
            STATUS(IOPTION)=IERR
            IOPTION=7
          ELSE
            IOPTION=7
          ENDIF
        ENDIF
C
C SUPER VOODOONIC
C
        IF (IERR.LT.0.AND.IOPTION.EQ.7) THEN
          IF (ITF.EQ.1.OR.ITF.EQ.4.OR.ITF.EQ.8) THEN
            G0(1)=0.0
            IF (ITF.EQ.1) G0(1)=30.0
            G0(2)=G0(1)
            ITF=5
            DO IG0=1,10
              
              G0(1)=G0(1)+REAL(SIGN(1,ITF-3)*IG0)*0.1250
              G0(2)=G0(1)
              IF (G0(1).LT.0.0) GOTO 1234
              CALL ABCFITKERNEL(1,NJET,NUP,P_REC,PARID,ITF
     &             ,CVAL,M0,G0,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)
              PRINT*,'IERR,G0',IERR,G0(1)
              IF (IERR.GE.0) GOTO 1234
            ENDDO
 1234       CONTINUE
            ITF=ITF_SAVE
            STATUS(IOPTION)=IERR
            IOPTION=8
          ELSE
            IOPTION=8
          ENDIF
        ENDIF
C        
C MEGA VOODOONIC
C        
        IF (IERR.LT.0.AND.IOPTION.EQ.8) THEN
          
          DO I=4,3,-1            
            IF (I.LE.3.AND.ABS(CVAL(I)).LE.ABS(CVAL(4))) CVAL(I)=
     &             -CVAL(4)-1.
            IF (I.EQ.4.AND.CVAL(I).GT.0) CVAL(I)=-1.
            
            PRINT*,'CVAL',I,CVAL
            CALL ABCFITKERNEL(1,NJET,NUP,P_REC,PARID,ITF
     &           ,CVAL,M0,G0,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)
            STATUS(IOPTION)=IERR
            CVAL(I)=CVAL_SAVE(I)
            IOPTION=9
            
            IF (IERR.GE.0) GOTO 12345
          ENDDO
        ENDIF
12345   CONTINUE
C
C DEFAULT EXPECTATION AND        
C
C NO NUPS IF ANY!
C        
        IF (IERR.LT.0.AND.IOPTION.EQ.9) THEN
          IF (NUP.GT.0) THEN
            DO I=NJET-NUP+1,NJET
              PARID(I)=001
            ENDDO
            I=ITEVOL
            ITEVOL=MOD(ITEVOL/10,10)*10+MOD(ITEVOL,10)
            NUP=0
            CALL ABCFITKERNEL(0,NJET,NUP,P_REC,PARID,ITF
     &           ,CVAL,M0,G0,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)
            NUP=NUP_SAVE
            ITEVOL=I
            DO I=NJET-NUP+1,NJET
              PARID(I)=PARID_SAVE(I)
            ENDDO
            STATUS(IOPTION)=IERR
            IOPTION=10
          ELSE
            IOPTION=10
          ENDIF
          I=1
        ENDIF
C
C SEMI-SOLUTION 1
C
        IF (IERR.LT.0.AND.IOPTION.EQ.10) THEN
          DO WHILE (IERR.LT.0.AND.I.LE.NJET-NUP)
            IF (MOD(PARID(I)/10,2).NE.1) THEN
              PARID(I)=PARID(I)+10
              CALL ABCFITKERNEL(1,NJET,NUP,P_REC,PARID,ITF
     &             ,CVAL,M0,G0,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)
              PARID(I)=PARID_SAVE(I)
            ENDIF
            print*,'changed: ',i,ierr
            I=I+1
          ENDDO
          STATUS(IOPTION)=IERR
          IOPTION=11
          J=1
        ENDIF
C
C SEMI-SOLUTION 2
C
        IF (IERR.LT.0.AND.IOPTION.EQ.11) THEN
          DO WHILE (IERR.LT.0.AND.J.LE.NJET-NUP)
            DO I=1,NJET-NUP
              PARID(I)=PARID_SAVE(J)+10
            ENDDO
            PARID(J)=PARID_SAVE(J)
            print*,'parid :',parid
            CALL ABCFITKERNEL(1,NJET,NUP,P_REC,PARID,ITF
     &           ,CVAL,M0,G0,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)
            print*,'CHANGED: ',J,ierr
            J=J+1
          ENDDO
          DO I=1,NJET
            PARID(I)=PARID_SAVE(I)
          ENDDO
          STATUS(IOPTION)=IERR
          IOPTION=12
        ENDIF

        IF (IERR.GE.0) THEN
          PRINT*,'DOING RECOVERY OPTION: ',IOPTION-1,IERR
        ELSE
          PRINT*,'TRIED RECOVERY OPTION: ',IOPTION,IERR
        ENDIF

      PRINT*,'BINGO staTUS: ',STATUS
      print*,'cval,nup,parid,njet: ',cval,nup,parid,njet,j
      print*,'p_fit: ',p_fit,ioption,ioptmax,ierr
        IF (IOPTION-1.LE.IOPTMAX.AND.IERR.GE.0) THEN
          CALL ABCFITKERNEL(2,NJET,NUP,P_REC,PARID,ITF
     &       ,CVAL,M0,G0,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)
          IF (IOPTION.EQ.11.AND.IERR.LT.0) IOPTION=10
          IF (IOPTION.EQ.12.AND.IERR.LT.0) IOPTION=11
        ELSEIF (IOPTION.EQ.IOPTMAX+1) THEN
          DO I=1,NJET
            DO J=1,4
              P_FIT(J,I)=P_REC(J,I)
            ENDDO
          ENDDO
          CALL ABCFITKERNEL(2,NJET,NUP,P_REC,PARID,ITF
     &         ,CVAL,M0,G0,P_COR,ITEVOL,P_FIT,CHI2T,NDF,IERR)
          STATUS(IOPTMAX+1)=IERR
        ENDIF

      ENDDO
      PRINT*,'STATUS: ',STATUS
 999  CONTINUE
      
      
      FILL(3)=REAL(IERR)
      FILL(4)=REAL(IOPTION)
      CALL HFN(10,FILL)
      RETURN
      END
C
C ISR AFTERBURNER FIT CODE
C
      SUBROUTINE ABCFIT_ISR(NJET,NUP,P_REC,ITF,CVAL,M0,G0,ITYPP,ITEVOL
     &     ,P_FIT,CHI2T,NDF,IERR)
*********************************************************************
*
* ISR fit using ABCFIT. The code is a stand alone afterburner code
* which calls ABCFIT several times to minimise the function 
*
*               CHI^2 = -2ln(prob(CHI^2(ABCFIT)))-2ln(CLisr)
*
* The ISR is assumed to be along the z-axis.
*
*********************************************************************
      IMPLICIT NONE

#include "abcfit_setup.inc"
#include "abcfit.inc"

c Input parameters
      INTEGER NJET,NUP,ITF,ITYPP,ITEVOL
      REAL    P_REC(4,NPARTICLES),CVAL(4),M0(3),G0(3)
c Output parameters
      INTEGER NDF,IERR
      REAL    P_FIT(4,NPARTICLES),CHI2T
C
C INTERNAL
C

C DON'T WORRY IT WILL BE THERE SOON....

      RETURN
      END
C
C Covariance matrix conversion
C
      SUBROUTINE CY02PM(INDEX,NDIM,ERRMTROUT,PULLOUT,COVP,COVM)
      IMPLICIT NONE
*****************************************************
* Routine to convert ABCFIT covariance matrix into
* momentum and mass covariance matrix by error
* propagation.
*
* To be called just after a successful ABCFIT call....
*
* Input:
*
*       LISTE(7,2) = List of particles forming up to two masses
*                    e.g.  (1,2,0,0,0,0,0),(3,4,0,0,0,0,0) for WW
*                    or    (1,2,3,4,0,0,0),(0,0,0,0,0,0,0) for root_s of WW
*
* Output:
*
*       NDIM             = Number of fit parameters
*       ERRMTROUT(25,25) = Covariance matrix of fit parameters (a,b,c...)
*       PULLOUT(25)      = Pull of the parameters
*       COVP(28,28)      = Momentum covariance matrix (4-momentum)
*       COVM(2,2)        = Mass covariance matrix
*
*****************************************************
#include "abcfit_setup.inc"
#include "abcfit.inc"
#include "abcfit_bmatrix.inc"

      INTEGER INDEX(NPARTICLES,2),NDIM
      REAL PULLOUT(PLEN),ERRMTROUT(PLEN,PLEN),COVP(4*NPARTICLES,4
     &     *NPARTICLES),COVM(2,2)
C
      INTEGER I,J,K
      DOUBLE PRECISION P(4,NPARTICLES),DPDY(4*NPARTICLES+3,4*NPARTICLES
     &     +3),DMDP(NBFMX,4*NPARTICLES+3),WORK1(4*NPARTICLES,4
     &     *NPARTICLES),WORK(4*NPARTICLES)
C
C DID WE HAVE A GOOD FIT???
C     
      IF(NUMJET.LE.0) RETURN
C YES
      NDIM=0
      DO I=1,NUMJET
        NDIM=NDIM+NPARAM(I)
      ENDDO

      DO I=1,NDIM
        PULLOUT(I)=REAL(DPULL(I))
        DO J=1,NDIM
          ERRMTROUT(I,J)=REAL(DERRMTR(I,J))
        ENDDO
      ENDDO
C
      CALL RMSET(4*NUMJET,4*NUMJET,0.0,COVP(1,1),COVP(1,2),COVP(2,1))
      CALL RMSET(2,2,0.0,COVM(1,1),COVM(1,2),COVM(2,1))
      
      CALL CVY2P(NUMJET,IPAM,DPARAM,P)
      DO I=1,NUMJET
        PIR(I)=P(1,I)**2+P(2,I)**2+P(3,I)**2
        EIR(I)=P(4,I)
        MIR(I)=SQRT(MAX(0.0D0,P(4,I)**2-PIR(I)))
        PIR(I)=SQRT(PIR(I))
        BIER(I)=DLOG(MAX(1.0D-32,MAX(1.0D0,EIR(I)/PIR(I))-1.0D0))
      ENDDO
      CALL DP_DPARAM(NUMJET,IPAM,P,DPDY)
C V*DPDY(T)
      CALL DMMLT(NDIM,NDIM,NUMJET*4,DERRMTR(1,1),DERRMTR(1,2),DERRMTR(2
     &  ,1),DPDY(1,1),DPDY(2,1),DPDY(1,2),WORK1(1,1),WORK1(1,2),WORK1(2
     &  ,1),WORK)
C COVP = DPDY*[V*DPDY(T)]
      CALL DMMLT(4*NUMJET,NDIM,4*NUMJET,DPDY(1,1),DPDY(1,2),DPDY(2,1)
     &     ,WORK1(1,1),WORK1(1,2),WORK1(2,1),WORK1(1,1),WORK1(1,2)
     &     ,WORK1(2,1),WORK)
      DO I=1,4*NUMJET
        DO J=1,4*NUMJET
          COVP(I,J)=REAL(WORK1(I,J))
        ENDDO
      ENDDO
      
      CALL VZERO(DMDP,NBFMX*(4*NPARTICLES+3)*2)
      I=1
      DO WHILE (INDEX(1,I).GT.0.AND.I.LE.2)
        CALL DM_DP(NUMJET,I,INDEX(1,I),P,DMDP)       
        I=I+1
      ENDDO
      I=I-1
C COVP(12,12)*DMDP(T)
      CALL DMMLT(4*NUMJET,4*NUMJET,I,WORK1(1,1),WORK1(1,2),WORK1(2,1)
     &     ,DMDP(1,1),DMDP(2,1),DMDP(1,2),WORK1(1,1),WORK1(1,2),WORK1(2
     &     ,1),WORK)
C COVM = DMDP*[COVP(12,12)*DMDP(T)]
      CALL DMMLT(I,4*NUMJET,I,DMDP(1,1),DMDP(1,2),DMDP(2,1),WORK1(1,1),
     &     WORK1(1,2),WORK1(2,1),WORK1(1,1),WORK1(1,2),WORK1(2,1),WORK)
      DO K=1,MIN(2,I)
        DO J=1,MIN(2,I)
          COVM(K,J)=REAL(WORK1(K,J))
        ENDDO
      ENDDO
 999  RETURN
      END
C***********************************************************************
C
C   Misc. service routines for ABCFIT
C
C***********************************************************************
      SUBROUTINE ABPUF(NJET,ITF,PARID,ITEVOL,ITFMN,ITEVOLKMN)
C***********************************************************************
      IMPLICIT NONE

#include "abcfit_setup.inc"
#include "abcfit.inc"

      INTEGER NJET,ITF,PARID(NPARTICLES),ITEVOL
      INTEGER ITFMN(2),ITEVOLKMN(3)

      LOGICAL FIRST
      INTEGER IF_ITF(10,2),IF_ITYPP(10,3,NPARTICLES),IF_ITEVOL(10,3)
     &     ,ITYPPKMN(3),I,LTAG,INIT,H
      CHARACTER*10 NUMBERS,LINE

      PARAMETER (INIT=NPARTICLES*30)
      DATA FIRST /.TRUE./
      DATA IF_ITF /20*0/, IF_ITYPP /INIT*0/, IF_ITEVOL /30*0/
      DATA NUMBERS /'0123456789'/

      SAVE

      IF (FIRST) THEN
        WRITE(6,1000)
 1000   FORMAT(//
     &'+-------------------------------------------------------------+'/
     &'|                                                             |'/
     &'|    AAAA                                                     |'/
     &'|   A    A                                                    |'/
     &'|   A    A                                                    |'/
     &'|   AAAAAA                                                    |'/
     &'|   A    A                                                    |'/
     &'|   A    A                                                    |'/
     &'|   A    A                                                    |'/
     &'|                                                             |'/
     &'|            BBBB                                             |'/
     &'|            B   B                                            |'/
     &'|            B   B                                            |'/
     &'|            BBBB  etter                                      |'/
     &'|            B   B                                            |'/
     &'|            B   B                                            |'/
     &'|            BBBB                                             |'/
     &'|                                                             |'/
     &'|                     CCCC     FFFFFF   IIIIII   TTTTTTTT  !! |'/
     &'|                    C    C    F          II        TT     !! |'/
     &'|                   C          F          II        TT     !! |'/
     &'|                   C          FFFFF      II        TT     !! |'/
     &'|                   C          F          II        TT     !! |'/
     &'|                    C    C    F          II        TT        |'/
     &'|                     CCCC     F        IIIIII      TT     !! |'/
     &'|                                                             |'/
     &'+-------------------------------------------------------------+'/
     &'                                                               '/
     &'                             BY                                '/
     &'                                                               '/
     &'                  OLIVER BUCHMULLER AND #5                     '/
     &'                                                               '/
     &'                 FORMERLY KNOWN AS MATHKINE                    '/
     &'                                                               '/
     &'                  Version 1.0 (August 1998)                    '/
     &'                                                               '/
     &'                                                               '/
     &'                                                               '/
     &'      +-----------------------------------------------+'/
     &'      | Original MATHKINE package Version 1.x         |'/
     &'      |          from  Alain  BONISSENT    (CPPM)     |'/
     &'      |                Fabien MOTSCH       (CPPM)     |'/
     &'      |                Andre  TILQUIN      (CPPM)     |'/
     &'      |                Adel   TRABELSI     (SACLAY)   |'/
     &'      |                MArie-Noelle MINARD (ANNECY)   |'/
     &'      +-----------------------------------------------+'//)
        FIRST=.FALSE.
      ENDIF

      ITFMN(1) = MOD(ITF,10)
      ITFMN(2) = ITF/10
      ITEVOLKMN(1) = MOD(ITEVOL,10)
      ITEVOLKMN(2) = MOD(ITEVOL,100)/10
      ITEVOLKMN(3) = ITEVOL/100
      
      IF (IF_ITF(ITFMN(1)+1,1).EQ.0.OR.IF_ITF(ITFMN(2)+1,2).EQ.0)
     &     THEN
        IF_ITF(ITFMN(1)+1,1)=1
        IF_ITF(ITFMN(2)+1,2)=1
        WRITE(6,*) 'ABCFIT called with Energy-momentum conservation'
        IF ( ITFMN(1).EQ.2 ) WRITE(6,*)
     &       'and with the constraint m_12 =m0'
        IF ( ITFMN(1).EQ.3 ) WRITE(6,*)
     &       'and with the constraint alpha_12*m_12 = m0'
        IF ( ITFMN(1).EQ.4 ) WRITE(6,*)
     &       'and with the constraint m_12=m_34'
        IF ( ITFMN(1).EQ.5 ) WRITE(6,*) 'and with the constraint'
     &       //' alpha_12*m_12 = alpha_34*m_34'
        IF ( ITFMN(1).EQ.6 ) WRITE(6,*) 'and with the constraints'
     &       //' m_12=m01   m_34=m02'
        IF ( ITFMN(1).EQ.7 ) WRITE(6,*) 'and with the constraints'
     &       //' alpha_12*m_12 = m01   alpha_34*m3_34=m02'
        IF ( ITFMN(1).EQ.8 ) WRITE(6,*)
     &       'and with the constraint e_12=e_34'
        IF ( ITFMN(1).EQ.3.OR.ITFMN(1).EQ.5.OR.ITFMN(1).EQ.7 ) THEN
          IF ( ITFMN(2).EQ.1 ) THEN
            WRITE(6,*)
     &            '        Using Breit Wigner probability distribution'
          ELSE
            WRITE(6,*)
     &            '        Using Gaussian probability distribution'
          ENDIF
        ENDIF
      ENDIF

      DO I=1,NJET
        ITYPPKMN(1) = MOD(PARID(I),10)
        ITYPPKMN(2) = MOD(PARID(I)/10,10)
        ITYPPKMN(3) = MOD(PARID(I)/100,10)
        IF (IF_ITYPP(ITYPPKMN(1)+1,1,I).EQ.0.OR.IF_ITYPP(ITYPPKMN(2)+1,2
     &    ,I).EQ.0.OR.IF_ITYPP(ITYPPKMN(3)+1,3,I).EQ.0) THEN
          IF_ITYPP(ITYPPKMN(1)+1,1,I)=1
          IF_ITYPP(ITYPPKMN(2)+1,2,I)=1
          IF_ITYPP(ITYPPKMN(3)+1,3,I)=1
          WRITE(6,*) ' Information for particle: ',I
          IF ( ITYPPKMN(1).EQ.0 ) THEN
            WRITE(6,*) ' Fit parameterisation a la ALEPH'
          ELSEIF ( ITYPPKMN(1).EQ.1 ) THEN
            WRITE(6,*) ' Fit parameterisation a la DELPHI'
          ELSEIF ( ITYPPKMN(1).EQ.2 ) THEN
            WRITE(6,*) ' Fit parameterisation using P,theta,phi'
          ELSEIF ( ITYPPKMN(1).EQ.3 ) THEN
            WRITE(6,*) ' Fit parameterisation using Px,Py,Pz'
          ELSEIF ( ITYPPKMN(1).EQ.4 ) THEN ! JUST A HACK
            WRITE(6,*) ' YOU ARE USING A HACK parameterisation!!!!!!!!'  
          ELSE            
            WRITE(6,*) ' Undefined fit parameterisation requested'
          ENDIF

          IF ( ITYPPKMN(3).EQ.0 ) THEN
            WRITE(6,*)
     &           ' With full momentum corrections'
          ELSEIF ( ITYPPKMN(3).EQ.1 ) THEN
            WRITE(6,*)
     &           ' Only longitudinal momentum corrections are allowed'
          ELSEIF ( ITYPPKMN(3).EQ.2 ) THEN
            WRITE(6,*)
     &           'Full correction on particle momentum, but the',
     &           'longitudinal parameter(s) is assumed to be unmeasured'
          ELSEIF ( ITYPPKMN(3).EQ.3 ) THEN
            WRITE(6,*)
     &           'Longitudinal correction on paricle momentum, and the',
     &           'longitudinal parameter(s) is assumed to be unmeasured'
          ELSE
            WRITE(6,*)
     &           ' Undefined ITYPP-k value specified'
          ENDIF

          IF ( ITYPPKMN(2).EQ.0 ) THEN
            WRITE(6,*)
     &           ' Fitted jet energies scale with fitted momenta'
          ELSEIF ( ITYPPKMN(2).EQ.1 ) THEN
            WRITE(6,*)
     &           ' Fitted jet energies is given by fixed input jet '
     &           //'mass and fitted momenta'
          ELSEIF ( ITYPPKMN(2).EQ.2 ) THEN
            WRITE(6,*) ' Fit jet mass using a 4th parameter: '
     &           //'M_fit = d*M_reco'
          ELSEIF ( ITYPPKMN(2).EQ.3 ) THEN
            WRITE(6,*) ' Fit jet mass using a 4th parameter: '
     &           //'M_fit = exp(d)*M_reco'
          ELSEIF ( ITYPPKMN(2).EQ.4 ) THEN
            WRITE(6,*) ' Fit jet mass using a 4th parameter: '
     &           //'E_fit = (exp(d+log(E_reco/P_reco+1))+1)*P_fit'
          ELSE
            WRITE(6,*) ' Undefined ITYPP-m value specified'
          ENDIF
        ENDIF
      ENDDO

      IF (IF_ITEVOL(ITEVOLKMN(1)+1,1).EQ.0.OR.IF_ITEVOL(ITEVOLKMN(2)+1,2
     &     ).EQ.0.OR.IF_ITEVOL(ITEVOLKMN(3)+1,3).EQ.0) THEN
        IF_ITEVOL(ITEVOLKMN(1)+1,1)=1
        IF_ITEVOL(ITEVOLKMN(2)+1,2)=1
        IF_ITEVOL(ITEVOLKMN(3)+1,3)=1
        IF (USER_AIBI_EVOL) THEN
          WRITE(6,*) 'Correction-factors are taken from subroutine:'
     &      //' AIBI_EVOL_USER'
        ELSE
          LINE=' '
          H=1
          DO I=1,NJET
            IF (MOD(PARID(I)/10,10).GT.1) THEN
              LINE=LINE(1:H)//NUMBERS(MOD(PARID(I)/10,10)
     &             +1:MOD(PARID(I)/10,10)+1)
              H=H+1
            ENDIF
            LINE=LINE(1:H)//NUMBERS(MOD(PARID(I),10)
     &           +1:MOD(PARID(I),10)+1)
            H=H+1
          ENDDO
          LTAG=MAX(0,INDEX(ABCTAG,' ')-1)
          IF (LTAG.GT.0) THEN
            WRITE(6,'(1x,a,i1,a,i1,a)')
     &        ' Correction-factors are taken from the file aibi_evol_'
     &        //ABCTAG(1:LTAG)//'_'//LINE(2:H)//'_',ITEVOLKMN(3),'0'
     &        ,ITEVOLKMN(1),'.dat'
          ELSE
            WRITE(6,'(1x,a,i1,a,i1,a)')
     &        ' Correction-factors are taken from the file aibi_evol_'
     &        //LINE(2:H)//'_',ITEVOLKMN(3),'0',ITEVOLKMN(1),'.dat'
          ENDIF
        ENDIF
        IF ( ITEVOLKMN(2).EQ.0 ) THEN
          WRITE(6,*)
     &          ' Correction-factors have diagonal covariance matrix'
        ELSE
          WRITE(6,*) ' Correction-factors have full non-diagonal'
     &          //' covariance matrix'
        ENDIF
        IF ( ITEVOLKMN(3).EQ.0 ) THEN
          WRITE(6,*)
     &          ' Correction-factors are binned in reco-values'
        ELSE
          WRITE(6,*)
     &          ' Correction-factors are binned in true-values'
        ENDIF
      ENDIF
     
      RETURN
      END
C***********************************************************************
      SUBROUTINE REMAP(IFLAG,YDIM,FDIM,FROM,TO,Y,Y0,DY0,B,V)
C***********************************************************************
      IMPLICIT NONE

#include "abcfit_setup.inc"
#include "abcfit.inc"

      INTEGER IFLAG,YDIM,FDIM,FROM,TO
      DOUBLE PRECISION Y(PLEN),Y0(PLEN),DY0(PLEN),B(NBFMX,PLEN),V(PLEN
     &     ,PLEN)

      INTEGER I,J,K
      DOUBLE PRECISION WORK(PLEN,PLEN)
C      RETURN
C      PRINT *,'ABC_DEBUG BEFORE: B'
C      DO I=1,FDIM
C        PRINT '(20f12.5)',(B(I,j),j=1,YDIM)
C      ENDDO
C      PRINT '(20f12.5)',(Y(j),j=1,YDIM)
      DO I=FROM,TO
        WORK(FDIM+1,I)=Y(I)
        IF (IFLAG.GE.0) THEN
          WORK(FDIM+2,I)=Y0(I)
          WORK(FDIM+3,I)=DY0(I)
          IF (ABS(IFLAG).LT.10) WORK(FDIM+4,I)=V(I,I)
          DO J=1,FDIM
            WORK(J,I)=B(J,I)
          ENDDO
        ENDIF
      ENDDO

      DO I=1,YDIM-TO
        Y(FROM+I-1)=Y(TO+I)
        IF (IFLAG.GE.0) THEN
          Y0(FROM+I-1)=Y0(TO+I)
          DY0(FROM+I-1)=DY0(TO+I)
          IF (ABS(IFLAG).LT.10) V(FROM+I-1,FROM+I-1)=V(TO+I,TO+I)
          DO J=1,FDIM
            B(J,FROM+I-1)=B(J,TO+I)
          ENDDO
        ENDIF
      ENDDO

      DO I=FROM,TO
        K=I+YDIM-TO
        Y(K)=WORK(FDIM+1,I)
        IF (IFLAG.GE.0) THEN
          Y0(K)=WORK(FDIM+2,I)
          DY0(K)=WORK(FDIM+3,I)
          IF (ABS(IFLAG).LT.10) V(K,K)=WORK(FDIM+4,I)
          DO J=1,FDIM
            B(J,K)=WORK(J,I)
          ENDDO
        ENDIF
      ENDDO
C
C NON DIAGONAL V MATRIX
C
      IF (ABS(IFLAG).GE.10) THEN
        DO I=FROM,YDIM
           DO J=1,YDIM
             WORK(J,I)=V(J,I)
           ENDDO
        ENDDO

        DO I=1,MAX(YDIM-TO,TO-FROM+1)
          DO J=1,FROM-1
            IF (I.LE.YDIM-TO) V(J,FROM+I-1)=WORK(J,TO+I)
            IF (I.LE.TO-FROM+1) V(J,FROM+YDIM-TO+I-1)=WORK(J,FROM+I-1)
          ENDDO
        ENDDO

        DO I=1,MAX(TO-FROM+1,YDIM-TO)
          DO J=1,I
            IF (I.LE.YDIM-TO) V(FROM+J-1,FROM+I-1)=WORK(TO+J,TO+I)
            IF (I.LE.TO-FROM+1) V(FROM+YDIM-TO+J-1,FROM+YDIM-TO+I-1)
     &        =WORK(FROM+J-1,FROM+I-1)
          ENDDO
        ENDDO

        DO I=1,TO-FROM+1
          DO J=1,YDIM-TO
            V(FROM+J-1,FROM+YDIM-TO+I-1)=WORK(TO+J,FROM+I-1)
          ENDDO
        ENDDO

        DO I=1,YDIM
          DO J=1,YDIM
            V(J,I)=V(I,J)
          ENDDO
        ENDDO
      ENDIF
C      PRINT *,'ABC_DEBUG AFTER: B'
C      DO I=1,FDIM
C        PRINT '(20f12.5)',(B(I,j),j=1,YDIM)
C      ENDDO
C      PRINT '(20f12.5)',(Y(j),j=1,YDIM)
      
 999  RETURN
      END












