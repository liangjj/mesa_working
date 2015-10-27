*     =================================================================
*       H L C
*     =================================================================
*
*
      DOUBLE PRECISION  FUNCTION HLC(EL,I,J,REL)
*
* ... COMPUTES HL(I,J) MODIFIED BY THE INTERACTIONS WITH THE CLOSED
*     SHELL
*
      IMPLICIT DOUBLE PRECISION (a-h,o,p,r-z)
      PARAMETER (nt=100,kx=8,nwd=10)
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      CHARACTER EL(*)*3
      LOGICAL REL
      HLC = HL(EL,I,J,REL)
      print *, ' hlc =',hlc
      DO 10 IP = 1,NCLOSD
         SUMIP = 4*L(IP)+2
         TM = RKy(I,IP,J,IP,0,REL)
         DO 20 K = IABS(L(I)-L(IP)),L(I)+L(IP),2
            TM = TM - CB(L(I),L(IP),K)*RKy(I,IP,IP,J,K,REL)
   20    CONTINUE
         HLC = HLC - 2.d0*SUMIP*TM
   10 CONTINUE
      END
