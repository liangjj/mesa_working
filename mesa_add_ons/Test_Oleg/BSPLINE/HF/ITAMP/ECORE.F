*    ==================================================================
*      E C O R E
*    ==================================================================
*
      SUBROUTINE ECORE(EL,EC,REL)
*
* ..  COMPUTES THE ENERGY OF THE COMMON CLOSED SHELLS
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
*
      EC = 0.d0
      DO 10 I = 1,NCLOSD
         SUMI = 4*L(I)+2
         TI   = FKy(I,I,0,REL)
         DO 20 K = 2,2*L(I),2
            TI = TI - CA(L(I),K)*FKy(I,I,K,REL)
   20    CONTINUE
         EC = EC + SUMI*((SUMI-1)*TI - HL(EL,I,I,REL))/2.d0
         DO 30 J = 1,I-1
            SUMJ = 4*L(J)+2
            TIJ = FKy(I,J,0,REL)
            DO 40 K=IABS(L(I)-L(J)),L(I)+L(J),2
               TIJ = TIJ -CB(L(I),L(J),K)*GKy(I,J,K,REL)
   40       CONTINUE
            EC = EC + SUMI*SUMJ*TIJ
   30    CONTINUE
   10 CONTINUE
      END
