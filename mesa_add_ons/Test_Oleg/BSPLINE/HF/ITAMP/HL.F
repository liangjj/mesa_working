*
*     ------------------------------------------------------------------
*              H L
*     ------------------------------------------------------------------
*
*       Returns the value of <i|L|j> .. assumes symmetry 
*
      DOUBLE PRECISION  FUNCTION HL(EL,i,j,REL)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,nwd=10)
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rk(nt,nt,kx,kx),lh,krk
*
      PARAMETER(IWRITE=6)
      CHARACTER EL(*)*3
      LOGICAL REL
*
      IF (IABS(L(i)-L(j)) .NE. 0) GO TO 999 
      if (lh .ne. l(i)) then
         call hlm(z,l(i),nt,kx,ks,ns,db2,rm1,rm2,h0)
	 lh = l(i)
      end if
      hl = 0.d0
      last = min0(max(i),max(j))
      do 100 k = 2,last
	hl = hl + p(k,i)*p(k,j)*h0(k,ks)
  100 continue
      if (i .eq. j) then
        do 110 kp = 2,ks
*         .. k'<k
	do   110 k = kp+1,last
	     hl = hl + 2.d0*p(k,i)*p(k-kp+1,j)*h0(k,ks-kp+1)
  110   continue
       else
        do 120 kp = 2,ks
*         .. k'<k
	  do 120 k = kp+1,last
	     hl = hl + 
     :         (p(k,i)*p(k-kp+1,j)+p(k-kp+1,i)*p(k,j))*h0(k,ks-kp+1)
  120   continue
      end if
*
      IF (REL) HL=HL-2.d0*RLSHFT(i,j)
      RETURN
  999 continue
      WRITE(IWRITE,4) EL(i),L(i),EL(j),L(j)
4     FORMAT(10X,'UNALLOWED L VALUES OCCURRED IN HL SUBROUTINE'/
     :   2(10X,A3,' HAS L = ',I3))
      stop
      END
