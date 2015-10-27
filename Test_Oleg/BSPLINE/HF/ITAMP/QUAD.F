*
*     ------------------------------------------------------------------
*              Q U A D
*     ------------------------------------------------------------------
*
*       Returns the value of <a, S b> where a and b are spline
*     expansion coefficients
*
      DOUBLE PRECISION  FUNCTION quad(ns,ks,a,b)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,nwd=10)
      DOUBLE  PRECISION a(nt), b(nt)
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rk(nt,nt,kx,kx),lh,krk
*
      quad = 0.d0
      do 100 i = 2,ns
	quad = quad + a(i)*b(i)*sb(i,ks)
  100 continue
      do 120 m = 2,ks
*       ..i-m+1  < m
	do 120 i = m+1,ns
	  quad = quad +(a(i)*b(i-m+1)+a(i-m+1)*b(i))*sb(i,ks-m+1)
  120 continue
      print *, 'quad', quad
*
      END
