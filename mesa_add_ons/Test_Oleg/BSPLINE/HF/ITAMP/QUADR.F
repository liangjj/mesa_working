*     ==================================================================
*              Q U A D R
*     ==================================================================
*
*                                   m
*       Evaluates the integral of  r  P (r) P (r) with respect to r
*                                       i     j
*
      DOUBLE PRECISION  FUNCTION QUADR(I,J,M)
      IMPLICIT DOUBLE PRECISION (a-h,o,p,r-z)
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
      DIMENSION RM(nt,kx),a(nt,kx)

      r = 0
      last = min0(max(i),max(j))
      jp = 1
      do 10 jj =1,last
        a(jj,jp) = p(jj,i)*p(jj,j)
  10  continue
      do 20 jp = 2,ks
	do 20 jj = 1,last-jp+1
	  a(jj,jp) = p(jj,i)*p(jj+jp-1,j) + p(jj,j)*p(jj+jp-1,i)
  20  continue

      IF ( m .eq. 1 ) THEN
        do 100 jj = 2,last
	  r = r + a(jj,1)*r1(jj,ks)
  100   continue
	do 110 jp = 2,ks
*         .. j' < jj
	  do 110 jj = jp+1,last
	    r = r + a(jj-jp+1,jp)*r1(jj,ks-jp+1)
  110     continue

      ELSE IF ( m .eq. 0 ) THEN
        do 120 jj = 2,last
	  r = r + a(jj,1)*sb(jj,ks)
  120   continue
	do 130 jp = 2,ks
*         .. j' < jj
	  do 130 jj = jp+1,last
	    r = r + a(jj-jp+1,jp)*sb(jj,ks-jp+1)
  130     continue

      ELSE IF ( m .eq.-1 ) THEN
        do 200 jj = 2,last
	  r = r + a(jj,1)*rm1(jj,ks)
  200   continue
	do 210 jp = 2,ks
*         .. j' < j
	  do 210 jj = jp+1,last
	    r = r + a(jj-jp+1,jp)*rm1(jj,ks-jp+1)
  210     continue

      ELSE IF ( m .eq.-2 ) THEN
        do 300 jj = 2,last
	  r = r + a(jj,1)*rm2(jj,ks)
  300   continue
	do 310 jp = 2,ks
*         .. j' <jj 
	  do 310 jj = jp+1,last
	    r = r + a(jj-jp+1,jp)*rm2(jj,ks-jp+1)
  310     continue

      ELSE 
	nint = ns -ks + 1
	IF ( m .gt. 0) THEN
 	  call mrm(m,nt,kx,ks,nint,gr,grw,bsp,rm)
	ELSE
	  call mrm(-m,nt,kx,ks,nint,grm,grw,bsp,rm)
	END IF
        do 400 jj = 2,last
	  r = r + a(jj,1)*rm(jj,ks)
  400   continue
	do 410 jp = 2,ks
*         .. j' < jj
	  do 410 jj = jp+1,last
	    r = r + a(jj-jp+1,jp)*rm(jj,ks-jp+1)
  410     continue
      END IF
      quadr = r
      END
