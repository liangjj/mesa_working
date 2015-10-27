*     ==================================================================
*              G R A D
*     ==================================================================
*
*  *****  THE GRAD FUNCTION SUBPROGRAM COMPUTES THE FOLLOWING DIRECTLY
*  *****         <P(J)^D + L(I)/R ^P(I)> WITH L(I) > L(J)
*
      DOUBLE PRECISION  FUNCTION GRAD(I,J,ier)
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
      LOGICAL ixch
      DIMENSION g(nt,kx),a(nt,kx)
      
      IF ( IABS(L(I) - L(J)) .NE. 1) THEN
	ier = 1
	grad = 0.0d0
      else
      LL = MAX0(L(I),L(J))
      i1 = i 
      i2 = j
      ixch = .false.
      IF ( L(I) .lt. L(J) ) THEN
	i1 = j
	i2 = i
	ixch = .true.
      end if
*
* ...  form the sums
*
      last = min0(max(i1),max(i2))
      grad = 0.d0
      do 100 jj = 2,last
         grad = grad + p(jj,i1)*p(jj,i2)*rm1(jj,ks)
  100 continue
      grad = ll*grad
      do 110 jp = 2,ks
*         .. j' < jj
        do 110 jj = 2,last-jp+1
          grad = grad +
     :     g(jj,ks-jp+1)*(p(jj,i2)*p(jj-jp+1,i1)-
     :                    p(jj,i1)*p(jj-jp+1,i2)) +
     :     ll*rm1(jj,ks-jp+1)*(p(jj,i2)*p(jj-jp+1,i1)+
     :                         p(jj,i1)*p(jj-jp+1,i2)) 
  110     continue

      IF ( ixch ) GRAD = - GRAD
      end if
      END
