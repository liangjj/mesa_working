*     ==================================================================
*              R E L - S H I F T
*     ==================================================================
*
*
*       Computes the mass velocity  and one-body Darwin term
*   corrections for the relativistic shift in the energy of the electron
*   including non-diagonal corrections
*
*
      DOUBLE PRECISION  FUNCTION RLSHFT(I1,I2)
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
      DIMENSION vc(nt,kx),a(nt,kx)
*
* *****  form the mass-velocity correction matrix in the spline basis
*
      LL = L(i1)
      print *, i1,i2,ll
      call mvc(LL,nt,kx,ks,ns,grm,grw,bsp,bspd,vc)
*
* .... store the product of coefficients
*
      last = min0(max(i1),max(i2))
      print *, last
      do 10 j =1,last
        a(j,1) = p(j,i1)*p(j,i2)
  10  continue
      do 20 jp = 2,ks
	do 20 j = 1,last-jp+1
	  a(j,jp) = p(j,i1)*p(j+jp-1,i2) + p(j,i2)*p(j+jp-1,i1)
  20  continue
*
* ...  form the sums
*
	r = 0.d0
        do 100 j = 2,last
	  r = r + a(j,1)*vc(j,ks)
  100   continue
	do 110 jp = 2,ks
*         .. j' < j
	  do 110 j = jp+1,last
	    r = r + a(j-jp+1,jp)*vc(j,ks-jp+1)
  110     continue

      print *, 'RLSHFT: relsh,z,h,ks,last',r,z,h,ks,last,fine
      IF ( L(I1) .EQ. 0 ) then
	RELSH = -r + Z*AZ(z,h,ks,LL+1,p(2,I1))*AZ(z,h,ks,LL+1,p(2,I2))
      else
	RELSH = -r
      end if
      print *, 'RLSHFT:', relsh
      RLSHFT = 0.5d0*RELSH*FINE
      RETURN
      END
