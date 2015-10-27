*     ------------------------------------------------------------------
*              R K Y
*     ------------------------------------------------------------------
*
*                   k
*       Evaluates  R (i1, j1; i2, j2)
*
*     The COMMON /SLATER/ is used to store results from previous
*     calls so that, when many integrals are evaluated, much of the
*     calculation can be avoided.  In this version, it is assumed
*     that the Yk function has been computed for j1,j2 and integration
*     occurs with i1,i2
*
      DOUBLE PRECISION FUNCTION RKy(I1,J1,I2,J2,K,REL)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      PARAMETER (nt=100,kx=8,ktx=3*kx-2,nwd=10)
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rkm(nt,nt,kx,kx),lh,krk

      COMMON /SLATER/fc1(nt,kx),fc2(nt,kx),yk(nt),fyk(nt,kx),
     :        dyk(ktx,nt),ipvtd(nt),fc(nt,kx),ic1,ic2,iy1,iy2,ky
*
      LOGICAL same, REL
*    
*     .. test if YK needs to be computed.
      print *, ' RKY: k,i1,j1,i2,j2',k,i1,j1,i2,j2
      print *, 'ky,iy1,iy2',ky,iy1,iy2, 'ic1,ic2', ic1,ic2
      if (k .ne. ky .or. (j1 .ne. iy1) .or. (j2 .ne. iy2)) then
	ky =k
	iy1 = j1
	iy2 = j2
        nint = ns-ks+1
        call ykf(j1,j2,k,rel)
        ic1 = 0
        ic2 = 0
      end if
      if (i1 .ne. ic1) then
	call yval(0,0,nt,kx,ks,nint,gr,bsp,p(1,i1),fc1)
	ic1 = i1
      end if
      if (i2 .ne. ic2) then
	call yval(0,0,nt,kx,ks,nint,gr,bsp,p(1,i2),fc2)
	ic2 = i2
      end if
      rky = 0.d0
      do 100 m = 1,ks
	 do 110 ii = 1,nint
	   fc(ii,m) = fc1(ii,m)*fc2(ii,m)*fyk(ii,m)
  110    continue
	 rky = rky + ddot(nint,fc(1,m),1,grw(1,m),1)
  100 continue
      end
