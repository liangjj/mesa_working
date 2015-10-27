*     ------------------------------------------------------------------
*    		I N I T S P
*     ------------------------------------------------------------------
*
*     Initialize arrays needed for spline basis calculations
*
      SUBROUTINE initsp
      IMPLICIT DOUBLE PRECISION (a-h,o,p,r-z)
      PARAMETER (nt=100,kx=8,nwd=10)
*
      COMMON /ORBITS/p(nt,nwd), n(nwd), l(nwd), max(nwd)
      COMMON /ORBSCF/sum(nwd), s(nwd), dpm(nwd), acc(nwd), meth(nwd),
     :               iord(nwd), e(nwd,nwd)
*
      COMMON /PARAM/ z,h,t(nt),fine,mass,ks,ns,nwf,nclosd,ncfg,ib,ic,
     :               id,nscf,no,kmax,lmax,noff
*
      COMMON /SPLINE/gr(nt,kx),grm(nt,kx),grw(nt,kx),
     :        bsp(nt,kx,kx),bspd(nt,kx,kx,2),db1(nt,kx),db2(nt,kx),
     :        sb(nt,kx),r1(nt,kx),rm1(nt,kx),rm2(nt,kx),
     :        h0(nt,kx), rkm(nt,nt,kx,kx), lh, krk
      LOGICAL ex
      dimension xg(kx),wg(kx),dbiatx(nt,kx,kx),w31(nt,kx,kx),
     :          w21(nt,kx),w22(nt,kx),w1(nt,kx),w2(nt,kx),ipvtd(nt)
*
*  ... generate the grid
*
      write(0,*) ' Enter order of spline and h'
      read(5,*) ks,h
      write(0,*) ' Enter rmax, hmax'
      read(5,*) rmax,hmax
      CALL MKGRID(z,h,hmax,rmax,ks,nt,ns,t,ier)
      if (ier .eq. 2) then
	write(iscw,*) ' NT too small even for linear grid'
	stop 
      end if
      nint=ns-ks+1
      nit = nwf - ib + 1
*
*  ... generate B-spline arrays
*
      write(0,*) 'Z,ns,h,hmax,rmax,ks,nit'
      write(0,*) Z,ns,h,hmax,rmax,ks,nit
      write(6,'(6f12.6)') (t(i), i=1,ns+ks)
      write(iscw,*)'nint =',nint, 'ns =',ns
      call gauss(ks,xg,wg)
      call initvb(nt,t,kx,ks,nint,xg,wg,dbiatx,w31,w21,w22,
     :            w1,w2,gr,grm,grw,bsp,bspd)
      call initas(nt,kx,ks,nint,gr,grm,grw,bsp,bspd,db1,db2,
     :            sb,r1,rm1,rm2)
      lh = -1
      krk = -1
      END
