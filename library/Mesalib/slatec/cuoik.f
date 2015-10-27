*deck cuoik
      subroutine cuoik (z, fnu, kode, ikflg, n, y, nuf, tol, elim, alim)
c***begin prologue  cuoik
c***subsidiary
c***purpose  subsidiary to cbesh, cbesi and cbesk
c***library   slatec
c***type      all (cuoik-a, zuoik-a)
c***author  amos, d. e., (snl)
c***description
c
c     cuoik computes the leading terms of the uniform asymptotic
c     expansions for the i and k functions and compares them
c     (in logarithmic form) to alim and elim for over and underflow
c     where alim.lt.elim. if the magnitude, based on the leading
c     exponential, is less than alim or greater than -alim, then
c     the result is on scale. if not, then a refined test using other
c     multipliers (in logarithmic form) is made based on elim. here
c     exp(-elim)=smallest machine number*1.0e+3 and exp(-alim)=
c     exp(-elim)/tol
c
c     ikflg=1 means the i sequence is tested
c          =2 means the k sequence is tested
c     nuf = 0 means the last member of the sequence is on scale
c         =-1 means an overflow would occur
c     ikflg=1 and nuf.gt.0 means the last nuf y values were set to zero
c             the first n-nuf values must be set by another routine
c     ikflg=2 and nuf.eq.n means all y values were set to zero
c     ikflg=2 and 0.lt.nuf.lt.n not considered. y must be set by
c             another routine
c
c***see also  cbesh, cbesi, cbesk
c***routines called  cuchk, cunhj, cunik, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cuoik
      complex arg, asum, bsum, cwrk, cz, czero, phi, sum, y, z, zb,
     * zeta1, zeta2, zn, zr
      real aarg, aic, alim, aphi, ascle, ax, ay, elim, fnn, fnu, gnn,
     * gnu, rcz, tol, x, yy, r1mach
      integer i, iform, ikflg, init, kode, n, nn, nuf, nw
      dimension y(n), cwrk(16)
      data czero / (0.0e0,0.0e0) /
      data aic / 1.265512123484645396e+00 /
c***first executable statement  cuoik
      nuf = 0
      nn = n
      x = real(z)
      zr = z
      if (x.lt.0.0e0) zr = -z
      zb = zr
      yy = aimag(zr)
      ax = abs(x)*1.7321e0
      ay = abs(yy)
      iform = 1
      if (ay.gt.ax) iform = 2
      gnu = max(fnu,1.0e0)
      if (ikflg.eq.1) go to 10
      fnn = nn
      gnn = fnu + fnn - 1.0e0
      gnu = max(gnn,fnn)
   10 continue
c-----------------------------------------------------------------------
c     only the magnitude of arg and phi are needed along with the
c     real parts of zeta1, zeta2 and zb. no attempt is made to get
c     the sign of the imaginary part correct.
c-----------------------------------------------------------------------
      if (iform.eq.2) go to 20
      init = 0
      call cunik(zr, gnu, ikflg, 1, tol, init, phi, zeta1, zeta2, sum,
     * cwrk)
      cz = -zeta1 + zeta2
      go to 40
   20 continue
      zn = -zr*cmplx(0.0e0,1.0e0)
      if (yy.gt.0.0e0) go to 30
      zn = conjg(-zn)
   30 continue
      call cunhj(zn, gnu, 1, tol, phi, arg, zeta1, zeta2, asum, bsum)
      cz = -zeta1 + zeta2
      aarg = abs(arg)
   40 continue
      if (kode.eq.2) cz = cz - zb
      if (ikflg.eq.2) cz = -cz
      aphi = abs(phi)
      rcz = real(cz)
c-----------------------------------------------------------------------
c     overflow test
c-----------------------------------------------------------------------
      if (rcz.gt.elim) go to 170
      if (rcz.lt.alim) go to 50
      rcz = rcz + alog(aphi)
      if (iform.eq.2) rcz = rcz - 0.25e0*alog(aarg) - aic
      if (rcz.gt.elim) go to 170
      go to 100
   50 continue
c-----------------------------------------------------------------------
c     underflow test
c-----------------------------------------------------------------------
      if (rcz.lt.(-elim)) go to 60
      if (rcz.gt.(-alim)) go to 100
      rcz = rcz + alog(aphi)
      if (iform.eq.2) rcz = rcz - 0.25e0*alog(aarg) - aic
      if (rcz.gt.(-elim)) go to 80
   60 continue
      do 70 i=1,nn
        y(i) = czero
   70 continue
      nuf = nn
      return
   80 continue
      ascle = 1.0e+3*r1mach(1)/tol
      cz = cz + clog(phi)
      if (iform.eq.1) go to 90
      cz = cz - cmplx(0.25e0,0.0e0)*clog(arg) - cmplx(aic,0.0e0)
   90 continue
      ax = exp(rcz)/tol
      ay = aimag(cz)
      cz = cmplx(ax,0.0e0)*cmplx(cos(ay),sin(ay))
      call cuchk(cz, nw, ascle, tol)
      if (nw.eq.1) go to 60
  100 continue
      if (ikflg.eq.2) return
      if (n.eq.1) return
c-----------------------------------------------------------------------
c     set underflows on i sequence
c-----------------------------------------------------------------------
  110 continue
      gnu = fnu + (nn-1)
      if (iform.eq.2) go to 120
      init = 0
      call cunik(zr, gnu, ikflg, 1, tol, init, phi, zeta1, zeta2, sum,
     * cwrk)
      cz = -zeta1 + zeta2
      go to 130
  120 continue
      call cunhj(zn, gnu, 1, tol, phi, arg, zeta1, zeta2, asum, bsum)
      cz = -zeta1 + zeta2
      aarg = abs(arg)
  130 continue
      if (kode.eq.2) cz = cz - zb
      aphi = abs(phi)
      rcz = real(cz)
      if (rcz.lt.(-elim)) go to 140
      if (rcz.gt.(-alim)) return
      rcz = rcz + alog(aphi)
      if (iform.eq.2) rcz = rcz - 0.25e0*alog(aarg) - aic
      if (rcz.gt.(-elim)) go to 150
  140 continue
      y(nn) = czero
      nn = nn - 1
      nuf = nuf + 1
      if (nn.eq.0) return
      go to 110
  150 continue
      ascle = 1.0e+3*r1mach(1)/tol
      cz = cz + clog(phi)
      if (iform.eq.1) go to 160
      cz = cz - cmplx(0.25e0,0.0e0)*clog(arg) - cmplx(aic,0.0e0)
  160 continue
      ax = exp(rcz)/tol
      ay = aimag(cz)
      cz = cmplx(ax,0.0e0)*cmplx(cos(ay),sin(ay))
      call cuchk(cz, nw, ascle, tol)
      if (nw.eq.1) go to 140
      return
  170 continue
      nuf = -1
      return
      end
