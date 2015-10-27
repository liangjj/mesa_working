*deck zuoik
      subroutine zuoik (zr, zi, fnu, kode, ikflg, n, yr, yi, nuf, tol,
     +   elim, alim)
c***begin prologue  zuoik
c***subsidiary
c***purpose  subsidiary to zbesh, zbesi and zbesk
c***library   slatec
c***type      all (cuoik-a, zuoik-a)
c***author  amos, d. e., (snl)
c***description
c
c     zuoik computes the leading terms of the uniform asymptotic
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
c***see also  zbesh, zbesi, zbesk
c***routines called  d1mach, zabs, zlog, zuchk, zunhj, zunik
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c   930122  added zlog to external statement.  (rwc)
c***end prologue  zuoik
c     complex arg,asum,bsum,cwrk,cz,czero,phi,sum,y,z,zb,zeta1,zeta2,zn,
c    *zr
      double precision aarg, aic, alim, aphi, argi, argr, asumi, asumr,
     * ascle, ax, ay, bsumi, bsumr, cwrki, cwrkr, czi, czr, elim, fnn,
     * fnu, gnn, gnu, phii, phir, rcz, str, sti, sumi, sumr, tol, yi,
     * yr, zbi, zbr, zeroi, zeror, zeta1i, zeta1r, zeta2i, zeta2r, zi,
     * zni, znr, zr, zri, zrr, d1mach, zabs
      integer i, idum, iform, ikflg, init, kode, n, nn, nuf, nw
      dimension yr(n), yi(n), cwrkr(16), cwrki(16)
      external zabs, zlog
      data zeror,zeroi / 0.0d0, 0.0d0 /
      data aic / 1.265512123484645396d+00 /
c***first executable statement  zuoik
      nuf = 0
      nn = n
      zrr = zr
      zri = zi
      if (zr.ge.0.0d0) go to 10
      zrr = -zr
      zri = -zi
   10 continue
      zbr = zrr
      zbi = zri
      ax = abs(zr)*1.7321d0
      ay = abs(zi)
      iform = 1
      if (ay.gt.ax) iform = 2
      gnu = max(fnu,1.0d0)
      if (ikflg.eq.1) go to 20
      fnn = nn
      gnn = fnu + fnn - 1.0d0
      gnu = max(gnn,fnn)
   20 continue
c-----------------------------------------------------------------------
c     only the magnitude of arg and phi are needed along with the
c     real parts of zeta1, zeta2 and zb. no attempt is made to get
c     the sign of the imaginary part correct.
c-----------------------------------------------------------------------
      if (iform.eq.2) go to 30
      init = 0
      call zunik(zrr, zri, gnu, ikflg, 1, tol, init, phir, phii,
     * zeta1r, zeta1i, zeta2r, zeta2i, sumr, sumi, cwrkr, cwrki)
      czr = -zeta1r + zeta2r
      czi = -zeta1i + zeta2i
      go to 50
   30 continue
      znr = zri
      zni = -zrr
      if (zi.gt.0.0d0) go to 40
      znr = -znr
   40 continue
      call zunhj(znr, zni, gnu, 1, tol, phir, phii, argr, argi, zeta1r,
     * zeta1i, zeta2r, zeta2i, asumr, asumi, bsumr, bsumi)
      czr = -zeta1r + zeta2r
      czi = -zeta1i + zeta2i
      aarg = zabs(argr,argi)
   50 continue
      if (kode.eq.1) go to 60
      czr = czr - zbr
      czi = czi - zbi
   60 continue
      if (ikflg.eq.1) go to 70
      czr = -czr
      czi = -czi
   70 continue
      aphi = zabs(phir,phii)
      rcz = czr
c-----------------------------------------------------------------------
c     overflow test
c-----------------------------------------------------------------------
      if (rcz.gt.elim) go to 210
      if (rcz.lt.alim) go to 80
      rcz = rcz + log(aphi)
      if (iform.eq.2) rcz = rcz - 0.25d0*log(aarg) - aic
      if (rcz.gt.elim) go to 210
      go to 130
   80 continue
c-----------------------------------------------------------------------
c     underflow test
c-----------------------------------------------------------------------
      if (rcz.lt.(-elim)) go to 90
      if (rcz.gt.(-alim)) go to 130
      rcz = rcz + log(aphi)
      if (iform.eq.2) rcz = rcz - 0.25d0*log(aarg) - aic
      if (rcz.gt.(-elim)) go to 110
   90 continue
      do 100 i=1,nn
        yr(i) = zeror
        yi(i) = zeroi
  100 continue
      nuf = nn
      return
  110 continue
      ascle = 1.0d+3*d1mach(1)/tol
      call zlog(phir, phii, str, sti, idum)
      czr = czr + str
      czi = czi + sti
      if (iform.eq.1) go to 120
      call zlog(argr, argi, str, sti, idum)
      czr = czr - 0.25d0*str - aic
      czi = czi - 0.25d0*sti
  120 continue
      ax = exp(rcz)/tol
      ay = czi
      czr = ax*cos(ay)
      czi = ax*sin(ay)
      call zuchk(czr, czi, nw, ascle, tol)
      if (nw.ne.0) go to 90
  130 continue
      if (ikflg.eq.2) return
      if (n.eq.1) return
c-----------------------------------------------------------------------
c     set underflows on i sequence
c-----------------------------------------------------------------------
  140 continue
      gnu = fnu + (nn-1)
      if (iform.eq.2) go to 150
      init = 0
      call zunik(zrr, zri, gnu, ikflg, 1, tol, init, phir, phii,
     * zeta1r, zeta1i, zeta2r, zeta2i, sumr, sumi, cwrkr, cwrki)
      czr = -zeta1r + zeta2r
      czi = -zeta1i + zeta2i
      go to 160
  150 continue
      call zunhj(znr, zni, gnu, 1, tol, phir, phii, argr, argi, zeta1r,
     * zeta1i, zeta2r, zeta2i, asumr, asumi, bsumr, bsumi)
      czr = -zeta1r + zeta2r
      czi = -zeta1i + zeta2i
      aarg = zabs(argr,argi)
  160 continue
      if (kode.eq.1) go to 170
      czr = czr - zbr
      czi = czi - zbi
  170 continue
      aphi = zabs(phir,phii)
      rcz = czr
      if (rcz.lt.(-elim)) go to 180
      if (rcz.gt.(-alim)) return
      rcz = rcz + log(aphi)
      if (iform.eq.2) rcz = rcz - 0.25d0*log(aarg) - aic
      if (rcz.gt.(-elim)) go to 190
  180 continue
      yr(nn) = zeror
      yi(nn) = zeroi
      nn = nn - 1
      nuf = nuf + 1
      if (nn.eq.0) return
      go to 140
  190 continue
      ascle = 1.0d+3*d1mach(1)/tol
      call zlog(phir, phii, str, sti, idum)
      czr = czr + str
      czi = czi + sti
      if (iform.eq.1) go to 200
      call zlog(argr, argi, str, sti, idum)
      czr = czr - 0.25d0*str - aic
      czi = czi - 0.25d0*sti
  200 continue
      ax = exp(rcz)/tol
      ay = czi
      czr = ax*cos(ay)
      czi = ax*sin(ay)
      call zuchk(czr, czi, nw, ascle, tol)
      if (nw.ne.0) go to 180
      return
  210 continue
      nuf = -1
      return
      end
