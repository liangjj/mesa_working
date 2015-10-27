*deck cbuni
      subroutine cbuni (z, fnu, kode, n, y, nz, nui, nlast, fnul, tol,
     +   elim, alim)
c***begin prologue  cbuni
c***subsidiary
c***purpose  subsidiary to cbesi and cbesk
c***library   slatec
c***type      all (cbuni-a, zbuni-a)
c***author  amos, d. e., (snl)
c***description
c
c     cbuni computes the i bessel function for large abs(z).gt.
c     fnul and fnu+n-1.lt.fnul. the order is increased from
c     fnu+n-1 greater than fnul by adding nui and computing
c     according to the uniform asymptotic expansion for i(fnu,z)
c     on iform=1 and the expansion for j(fnu,z) on iform=2
c
c***see also  cbesi, cbesk
c***routines called  cuni1, cuni2, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  cbuni
      complex cscl, cscr, cy, rz, st, s1, s2, y, z
      real alim, ax, ay, dfnu, elim, fnu, fnui, fnul, gnu, tol, xx, yy,
     * ascle, bry, str, sti, stm, r1mach
      integer i, iflag, iform, k, kode, n, nl, nlast, nui, nw, nz
      dimension y(n), cy(2), bry(3)
c***first executable statement  cbuni
      nz = 0
      xx = real(z)
      yy = aimag(z)
      ax = abs(xx)*1.7321e0
      ay = abs(yy)
      iform = 1
      if (ay.gt.ax) iform = 2
      if (nui.eq.0) go to 60
      fnui = nui
      dfnu = fnu + (n-1)
      gnu = dfnu + fnui
      if (iform.eq.2) go to 10
c-----------------------------------------------------------------------
c     asymptotic expansion for i(fnu,z) for large fnu applied in
c     -pi/3.le.arg(z).le.pi/3
c-----------------------------------------------------------------------
      call cuni1(z, gnu, kode, 2, cy, nw, nlast, fnul, tol, elim, alim)
      go to 20
   10 continue
c-----------------------------------------------------------------------
c     asymptotic expansion for j(fnu,z*exp(m*hpi)) for large fnu
c     applied in pi/3.lt.abs(arg(z)).le.pi/2 where m=+i or -i
c     and hpi=pi/2
c-----------------------------------------------------------------------
      call cuni2(z, gnu, kode, 2, cy, nw, nlast, fnul, tol, elim, alim)
   20 continue
      if (nw.lt.0) go to 50
      if (nw.ne.0) go to 90
      ay = abs(cy(1))
c----------------------------------------------------------------------
c     scale backward recurrence, bry(3) is defined but never used
c----------------------------------------------------------------------
      bry(1) = 1.0e+3*r1mach(1)/tol
      bry(2) = 1.0e0/bry(1)
      bry(3) = bry(2)
      iflag = 2
      ascle = bry(2)
      ax = 1.0e0
      cscl = cmplx(ax,0.0e0)
      if (ay.gt.bry(1)) go to 21
      iflag = 1
      ascle = bry(1)
      ax = 1.0e0/tol
      cscl = cmplx(ax,0.0e0)
      go to 25
   21 continue
      if (ay.lt.bry(2)) go to 25
      iflag = 3
      ascle = bry(3)
      ax = tol
      cscl = cmplx(ax,0.0e0)
   25 continue
      ay = 1.0e0/ax
      cscr = cmplx(ay,0.0e0)
      s1 = cy(2)*cscl
      s2 = cy(1)*cscl
      rz = cmplx(2.0e0,0.0e0)/z
      do 30 i=1,nui
        st = s2
        s2 = cmplx(dfnu+fnui,0.0e0)*rz*s2 + s1
        s1 = st
        fnui = fnui - 1.0e0
        if (iflag.ge.3) go to 30
        st = s2*cscr
        str = real(st)
        sti = aimag(st)
        str = abs(str)
        sti = abs(sti)
        stm = max(str,sti)
        if (stm.le.ascle) go to 30
        iflag = iflag+1
        ascle = bry(iflag)
        s1 = s1*cscr
        s2 = st
        ax = ax*tol
        ay = 1.0e0/ax
        cscl = cmplx(ax,0.0e0)
        cscr = cmplx(ay,0.0e0)
        s1 = s1*cscl
        s2 = s2*cscl
   30 continue
      y(n) = s2*cscr
      if (n.eq.1) return
      nl = n - 1
      fnui = nl
      k = nl
      do 40 i=1,nl
        st = s2
        s2 = cmplx(fnu+fnui,0.0e0)*rz*s2 + s1
        s1 = st
        st = s2*cscr
        y(k) = st
        fnui = fnui - 1.0e0
        k = k - 1
        if (iflag.ge.3) go to 40
        str = real(st)
        sti = aimag(st)
        str = abs(str)
        sti = abs(sti)
        stm = max(str,sti)
        if (stm.le.ascle) go to 40
        iflag = iflag+1
        ascle = bry(iflag)
        s1 = s1*cscr
        s2 = st
        ax = ax*tol
        ay = 1.0e0/ax
        cscl = cmplx(ax,0.0e0)
        cscr = cmplx(ay,0.0e0)
        s1 = s1*cscl
        s2 = s2*cscl
   40 continue
      return
   50 continue
      nz = -1
      if(nw.eq.(-2)) nz=-2
      return
   60 continue
      if (iform.eq.2) go to 70
c-----------------------------------------------------------------------
c     asymptotic expansion for i(fnu,z) for large fnu applied in
c     -pi/3.le.arg(z).le.pi/3
c-----------------------------------------------------------------------
      call cuni1(z, fnu, kode, n, y, nw, nlast, fnul, tol, elim, alim)
      go to 80
   70 continue
c-----------------------------------------------------------------------
c     asymptotic expansion for j(fnu,z*exp(m*hpi)) for large fnu
c     applied in pi/3.lt.abs(arg(z)).le.pi/2 where m=+i or -i
c     and hpi=pi/2
c-----------------------------------------------------------------------
      call cuni2(z, fnu, kode, n, y, nw, nlast, fnul, tol, elim, alim)
   80 continue
      if (nw.lt.0) go to 50
      nz = nw
      return
   90 continue
      nlast = n
      return
      end
