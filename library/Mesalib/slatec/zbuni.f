*deck zbuni
      subroutine zbuni (zr, zi, fnu, kode, n, yr, yi, nz, nui, nlast,
     +   fnul, tol, elim, alim)
c***begin prologue  zbuni
c***subsidiary
c***purpose  subsidiary to zbesi and zbesk
c***library   slatec
c***type      all (cbuni-a, zbuni-a)
c***author  amos, d. e., (snl)
c***description
c
c     zbuni computes the i bessel function for large abs(z).gt.
c     fnul and fnu+n-1.lt.fnul. the order is increased from
c     fnu+n-1 greater than fnul by adding nui and computing
c     according to the uniform asymptotic expansion for i(fnu,z)
c     on iform=1 and the expansion for j(fnu,z) on iform=2
c
c***see also  zbesi, zbesk
c***routines called  d1mach, zabs, zuni1, zuni2
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zbuni
c     complex cscl,cscr,cy,rz,st,s1,s2,y,z
      double precision alim, ax, ay, csclr, cscrr, cyi, cyr, dfnu,
     * elim, fnu, fnui, fnul, gnu, raz, rzi, rzr, sti, str, s1i, s1r,
     * s2i, s2r, tol, yi, yr, zi, zr, zabs, ascle, bry, c1r, c1i, c1m,
     * d1mach
      integer i, iflag, iform, k, kode, n, nl, nlast, nui, nw, nz
      dimension yr(n), yi(n), cyr(2), cyi(2), bry(3)
      external zabs
c***first executable statement  zbuni
      nz = 0
      ax = abs(zr)*1.7321d0
      ay = abs(zi)
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
      call zuni1(zr, zi, gnu, kode, 2, cyr, cyi, nw, nlast, fnul, tol,
     * elim, alim)
      go to 20
   10 continue
c-----------------------------------------------------------------------
c     asymptotic expansion for j(fnu,z*exp(m*hpi)) for large fnu
c     applied in pi/3.lt.abs(arg(z)).le.pi/2 where m=+i or -i
c     and hpi=pi/2
c-----------------------------------------------------------------------
      call zuni2(zr, zi, gnu, kode, 2, cyr, cyi, nw, nlast, fnul, tol,
     * elim, alim)
   20 continue
      if (nw.lt.0) go to 50
      if (nw.ne.0) go to 90
      str = zabs(cyr(1),cyi(1))
c----------------------------------------------------------------------
c     scale backward recurrence, bry(3) is defined but never used
c----------------------------------------------------------------------
      bry(1)=1.0d+3*d1mach(1)/tol
      bry(2) = 1.0d0/bry(1)
      bry(3) = bry(2)
      iflag = 2
      ascle = bry(2)
      csclr = 1.0d0
      if (str.gt.bry(1)) go to 21
      iflag = 1
      ascle = bry(1)
      csclr = 1.0d0/tol
      go to 25
   21 continue
      if (str.lt.bry(2)) go to 25
      iflag = 3
      ascle=bry(3)
      csclr = tol
   25 continue
      cscrr = 1.0d0/csclr
      s1r = cyr(2)*csclr
      s1i = cyi(2)*csclr
      s2r = cyr(1)*csclr
      s2i = cyi(1)*csclr
      raz = 1.0d0/zabs(zr,zi)
      str = zr*raz
      sti = -zi*raz
      rzr = (str+str)*raz
      rzi = (sti+sti)*raz
      do 30 i=1,nui
        str = s2r
        sti = s2i
        s2r = (dfnu+fnui)*(rzr*str-rzi*sti) + s1r
        s2i = (dfnu+fnui)*(rzr*sti+rzi*str) + s1i
        s1r = str
        s1i = sti
        fnui = fnui - 1.0d0
        if (iflag.ge.3) go to 30
        str = s2r*cscrr
        sti = s2i*cscrr
        c1r = abs(str)
        c1i = abs(sti)
        c1m = max(c1r,c1i)
        if (c1m.le.ascle) go to 30
        iflag = iflag+1
        ascle = bry(iflag)
        s1r = s1r*cscrr
        s1i = s1i*cscrr
        s2r = str
        s2i = sti
        csclr = csclr*tol
        cscrr = 1.0d0/csclr
        s1r = s1r*csclr
        s1i = s1i*csclr
        s2r = s2r*csclr
        s2i = s2i*csclr
   30 continue
      yr(n) = s2r*cscrr
      yi(n) = s2i*cscrr
      if (n.eq.1) return
      nl = n - 1
      fnui = nl
      k = nl
      do 40 i=1,nl
        str = s2r
        sti = s2i
        s2r = (fnu+fnui)*(rzr*str-rzi*sti) + s1r
        s2i = (fnu+fnui)*(rzr*sti+rzi*str) + s1i
        s1r = str
        s1i = sti
        str = s2r*cscrr
        sti = s2i*cscrr
        yr(k) = str
        yi(k) = sti
        fnui = fnui - 1.0d0
        k = k - 1
        if (iflag.ge.3) go to 40
        c1r = abs(str)
        c1i = abs(sti)
        c1m = max(c1r,c1i)
        if (c1m.le.ascle) go to 40
        iflag = iflag+1
        ascle = bry(iflag)
        s1r = s1r*cscrr
        s1i = s1i*cscrr
        s2r = str
        s2i = sti
        csclr = csclr*tol
        cscrr = 1.0d0/csclr
        s1r = s1r*csclr
        s1i = s1i*csclr
        s2r = s2r*csclr
        s2i = s2i*csclr
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
      call zuni1(zr, zi, fnu, kode, n, yr, yi, nw, nlast, fnul, tol,
     * elim, alim)
      go to 80
   70 continue
c-----------------------------------------------------------------------
c     asymptotic expansion for j(fnu,z*exp(m*hpi)) for large fnu
c     applied in pi/3.lt.abs(arg(z)).le.pi/2 where m=+i or -i
c     and hpi=pi/2
c-----------------------------------------------------------------------
      call zuni2(zr, zi, fnu, kode, n, yr, yi, nw, nlast, fnul, tol,
     * elim, alim)
   80 continue
      if (nw.lt.0) go to 50
      nz = nw
      return
   90 continue
      nlast = n
      return
      end
