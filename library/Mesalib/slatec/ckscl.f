*deck ckscl
      subroutine ckscl (zr, fnu, n, y, nz, rz, ascle, tol, elim)
c***begin prologue  ckscl
c***subsidiary
c***purpose  subsidiary to cbknu, cunk1 and cunk2
c***library   slatec
c***type      all (ckscl-a, zkscl-a)
c***author  amos, d. e., (snl)
c***description
c
c     set k functions to zero on underflow, continue recurrence
c     on scaled functions until two members come on scale, then
c     return with min(nz+2,n) values scaled by 1/tol.
c
c***see also  cbknu, cunk1, cunk2
c***routines called  cuchk
c***revision history  (yymmdd)
c   ??????  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  ckscl
      complex ck, cs, cy, czero, rz, s1, s2, y, zr, zd, celm
      real aa, ascle, acs, as, csi, csr, elim, fn, fnu, tol, xx, zri,
     * elm, alas, helim
      integer i, ic, k, kk, n, nn, nw, nz
      dimension y(n), cy(2)
      data czero / (0.0e0,0.0e0) /
c***first executable statement  cuchk
      nz = 0
      ic = 0
      xx = real(zr)
      nn = min(2,n)
      do 10 i=1,nn
        s1 = y(i)
        cy(i) = s1
        as = abs(s1)
        acs = -xx + alog(as)
        nz = nz + 1
        y(i) = czero
        if (acs.lt.(-elim)) go to 10
        cs = -zr + clog(s1)
        csr = real(cs)
        csi = aimag(cs)
        aa = exp(csr)/tol
        cs = cmplx(aa,0.0e0)*cmplx(cos(csi),sin(csi))
        call cuchk(cs, nw, ascle, tol)
        if (nw.ne.0) go to 10
        y(i) = cs
        nz = nz - 1
        ic = i
   10 continue
      if (n.eq.1) return
      if (ic.gt.1) go to 20
      y(1) = czero
      nz = 2
   20 continue
      if (n.eq.2) return
      if (nz.eq.0) return
      fn = fnu + 1.0e0
      ck = cmplx(fn,0.0e0)*rz
      s1 = cy(1)
      s2 = cy(2)
      helim = 0.5e0*elim
      elm = exp(-elim)
      celm = cmplx(elm,0.0e0)
      zri =aimag(zr)
      zd = zr
c
c     find two consecutive y values on scale. scale recurrence if
c     s2 gets larger than exp(elim/2)
c
      do 30 i=3,n
        kk = i
        cs = s2
        s2 = ck*s2 + s1
        s1 = cs
        ck = ck + rz
        as = abs(s2)
        alas = alog(as)
        acs = -xx + alas
        nz = nz + 1
        y(i) = czero
        if (acs.lt.(-elim)) go to 25
        cs = -zd + clog(s2)
        csr = real(cs)
        csi = aimag(cs)
        aa = exp(csr)/tol
        cs = cmplx(aa,0.0e0)*cmplx(cos(csi),sin(csi))
        call cuchk(cs, nw, ascle, tol)
        if (nw.ne.0) go to 25
        y(i) = cs
        nz = nz - 1
        if (ic.eq.(kk-1)) go to 40
        ic = kk
        go to 30
   25   continue
        if(alas.lt.helim) go to 30
        xx = xx-elim
        s1 = s1*celm
        s2 = s2*celm
        zd = cmplx(xx,zri)
   30 continue
      nz = n
      if(ic.eq.n) nz=n-1
      go to 45
   40 continue
      nz = kk - 2
   45 continue
      do 50 k=1,nz
        y(k) = czero
   50 continue
      return
      end
