*deck zkscl
      subroutine zkscl (zrr, zri, fnu, n, yr, yi, nz, rzr, rzi, ascle,
     +   tol, elim)
c***begin prologue  zkscl
c***subsidiary
c***purpose  subsidiary to zbesk
c***library   slatec
c***type      all (ckscl-a, zkscl-a)
c***author  amos, d. e., (snl)
c***description
c
c     set k functions to zero on underflow, continue recurrence
c     on scaled functions until two members come on scale, then
c     return with min(nz+2,n) values scaled by 1/tol.
c
c***see also  zbesk
c***routines called  zabs, zlog, zuchk
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c   930122  added zlog to external statement.  (rwc)
c***end prologue  zkscl
c     complex ck,cs,cy,czero,rz,s1,s2,y,zr,zd,celm
      double precision acs, as, ascle, cki, ckr, csi, csr, cyi,
     * cyr, elim, fn, fnu, rzi, rzr, str, s1i, s1r, s2i,
     * s2r, tol, yi, yr, zeroi, zeror, zri, zrr, zabs,
     * zdr, zdi, celmr, elm, helim, alas
      integer i, ic, idum, kk, n, nn, nw, nz
      dimension yr(n), yi(n), cyr(2), cyi(2)
      external zabs, zlog
      data zeror,zeroi / 0.0d0 , 0.0d0 /
c***first executable statement  zkscl
      nz = 0
      ic = 0
      nn = min(2,n)
      do 10 i=1,nn
        s1r = yr(i)
        s1i = yi(i)
        cyr(i) = s1r
        cyi(i) = s1i
        as = zabs(s1r,s1i)
        acs = -zrr + log(as)
        nz = nz + 1
        yr(i) = zeror
        yi(i) = zeroi
        if (acs.lt.(-elim)) go to 10
        call zlog(s1r, s1i, csr, csi, idum)
        csr = csr - zrr
        csi = csi - zri
        str = exp(csr)/tol
        csr = str*cos(csi)
        csi = str*sin(csi)
        call zuchk(csr, csi, nw, ascle, tol)
        if (nw.ne.0) go to 10
        yr(i) = csr
        yi(i) = csi
        ic = i
        nz = nz - 1
   10 continue
      if (n.eq.1) return
      if (ic.gt.1) go to 20
      yr(1) = zeror
      yi(1) = zeroi
      nz = 2
   20 continue
      if (n.eq.2) return
      if (nz.eq.0) return
      fn = fnu + 1.0d0
      ckr = fn*rzr
      cki = fn*rzi
      s1r = cyr(1)
      s1i = cyi(1)
      s2r = cyr(2)
      s2i = cyi(2)
      helim = 0.5d0*elim
      elm = exp(-elim)
      celmr = elm
      zdr = zrr
      zdi = zri
c
c     find two consecutive y values on scale. scale recurrence if
c     s2 gets larger than exp(elim/2)
c
      do 30 i=3,n
        kk = i
        csr = s2r
        csi = s2i
        s2r = ckr*csr - cki*csi + s1r
        s2i = cki*csr + ckr*csi + s1i
        s1r = csr
        s1i = csi
        ckr = ckr + rzr
        cki = cki + rzi
        as = zabs(s2r,s2i)
        alas = log(as)
        acs = -zdr + alas
        nz = nz + 1
        yr(i) = zeror
        yi(i) = zeroi
        if (acs.lt.(-elim)) go to 25
        call zlog(s2r, s2i, csr, csi, idum)
        csr = csr - zdr
        csi = csi - zdi
        str = exp(csr)/tol
        csr = str*cos(csi)
        csi = str*sin(csi)
        call zuchk(csr, csi, nw, ascle, tol)
        if (nw.ne.0) go to 25
        yr(i) = csr
        yi(i) = csi
        nz = nz - 1
        if (ic.eq.kk-1) go to 40
        ic = kk
        go to 30
   25   continue
        if(alas.lt.helim) go to 30
        zdr = zdr - elim
        s1r = s1r*celmr
        s1i = s1i*celmr
        s2r = s2r*celmr
        s2i = s2i*celmr
   30 continue
      nz = n
      if(ic.eq.n) nz=n-1
      go to 45
   40 continue
      nz = kk - 2
   45 continue
      do 50 i=1,nz
        yr(i) = zeror
        yi(i) = zeroi
   50 continue
      return
      end
